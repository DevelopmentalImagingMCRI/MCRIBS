/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "mirtk/ConformalSurfaceFlattening.h"

#include "mirtk/VtkMath.h"
#include "mirtk/Triangle.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/SurfaceCurvature.h"
#include "mirtk/PiecewiseLinearMap.h"
#include "mirtk/PointSetUtils.h"

#include "vtkIdList.h"
#include "vtkPointData.h"
#include "vtkCellData.h"

#if MIRTK_USE_FLOAT_BY_DEFAULT
  #include "vtkFloatArray.h"
#else
  #include "vtkDoubleArray.h"
#endif

#include "Eigen/SparseCore"
#include "Eigen/SparseLU"
#include "Eigen/OrderingMethods"
#include "Eigen/IterativeLinearSolvers"


namespace mirtk {


// Global flags (cf. mirtk/Options.h)
MIRTK_Common_EXPORT extern int verbose;


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void ConformalSurfaceFlattening
::CopyAttributes(const ConformalSurfaceFlattening &other)
{
  _PolarCellId        = other._PolarCellId;
  _MapToSphere        = other._MapToSphere;
  _Scale              = other._Scale;
  _Radius             = other._Radius;
  _NumberOfIterations = other._NumberOfIterations;
  _Tolerance          = other._Tolerance;
}

// -----------------------------------------------------------------------------
ConformalSurfaceFlattening::ConformalSurfaceFlattening()
:
  _PolarCellId(-1),
  _MapToSphere(true),
  _Scale(0.),
  _Radius(1.),
  _NumberOfIterations(-1),
  _Tolerance(-1.)
{
}

// -----------------------------------------------------------------------------
ConformalSurfaceFlattening::ConformalSurfaceFlattening(
  const ConformalSurfaceFlattening &other
) :
  SphericalSurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
ConformalSurfaceFlattening &ConformalSurfaceFlattening
::operator =(const ConformalSurfaceFlattening &other)
{
  if (this != &other) {
    SphericalSurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
ConformalSurfaceFlattening::~ConformalSurfaceFlattening()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void ConformalSurfaceFlattening::Initialize()
{
  // Initialize base class
  SphericalSurfaceMapper::Initialize();

  // Check that fixed point cell ID is valid
  if (_PolarCellId >= _Surface->GetNumberOfCells()) {
    cerr << this->NameOfType() << "::Initialize: Invalid fixed point cell ID!" << endl;
    exit(1);
  }

  // Select cell with lowest average curvature as cell containing the polar point
  if (_PolarCellId < 0) {
    SurfaceCurvature filter(SurfaceCurvature::Mean);
    filter.Input(_Surface);
    filter.Run();
    vtkDataArray *curvature = filter.GetMeanCurvature();
    for (vtkIdType ptId = 1; ptId < curvature->GetNumberOfTuples(); ++ptId) {
      curvature->SetComponent(ptId, 0, abs(curvature->GetComponent(ptId, 0)));
    }
    vtkIdType polarPoint = 0;
    double    polarValue = curvature->GetComponent(0, 0);
    for (vtkIdType ptId = 1; ptId < curvature->GetNumberOfTuples(); ++ptId) {
      if (curvature->GetComponent(ptId, 0) < polarValue) {
        polarValue = curvature->GetComponent(ptId, 0);
        polarPoint = ptId;
      }
    }
    unsigned short ncells;
    vtkIdType      npts, *cells, *pts;
    _Surface->GetPointCells(polarPoint, ncells, cells);
    _PolarCellId = cells[0];
    double avgValue, minValue = .0;
    for (unsigned short i = 0; i < ncells; ++i) {
      _Surface->GetCellPoints(cells[i], npts, pts);
      avgValue = 0.;
      for (vtkIdType j = 0; j < npts; ++j) {
        avgValue += curvature->GetComponent(pts[j], 0);
      }
      avgValue /= npts;
      if (i == 0 || avgValue < minValue) {
        _PolarCellId = cells[i];
        minValue     = avgValue;
      }
    }
  }

  // Initialize map values
  #if MIRTK_USE_FLOAT_BY_DEFAULT
    _Values = vtkSmartPointer<vtkFloatArray>::New();
  #else
    _Values = vtkSmartPointer<vtkDoubleArray>::New();
  #endif
  _Values->SetName("SurfaceMap");
  _Values->SetNumberOfComponents(_MapToSphere ? 3 : 2);
  _Values->SetNumberOfTuples(_Surface->GetNumberOfPoints());
  for (int j = 0; j < _Values->GetNumberOfComponents(); ++j) {
    _Values->FillComponent(j, 0.);
  }
}

// -----------------------------------------------------------------------------
void ConformalSurfaceFlattening::ComputeMap()
{
  const double use_direct_solver = (_NumberOfIterations < 0 || _NumberOfIterations == 1);

  MIRTK_START_TIMING();

  typedef Eigen::MatrixXd             Values;
  typedef Eigen::SparseMatrix<double> Matrix;

  const int n = NumberOfPoints();
  const int m = 2;

  // Calculate matrix D
  Matrix D(n, n);
  {
    typedef Eigen::VectorXi Vector;

    vtkIdType npts, *pts;                // cell links list reference
    double posA [3], posB [3], posC [3]; // position of cell corner points
    double cotABC, cotBCA, cotCAB;       // cotangent of cell angles

    D.reserve(Vector::Constant(n, _EdgeTable->MaxNumberOfAdjacentPoints() + 1));
    for (vtkIdType cellId = 0; cellId < _Surface->GetNumberOfCells(); ++cellId) {
      _Surface->GetCellPoints(cellId, npts, pts);
      if (npts != 3) {
        cerr << this->NameOfType() << "::ComputeMap: Surface mesh must be triangulated!" << endl;
        exit(1);
      }

      const vtkIdType &ptIdA = pts[0];
      const vtkIdType &ptIdB = pts[1];
      const vtkIdType &ptIdC = pts[2];

      _Surface->GetPoint(ptIdA, posA);
      _Surface->GetPoint(ptIdB, posB);
      _Surface->GetPoint(ptIdC, posC);

      cotABC = Triangle::Cotangent(posA, posB, posC);
      cotBCA = Triangle::Cotangent(posB, posC, posA);
      cotCAB = Triangle::Cotangent(posC, posA, posB);

      D.coeffRef(ptIdA, ptIdA) += cotABC + cotBCA;
      D.coeffRef(ptIdA, ptIdB) -= cotBCA;
      D.coeffRef(ptIdA, ptIdC) -= cotABC;

      D.coeffRef(ptIdB, ptIdB) += cotBCA + cotCAB;
      D.coeffRef(ptIdB, ptIdA) -= cotBCA;
      D.coeffRef(ptIdB, ptIdC) -= cotCAB;

      D.coeffRef(ptIdC, ptIdC) += cotCAB + cotABC;
      D.coeffRef(ptIdC, ptIdB) -= cotCAB;
      D.coeffRef(ptIdC, ptIdA) -= cotABC;
    }
    D.makeCompressed();
  }

  // Calculate (complex) right hand side vector b
  Values b(n, m);
  {
    vtkIdType npts, *pts;             // polar cell links list reference
    double posA[3], posB[3], posC[3]; // position of cell corner points
    double vecAB[3], vecCA[3];        // edge vectors AB and AC
    double posE[3];                   // projection of C onto AB
    double vecEC[3];                  // vector from E to C

    // Get corners of cell containing polar point
    _Surface->GetCellPoints(_PolarCellId, npts, pts);
    const vtkIdType &ptIdA = pts[0];
    const vtkIdType &ptIdB = pts[1];
    const vtkIdType &ptIdC = pts[2];

    _Surface->GetPoint(ptIdA, posA);
    _Surface->GetPoint(ptIdB, posB);
    _Surface->GetPoint(ptIdC, posC);

    // Edge vectors AB and AC
    vtkMath::Subtract(posB, posA, vecAB);
    vtkMath::Subtract(posA, posC, vecCA);

    // Orthogonal projection of C onto AB
    const double normAB2 = vtkMath::Dot(vecAB, vecAB);
    const double theta   = - vtkMath::Dot(vecCA, vecAB) / normAB2;
    posE[0] = posA[0] + theta * vecAB[0];
    posE[1] = posA[1] + theta * vecAB[1];
    posE[2] = posA[2] + theta * vecAB[2];
    vtkMath::Subtract(posC, posE, vecEC);

    // Constraints of corner points
    b.setZero();

    const double x = 2. / sqrt(normAB2);
    b(ptIdA, 0) = -x;
    b(ptIdB, 0) =  x;

    const double y = 2. / vtkMath::Norm(vecEC);
    b(ptIdA, 1) =  y * (1. - theta);
    b(ptIdB, 1) =  y * theta;
    b(ptIdC, 1) = -y;
  }

  MIRTK_DEBUG_TIMING(1, "building sparse linear system");

  // Solve linear system
  if (verbose) {
    cout << "\n";
    cout << "  No. of surface points  = " << NumberOfPoints() << "\n";
    cout << "  No. of fixed points    = 3\n";
    cout << "  No. of non-zero values = " << D.nonZeros() << "\n";
    cout << "  Dimension of codomain  = " << m << "\n";
    cout.flush();
  }

  MIRTK_RESET_TIMING();

  Values x(n, m);
  int    niter = 0;
  double error = nan;

  if (use_direct_solver) {
    Eigen::SparseLU<Matrix> solver;
    solver.analyzePattern(D);
    solver.factorize(D);
    x = solver.solve(b);
  } else {
    Eigen::ConjugateGradient<Matrix> solver(D);
    if (_NumberOfIterations > 0 ) solver.setMaxIterations(_NumberOfIterations);
    if (_Tolerance          > 0.) solver.setTolerance(_Tolerance);
    x = solver.solve(b);
    niter = solver.iterations();
    error = solver.error();
  }

  for (int i = 0; i < n; ++i) {
    for (int l = 0; l < m; ++l) {
      _Values->SetComponent(static_cast<vtkIdType>(i), l, x(i, l));
    }
  }

  MIRTK_DEBUG_TIMING(1, "solving sparse linear system");

  if (verbose && !use_direct_solver) {
    cout << "  No. of iterations      = " << niter << "\n";
    cout << "  Estimated error        = " << error << "\n";
    cout.flush();
  }
}

// -----------------------------------------------------------------------------
void ConformalSurfaceFlattening::Finalize()
{
  // Inverse stereographic projection
  if (_MapToSphere) {
    const int n = NumberOfPoints();
    double x, y, r2, scale = _Scale;

    // Choose scale for inverse stereographic projection automatically s.t.
    // afterwards, upper and lower hemisphere have same number of points
    if (scale <= .0) {
      Array<double> v_r2(n);
      auto v_r2_it = v_r2.begin();
      for (int i = 0; i < n; ++i, ++v_r2_it) {
        x = _Values->GetComponent(static_cast<vtkIdType>(i), 0);
        y = _Values->GetComponent(static_cast<vtkIdType>(i), 1);
        *v_r2_it = x*x + y*y;
      }
      sort(v_r2.begin(), v_r2.end());
      int i = ((n % 2) == 0 ? n : n - 1) / 2;
      scale = 1.0 / sqrt(v_r2[i]);
    }

    // Perform inverse stereographic projection
    const double s = 2. * _Radius;
    for (int i = 0; i < n; ++i) {
      x = scale * _Values->GetComponent(static_cast<vtkIdType>(i), 0);
      y = scale * _Values->GetComponent(static_cast<vtkIdType>(i), 1);
      r2 = x*x + y*y;
      _Values->SetComponent(static_cast<vtkIdType>(i), 0, s *  x / (1. + r2));
      _Values->SetComponent(static_cast<vtkIdType>(i), 1, s *  y / (1. + r2));
      _Values->SetComponent(static_cast<vtkIdType>(i), 2, s * r2 / (1. + r2) - 1.);
    }
  }

  // Set output surface map
  SharedPtr<PiecewiseLinearMap> map = NewShared<PiecewiseLinearMap>();
  vtkSmartPointer<vtkPolyData> domain;
  domain.TakeReference(_Surface->NewInstance());
  domain->ShallowCopy(_Surface);
  domain->GetPointData()->Initialize();
  domain->GetCellData()->Initialize();
  map->Domain(domain);
  map->Values(_Values);
  _Output = map;

  // Finalize base class
  SphericalSurfaceMapper::Finalize();
}


} // namespace mirtk
