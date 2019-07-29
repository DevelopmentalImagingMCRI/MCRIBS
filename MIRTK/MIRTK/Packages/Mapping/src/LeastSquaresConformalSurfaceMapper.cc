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

#include "mirtk/LeastSquaresConformalSurfaceMapper.h"

#include "mirtk/Algorithm.h"
#include "mirtk/Triangle.h"
#include "mirtk/PiecewiseLinearMap.h"
#include "mirtk/ChordLengthBoundarySegmentParameterizer.h"

#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"

#include "Eigen/Sparse"


namespace mirtk {


// Global flags (cf. mirtk/Options.h)
MIRTK_Common_EXPORT extern int verbose;


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void LeastSquaresConformalSurfaceMapper
::CopyAttributes(const LeastSquaresConformalSurfaceMapper &other)
{
  _NumberOfIterations = other._NumberOfIterations;
  _Tolerance          = other._Tolerance;
  _PointIndex         = other._PointIndex;
  _FreePoints         = other._FreePoints;
  _FixedPoints        = other._FixedPoints;
  _FixedValues        = other._FixedValues;

  if (other._Values) {
    _Values.TakeReference(other._Values->NewInstance());
    _Values->DeepCopy(other._Values);
  } else {
    _Values = nullptr;
  }
}

// -----------------------------------------------------------------------------
LeastSquaresConformalSurfaceMapper::LeastSquaresConformalSurfaceMapper()
:
  _NumberOfIterations(-1),
  _Tolerance(-1.)
{
}

// -----------------------------------------------------------------------------
LeastSquaresConformalSurfaceMapper
::LeastSquaresConformalSurfaceMapper(const LeastSquaresConformalSurfaceMapper &other)
:
  FreeBoundarySurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
LeastSquaresConformalSurfaceMapper &LeastSquaresConformalSurfaceMapper
::operator =(const LeastSquaresConformalSurfaceMapper &other)
{
  if (this != &other) {
    FreeBoundarySurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
LeastSquaresConformalSurfaceMapper::~LeastSquaresConformalSurfaceMapper()
{
}

// =============================================================================
// Constraints
// =============================================================================

// -----------------------------------------------------------------------------
void LeastSquaresConformalSurfaceMapper::AddFixedPoint(int i, double u, double v)
{
  auto it = Find(_FixedPoints, i);
  if (it == _FixedPoints.end()) {
    int n = NumberOfFixedPoints();
    _FixedPoints.push_back(i);
    if (n == _FixedValues.Cols()) {
      _FixedValues.Resize(NumberOfComponents(), n + 10);
    }
    auto col = _FixedValues.Col(n);
    col[0] = u;
    col[1] = v;
  } else {
    cerr << this->NameOfType() << "::AddFixedPoint: Fixed constraint for point " << i << " already set" << endl;
    exit(1);
  }
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
double LeastSquaresConformalSurfaceMapper::Weight(int i, int j) const
{
  int    k, l;
  double a[3], b[3], c[3], w;

  const int n = GetEdgeNeighborPoints(i, j, k, l);
  if (n == 0 || n > 2) {
    cerr << this->NameOfType() << "::Weight: Surface mesh must be triangulated!" << endl;
    exit(1);
  }

  _Surface->GetPoint(static_cast<vtkIdType>(i), a);
  _Surface->GetPoint(static_cast<vtkIdType>(k), b);
  _Surface->GetPoint(static_cast<vtkIdType>(j), c);
  w = Triangle::Cotangent(a, b, c);
  if (n == 2) {
    _Surface->GetPoint(static_cast<vtkIdType>(l), b);
    w += Triangle::Cotangent(a, b, c);
  }

  return w;
}

// -----------------------------------------------------------------------------
void LeastSquaresConformalSurfaceMapper::Initialize()
{
  // Initialize base class
  FreeBoundarySurfaceMapper::Initialize();

  // Check input
  if (_Surface->GetNumberOfPolys() == 0) {
    cerr << this->NameOfType() << "::Initialize: Input point set must be a surface mesh" << endl;
    exit(1);
  }

  // Choose fixed points
  int s, i;
  if (NumberOfFixedPoints() == 0) {
    i = 0;
    s = _Boundary->FindLongestSegment();
    AddFixedPoint(_Boundary->Segment(s).PointId(i), 0., 0.);
  } else {
    s = _Boundary->FindSegment(FixedPointId(0), &i);
    if (s == -1) {
      cerr << this->NameOfType() << "::Initialize: Fixed point must be on boundary" << endl;
      exit(1);
    }
  }
  if (NumberOfFixedPoints() == 1) {
    ChordLengthBoundarySegmentParameterizer param;
    BoundarySegment segment = _Boundary->Segment(s);
    segment.ClearSelection();
    segment.SelectPoint(i);
    param.Boundary(segment);
    param.Run();
    for (i = 1; i < segment.NumberOfPoints(); ++i) {
      if (param.Value(i) >= .5) {
        AddFixedPoint(segment.PointId(i), param.Value(i) * segment.Length(), 0.);
        break;
      }
    }
  } else {
    for (i = 1; i < NumberOfFixedPoints(); ++i) {
      if (_Boundary->FindSegment(FixedPointId(i)) != s) {
        cerr << this->NameOfType() << "::Initialize: All fixed points must be on the same boundary segment" << endl;
        exit(1);
      }
    }
  }

  // Initialize map values and determine points with free map values
  const int m = NumberOfComponents();
  const int n = NumberOfPoints();

  #if MIRTK_USE_FLOAT_BY_DEFAULT
    _Values = vtkSmartPointer<vtkFloatArray>::New();
  #else
    _Values = vtkSmartPointer<vtkDoubleArray>::New();
  #endif
  _Values->SetName("SurfaceMap");
  _Values->SetNumberOfComponents(m);
  _Values->SetNumberOfTuples(n);

  _FreePoints.clear();
  _FreePoints.reserve(n - NumberOfFixedPoints());
  _PointIndex.resize(n);

  const double zero[2] = {0.};

  for (int i = 0, j; i < n; ++i) {
    j = FindIndex(_FixedPoints, i);
    if (j < 0) {
      j = static_cast<int>(_FreePoints.size());
      _FreePoints.push_back(i);
      _PointIndex[i] = j;
      _Values->SetTuple(i, zero);
    } else {
      _PointIndex[i] = -(j + 1);
      _Values->SetTuple(i, _FixedValues.Col(j));
    }
  }
}

// -----------------------------------------------------------------------------
void LeastSquaresConformalSurfaceMapper::ComputeMap()
{
  const bool use_direct_solver = (_NumberOfIterations == 1 || _NumberOfIterations < 0);

  MIRTK_START_TIMING();

  typedef Eigen::VectorXd             Vector;
  typedef Eigen::SparseMatrix<double> Matrix;
  typedef Eigen::Triplet<double>      NZEntry;

  const int n = NumberOfFreePoints();
  const int m = 2;

  int i, j, ui, vi, uj, vj;

  Matrix A(m * n, m * n);
  Vector b(m * n);
  {
    Array<NZEntry> w;
    Array<double>  w_ii(n, .0);
    double         w_ij;

    w.reserve(m * n * (_EdgeTable->MaxNumberOfAdjacentPoints() + 1));
    b.setZero();

    EdgeIterator edgeIt(*_EdgeTable);
    for (edgeIt.InitTraversal(); edgeIt.GetNextEdge(i, j) != -1;) {
      ui = FreePointIndex(i), vi = ui + n;
      uj = FreePointIndex(j), vj = uj + n;
      if (ui >= 0 || uj >= 0) {
        w_ij  = - this->Weight(i, j);
        if (ui >= 0 && uj >= 0) {
          w.push_back(NZEntry(ui, uj, w_ij));
          w.push_back(NZEntry(uj, ui, w_ij));
          w.push_back(NZEntry(vi, vj, w_ij));
          w.push_back(NZEntry(vj, vi, w_ij));
        } else if (ui >= 0) {
          b(ui) -= w_ij * GetValue(j, 0);
          b(vi) -= w_ij * GetValue(j, 1);
        } else if (uj >= 0) {
          b(uj) -= w_ij * GetValue(i, 0);
          b(vj) -= w_ij * GetValue(i, 1);
        }
        if (ui >= 0) {
          w_ii[ui] -= w_ij;
        }
        if (uj >= 0) {
          w_ii[uj] -= w_ij;
        }
      }
    }
    for (ui = 0, vi = n; ui < n; ++ui, ++vi) {
      w.push_back(NZEntry(ui, ui, w_ii[ui]));
      w.push_back(NZEntry(vi, vi, w_ii[ui]));
    }

    A.setFromTriplets(w.begin(), w.end());

    const int s = _Boundary->FindSegment(FixedPointId(0));
    const auto &segment = _Boundary->Segment(s);
    for (int k = 0; k < segment.NumberOfPoints(); ++k) {
      i = segment.PointId(k);
      j = segment.PointId(k+1);
      ui = FreePointIndex(i), vi = ui + n;
      uj = FreePointIndex(j), vj = uj + n;
      if (ui >= 0 && uj >= 0) {
        A.coeffRef(ui, vj) -= 1.;
        A.coeffRef(vj, ui) -= 1.;
        A.coeffRef(uj, vi) += 1.;
        A.coeffRef(vi, uj) += 1.;
      } else if (ui >= 0) {
        b(ui) += GetValue(j, 1);
        b(vi) -= GetValue(j, 0);
      } else if (uj >= 0) {
        b(uj) -= GetValue(i, 1);
        b(vj) += GetValue(i, 0);
      }
    }

    A.makeCompressed();
  }

  MIRTK_DEBUG_TIMING(1, "building sparse linear system");

  if (verbose) {
    cout << "\n";
    cout << "  No. of surface points        = " << NumberOfPoints() << "\n";
    cout << "  No. of fixed points          = " << NumberOfFixedPoints() << "\n";
    cout << "  No. of free points           = " << n << "\n";
    cout << "  No. of non-zero coefficients = " << A.nonZeros() << "\n";
    cout << "  Dimension of map codomain    = " << m << "\n";
    cout.flush();
  }

  MIRTK_RESET_TIMING();

  Vector x(m * n);
  int    niter = 0;
  double error = nan;

  if (use_direct_solver) {
    Eigen::SparseLU<Matrix> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    x = solver.solve(b);
  } else {
    Eigen::ConjugateGradient<Matrix> solver(A);
    if (_NumberOfIterations > 0 ) solver.setMaxIterations(_NumberOfIterations);
    if (_Tolerance          > 0.) solver.setTolerance(_Tolerance);
    x = solver.solve(b);
    niter = solver.iterations();
    error = solver.error();
  }

  for (ui = 0, vi = n; ui < n; ++ui, ++vi) {
    i = FreePointId(ui);
    SetValue(i, 0, x(ui));
    SetValue(i, 1, x(vi));
  }

  MIRTK_DEBUG_TIMING(1, "solving sparse linear system");

  if (verbose && !use_direct_solver) {
    cout << "  No. of iterations            = " << niter << "\n";
    cout << "  Estimated error              = " << error << "\n";
    cout.flush();
  }
}

// -----------------------------------------------------------------------------
void LeastSquaresConformalSurfaceMapper::Finalize()
{
  // Assemble surface map
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
  FreeBoundarySurfaceMapper::Finalize();
}


} // namespace mirtk
