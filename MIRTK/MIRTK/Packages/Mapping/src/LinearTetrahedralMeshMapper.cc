/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
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

#include "mirtk/LinearTetrahedralMeshMapper.h"

#include "mirtk/Array.h"
#include "mirtk/Parallel.h"
#include "mirtk/Matrix3x3.h"
#include "mirtk/VtkMath.h"

#include "vtkSmartPointer.h"
#include "vtkPointSet.h"
#include "vtkPointData.h"
#include "vtkIdList.h"
#include "vtkTetra.h"

#include "Eigen/SparseCore"
#include "Eigen/IterativeLinearSolvers"


namespace mirtk {


// Global flags (cf. mirtk/Options.h)
MIRTK_Common_EXPORT extern int verbose;


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace LinearTetrahedralMeshMapperUtils {


// -----------------------------------------------------------------------------
template <class Scalar>
class LinearSystem
{
  const LinearTetrahedralMeshMapper        *_Filter;
  const LinearTetrahedralMeshMapper        *_Operator;
  Array<Eigen::Triplet<Scalar> >            _Coefficients;
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1>  _RightHandSide;

public:

  // ---------------------------------------------------------------------------
  LinearSystem(const LinearTetrahedralMeshMapper *filter,
               const LinearTetrahedralMeshMapper *map, int n)
  :
    _Filter(filter), _Operator(map ? map : filter), _RightHandSide(n)
  {
    _RightHandSide.setZero();
  }

  LinearSystem(const LinearSystem &other, split)
  :
    _Filter(other._Filter), _Operator(other._Operator),
    _RightHandSide(other._RightHandSide.rows())
  {
    _RightHandSide.setZero();
  }

  // ---------------------------------------------------------------------------
  void join(const LinearSystem &other)
  {
    _Coefficients.insert(_Coefficients.end(), other._Coefficients.begin(), other._Coefficients.end());
    _RightHandSide += other._RightHandSide;
  }

  // ---------------------------------------------------------------------------
  void AddWeight(vtkIdType ptId0, bool isBoundary0, vtkIdType ptId1, bool isBoundary1, const Matrix3x3 &weight)
  {
    if (isBoundary0 && isBoundary1) {

      // Unused coefficients

    } else if (isBoundary0) {

      // Point variables base index
      const int c = _Filter->InteriorPointPos()[ptId1];

      // Pre-multiply coefficient by constant boundary coordinates
      for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        const double &w = weight[i][j];
        if (w == .0) continue;
        _RightHandSide(c + j) -= w * _Filter->Coords()->GetComponent(ptId0, i);
        _Coefficients.push_back(Eigen::Triplet<Scalar>(c + j, c + i, -w));
      }

    } else if (isBoundary1) {

      // Point variables base index
      const int r = _Filter->InteriorPointPos()[ptId0];

      // Pre-multiply coefficient by constant boundary coordinates
      for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        const double &w = weight[i][j];
        if (w == .0) continue;
        _RightHandSide(r + i) -= w * _Filter->Coords()->GetComponent(ptId1, j);
        _Coefficients.push_back(Eigen::Triplet<Scalar>(r + i, r + j, -w));
      }

    } else {

      // Point variables base indices
      const int r = _Filter->InteriorPointPos()[ptId0];
      const int c = _Filter->InteriorPointPos()[ptId1];

      // Add symmetric coefficients
      for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        const double &w = weight[i][j];
        if (w == .0) continue;
        _Coefficients.push_back(Eigen::Triplet<Scalar>(r + i, c + j,  w));
        _Coefficients.push_back(Eigen::Triplet<Scalar>(r + i, r + j, -w));
        _Coefficients.push_back(Eigen::Triplet<Scalar>(c + j, r + i,  w));
        _Coefficients.push_back(Eigen::Triplet<Scalar>(c + j, c + i, -w));
      }

    }
  }

  // ---------------------------------------------------------------------------
  void operator ()(const blocked_range<vtkIdType> &cellIds)
  {
    vtkIdType i0, i1, i2, i3;
    bool      b0, b1, b2, b3;
    double    v0[3], v1[3], v2[3], v3[3], volume;

    vtkPointSet * const pointset = _Filter->Volume();
    vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

    for (vtkIdType cellId = cellIds.begin(); cellId != cellIds.end(); ++cellId) {
      pointset->GetCellPoints(cellId, ptIds);

      i0 = ptIds->GetId(0);
      i1 = ptIds->GetId(1);
      i2 = ptIds->GetId(2);
      i3 = ptIds->GetId(3);

      b0 = _Filter->IsBoundaryPoint(i0);
      b1 = _Filter->IsBoundaryPoint(i1);
      b2 = _Filter->IsBoundaryPoint(i2);
      b3 = _Filter->IsBoundaryPoint(i3);

      pointset->GetPoint(i0, v0);
      pointset->GetPoint(i1, v1);
      pointset->GetPoint(i2, v2);
      pointset->GetPoint(i3, v3);

      volume = vtkTetra::ComputeVolume(v0, v1, v2, v3);

      AddWeight(i0, b0, i1, b1, _Operator->GetWeight(cellId, v0, v1, v2, v3, volume));
      AddWeight(i0, b0, i2, b2, _Operator->GetWeight(cellId, v0, v2, v3, v1, volume));
      AddWeight(i0, b0, i3, b3, _Operator->GetWeight(cellId, v0, v3, v1, v2, volume));
      AddWeight(i1, b1, i2, b2, _Operator->GetWeight(cellId, v1, v2, v0, v3, volume));
      AddWeight(i1, b1, i3, b3, _Operator->GetWeight(cellId, v1, v3, v2, v0, volume));
      AddWeight(i2, b2, i3, b3, _Operator->GetWeight(cellId, v2, v3, v0, v1, volume));
    }
  }

  // ---------------------------------------------------------------------------
  static void Build(const LinearTetrahedralMeshMapper        *filter,
                    const LinearTetrahedralMeshMapper        *mapop,
                    Eigen::SparseMatrix<Scalar>              &A,
                    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &b,
                    int                                       n)
  {
    LinearSystem problem(filter, mapop, n);
    blocked_range<vtkIdType> cellIds(0, filter->Volume()->GetNumberOfCells());
    parallel_reduce(cellIds, problem);
    A.resize(n, n);
    A.setFromTriplets(problem._Coefficients.begin(), problem._Coefficients.end());
    b = problem._RightHandSide;
  }
};


} // namespace LinearTetrahedralMeshMapperUtils
using namespace LinearTetrahedralMeshMapperUtils;

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void LinearTetrahedralMeshMapper
::CopyAttributes(const LinearTetrahedralMeshMapper &other)
{
  _NumberOfIterations = other._NumberOfIterations;
  _Tolerance          = other._Tolerance;
  _RelaxationFactor   = other._RelaxationFactor;
  _InteriorPointId    = other._InteriorPointId;
  _InteriorPointPos   = other._InteriorPointPos;
}

// -----------------------------------------------------------------------------
LinearTetrahedralMeshMapper::LinearTetrahedralMeshMapper()
:
  _NumberOfIterations(0),
  _Tolerance(.0),
  _RelaxationFactor(1.0)
{
}

// -----------------------------------------------------------------------------
LinearTetrahedralMeshMapper::LinearTetrahedralMeshMapper(const LinearTetrahedralMeshMapper &other)
:
  TetrahedralMeshMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
LinearTetrahedralMeshMapper &LinearTetrahedralMeshMapper
::operator =(const LinearTetrahedralMeshMapper &other)
{
  if (this != &other) {
    TetrahedralMeshMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
LinearTetrahedralMeshMapper::~LinearTetrahedralMeshMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void LinearTetrahedralMeshMapper::Initialize()
{
  const int dim = 3; // Dimension of output domain

  // Initialize base class
  TetrahedralMeshMapper::Initialize();

  // Force number of output map components to be equal to dim
  if (_Coords->GetNumberOfComponents() != dim) {
    vtkSmartPointer<vtkDataArray> coords;
    coords.TakeReference(_Coords->NewInstance());
    coords->SetNumberOfComponents(dim);
    coords->SetNumberOfTuples(_Coords->GetNumberOfTuples());
    coords->SetName(_Coords->GetName());
    for (vtkIdType i = 0; i < coords->GetNumberOfTuples(); ++i) {
      for (int j = 0; j < dim; ++j) {
        coords->SetComponent(i, j, _Coords->GetComponent(i, j));
      }
    }
    vtkPointData *volumePD = _Volume->GetPointData();
    for (int i = 0; i < volumePD->GetNumberOfArrays(); ++i) {
      if (volumePD->GetArray(i) == _Coords) {
        volumePD->RemoveArray(i);
        break;
      }
    }
    volumePD->AddArray(coords);
    _Coords = coords;
  }

  // Pre-compute maps from interior point index to global point ID as well as
  // the position of interior points (given its global ID) in the linear system
  _InteriorPointId .resize(_NumberOfInteriorPoints, -1);
  _InteriorPointPos.resize(_NumberOfPoints,         -1);
  for (int ptId = 0, i = 0; ptId < _NumberOfPoints; ++ptId) {
    if (IsBoundaryPoint(ptId)) continue;
    _InteriorPointId [i   ] = ptId;
    _InteriorPointPos[ptId] = dim * i;
    ++i;
  }
}

// -----------------------------------------------------------------------------
void LinearTetrahedralMeshMapper::Solve()
{
  Solve(this);
}

// -----------------------------------------------------------------------------
void LinearTetrahedralMeshMapper
::Solve(const LinearTetrahedralMeshMapper *mapop)
{
  typedef double                                   Scalar;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
  typedef Eigen::SparseMatrix<Scalar>              Matrix;
  typedef Eigen::DiagonalPreconditioner<Scalar>    Preconditioner;

  const int dim = 3;                             // Dimension of output domain
  const int n   = dim * _NumberOfInteriorPoints; // Size of linear system

  Matrix A;
  Vector x, b;

  // Use current parameterization of interior points as initial guess
  x.resize(n);
  for (int i = 0, r = 0; i < _NumberOfInteriorPoints; ++i) {
    for (int j = 0; j < dim; ++j, ++r) {
      x(r) = static_cast<Scalar>(_Coords->GetComponent(_InteriorPointId[i], j));
    }
  }

  // Build linear system
  if (verbose) cout << "\nBuilding linear system...", cout.flush();
  LinearSystem<Scalar>::Build(this, mapop, A, b, n);
  if (verbose) cout << " done" << endl;

  // Solve linear system
  if (verbose) cout << "Solve system using conjugate gradient...", cout.flush();
  Eigen::ConjugateGradient<Matrix, Eigen::Upper|Eigen::Lower, Preconditioner> solver(A);
  if (_NumberOfIterations >  0) solver.setMaxIterations(_NumberOfIterations);
  if (_Tolerance          > .0) solver.setTolerance(_Tolerance);
  x = solver.solveWithGuess(b, x);
  if (verbose) {
    cout << " done" << endl;
    cout << "\nNo. of iterations = " << solver.iterations();
    cout << "\nEstimated error   = " << solver.error();
    cout << endl;
  }

  // Update parameterization of interior points
  for (int i = 0, r = 0; i < _NumberOfInteriorPoints; ++i, r += dim) {
    for (int j = 0; j < dim; ++j) {
      _Coords->SetComponent(_InteriorPointId[i], j, static_cast<double>(x(r + j)));
    }
  }
}


} // namespace mirtk
