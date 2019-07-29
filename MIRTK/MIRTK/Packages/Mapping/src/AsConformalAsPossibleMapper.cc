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

#include "mirtk/AsConformalAsPossibleMapper.h"

#include "mirtk/Math.h"
#include "mirtk/Array.h"
#include "mirtk/Matrix3x3.h"
#include "mirtk/Parallel.h"
#include "mirtk/HarmonicTetrahedralMeshMapper.h" // used to obtain initial map
#include "mirtk/VtkMath.h"

#include "vtkSmartPointer.h"
#include "vtkTetra.h"
#include "vtkIdList.h"

#include "Eigen/SVD"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

namespace AsConformalAsPossibleMapperUtils {


// -----------------------------------------------------------------------------
/// Compute orientation of each tetrahedron using SVD of map Jacobian
class ComputeOrientationOfTetrahedra
{
public:

  typedef Eigen::Matrix<double, 3, 3>                              EigenMatrix;
  typedef Eigen::JacobiSVD<EigenMatrix, Eigen::NoQRPreconditioner> SVDSolver;

  vtkPointSet      *_PointSet;
  vtkDataArray     *_Coords;
  Array<Matrix3x3> *_Orientation;

  // ---------------------------------------------------------------------------
  /// Add tensor product of 3D vectors to given 3x3 matrix
  static inline void AddTensorProduct(const double a[3], const double b[3], EigenMatrix &m)
  {
    m(0, 0) += a[0] * b[0];
    m(0, 1) += a[0] * b[1];
    m(0, 2) += a[0] * b[2];
    m(1, 0) += a[1] * b[0];
    m(1, 1) += a[1] * b[1];
    m(1, 2) += a[1] * b[2];
    m(2, 0) += a[2] * b[0];
    m(2, 1) += a[2] * b[1];
    m(2, 2) += a[2] * b[2];
  }

  // ---------------------------------------------------------------------------
  /// Determinant of 3x3 matrix
  static inline double Det(const EigenMatrix &m)
  {
    return m(0, 0) * m(1, 1) * m(2, 2) +
           m(0, 1) * m(1, 2) * m(2, 0) +
           m(0, 2) * m(1, 0) * m(2, 1) -
           m(0, 2) * m(1, 1) * m(2, 0) -
           m(0, 1) * m(1, 0) * m(2, 2) -
           m(0, 0) * m(1, 2) * m(2, 1);
  }

  // ---------------------------------------------------------------------------
  /// Compute local orientation of each tetrahedron via SVD of Jacobian
  void operator()(const blocked_range<vtkIdType> &re) const
  {
    vtkIdType   i0, i1, i2, i3;
    double      v0[3], v1[3], v2[3], v3[3], a[3], b[3], c[3], n[3];
    EigenMatrix jac, rot, inv;
    SVDSolver   svd;

    inv.setIdentity();
    inv(2, 2) = -1.0;

    vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
    for (vtkIdType cellId = re.begin(); cellId != re.end(); ++cellId) {

      // Get indices of cell points
      _PointSet->GetCellPoints(cellId, ptIds);
      i0 = ptIds->GetId(0);
      i1 = ptIds->GetId(1);
      i2 = ptIds->GetId(2);
      i3 = ptIds->GetId(3);

      // Get input domain coordinates of cell points
      _PointSet->GetPoint(i0, v0);
      _PointSet->GetPoint(i1, v1);
      _PointSet->GetPoint(i2, v2);
      _PointSet->GetPoint(i3, v3);

      // Compute Jacobian of volumetric map
      jac.setZero();

      c[0] = _Coords->GetComponent(i0, 0);
      c[1] = _Coords->GetComponent(i0, 1);
      c[2] = _Coords->GetComponent(i0, 2);

      vtkMath::Subtract(v2, v1, a);
      vtkMath::Subtract(v3, v1, b);
      vtkMath::Cross(a, b, n);
      vtkMath::MultiplyScalar(n, .5);
      AddTensorProduct(c, n, jac);

      c[0] = _Coords->GetComponent(i1, 0);
      c[1] = _Coords->GetComponent(i1, 1);
      c[2] = _Coords->GetComponent(i1, 2);

      vtkMath::Subtract(v0, v2, a);
      vtkMath::Subtract(v3, v2, b);
      vtkMath::Cross(a, b, n);
      vtkMath::MultiplyScalar(n, .5);
      AddTensorProduct(c, n, jac);

      c[0] = _Coords->GetComponent(i2, 0);
      c[1] = _Coords->GetComponent(i2, 1);
      c[2] = _Coords->GetComponent(i2, 2);

      vtkMath::Subtract(v0, v3, a);
      vtkMath::Subtract(v1, v3, b);
      vtkMath::Cross(a, b, n);
      vtkMath::MultiplyScalar(n, .5);
      AddTensorProduct(c, n, jac);

      c[0] = _Coords->GetComponent(i3, 0);
      c[1] = _Coords->GetComponent(i3, 1);
      c[2] = _Coords->GetComponent(i3, 2);

      vtkMath::Subtract(v2, v0, a);
      vtkMath::Subtract(v1, v0, b);
      vtkMath::Cross(a, b, n);
      vtkMath::MultiplyScalar(n, .5);
      AddTensorProduct(c, n, jac);

      jac /= -3.0 * vtkTetra::ComputeVolume(v0, v1, v2, v3);

      // Compute principle directions using SVD
      svd.compute(jac, Eigen::ComputeFullU | Eigen::ComputeFullV);
      rot = svd.matrixU() * svd.matrixV().transpose();
      if (Det(rot) < 0) rot = svd.matrixU() * inv * svd.matrixV().transpose();

      // Set local orientation matrix
      (*_Orientation)[cellId] = Matrix3x3(rot(0, 0), rot(0, 1), rot(0, 2),
                                          rot(1, 0), rot(1, 1), rot(1, 2),
                                          rot(2, 0), rot(2, 1), rot(2, 2));
    }
  }
};


}
// namespace AsConformalAsPossibleMapperUtils
using namespace AsConformalAsPossibleMapperUtils;

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void AsConformalAsPossibleMapper
::CopyAttributes(const AsConformalAsPossibleMapper &other)
{
  _UniformWeight = other._UniformWeight;
  _Orientation   = other._Orientation;
}

// -----------------------------------------------------------------------------
AsConformalAsPossibleMapper::AsConformalAsPossibleMapper()
:
  _UniformWeight(.7)
{
}

// -----------------------------------------------------------------------------
AsConformalAsPossibleMapper
::AsConformalAsPossibleMapper(const AsConformalAsPossibleMapper &other)
:
  LinearTetrahedralMeshMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
AsConformalAsPossibleMapper &AsConformalAsPossibleMapper
::operator =(const AsConformalAsPossibleMapper &other)
{
  if (this != &other) {
    LinearTetrahedralMeshMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
AsConformalAsPossibleMapper::~AsConformalAsPossibleMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void AsConformalAsPossibleMapper::Initialize()
{
  const int d = 3; // dimension of output domain

  // Initialize base class
  LinearTetrahedralMeshMapper::Initialize();

  // Limit uniform scaling weight to interval [0, 1]
  _UniformWeight = max(.0, min(1.0, _UniformWeight));

  // Check if initial volumetric map is given
  bool no_initial_map = true;
  for (int i = 0; i < _NumberOfInteriorPoints; ++i) {
    for (int j = 0; j < d; ++j) {
      if (!fequal(_Coords->GetComponent(_InteriorPointId[i], j), .0)) {
        no_initial_map = false;
        i = _NumberOfInteriorPoints; // break outer loop
        break;
      }
    }
  }

  // Obtain initial harmonic map
  if (no_initial_map) {
    HarmonicTetrahedralMeshMapper harmonic_map;
    Solve(&harmonic_map);
  }

  // Compute local orientation of each tetrahedron
  _Orientation.resize(_Volume->GetNumberOfCells());
  ComputeOrientationOfTetrahedra eval;
  eval._PointSet    = _Volume;
  eval._Coords      = _Coords;
  eval._Orientation = &_Orientation;
  blocked_range<vtkIdType> cellIds(0, _Volume->GetNumberOfCells());
  parallel_for(cellIds, eval);
}

// -----------------------------------------------------------------------------
void AsConformalAsPossibleMapper::Finalize()
{
  // Finalize base class
  LinearTetrahedralMeshMapper::Finalize();

  // Clear orientation matrices
  _Orientation.clear();
}

// =============================================================================
// Auxiliary functions
// =============================================================================

// -----------------------------------------------------------------------------
static inline Matrix3x3 GetScaleMatrix(const double n[3], double w)
{
  const double wn[3] = { w * n[0], w * n[1], w * n[2] };
  return Matrix3x3(    .0,  wn[1], -wn[2],
                   -wn[0],     .0,  wn[2],
                    wn[0], -wn[1],     .0);
}

// -----------------------------------------------------------------------------
static inline Matrix3x3 GetAngleMatrix(const double n[3], double w)
{
  const double wn[3] = { w * n[0], w * n[1], w * n[2] };
  return Matrix3x3(   .0, wn[2], wn[1],
                   wn[2],    .0, wn[0],
                   wn[1], wn[0],    .0);
}

// -----------------------------------------------------------------------------
Matrix3x3 AsConformalAsPossibleMapper
::GetWeight(vtkIdType cellId, const double v0[3], const double v1[3],
                              const double v2[3], const double v3[3], double volume) const
{
  // Note: Factor 2 is "pre-multiplied" by 1/2 factor of edge cross product
  const double &scale_weight = _UniformWeight;
  const double  angle_weight = 1.0 - _UniformWeight;

  double a[3], b[3], n0[3], n1[3];

  vtkMath::Subtract(v2, v1, a);
  vtkMath::Subtract(v3, v1, b);
  vtkMath::Cross(a, b, n0);

  vtkMath::Subtract(v3, v0, a);
  vtkMath::Subtract(v2, v0, b);
  vtkMath::Cross(a, b, n1);

  Matrix3x3 scale0 = GetScaleMatrix(n0, scale_weight);
  Matrix3x3 scale1 = GetScaleMatrix(n1, scale_weight);
  Matrix3x3 angle0 = GetAngleMatrix(n0, angle_weight);
  Matrix3x3 angle1 = GetAngleMatrix(n1, angle_weight);

  Matrix3x3 dTd = (scale0.Transpose() * scale1 + angle0.Transpose() * angle1);
  dTd /= 9.0 * volume;

  return _Orientation[cellId] * dTd * _Orientation[cellId].Transpose();
}


} // namespace mirtk
