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

#ifndef MIRTK_LeastSquaresConformalSurfaceMapper_H
#define MIRTK_LeastSquaresConformalSurfaceMapper_H

#include "mirtk/FreeBoundarySurfaceMapper.h"

#include "mirtk/Array.h"

#include "vtkSmartPointer.h"
#include "vtkDataArray.h"


namespace mirtk {


/**
 * Least squares conformal parameterization
 *
 * This filter implements the least squares conformal map (LSCM) proposed by
 * Levy (2002) and the identical discrete natural conformal parameterization
 * (DNCP) proposed by Meyer et al. (2002). It follows the formulation of in
 * Mullen et al. (2002), where the LSCM is defined as the minimizer of the
 * conformal energy \f$E_C(\vec{u}) = E_LSCM\vec{u} = E_D(\vec{u}) - A(\vec{u})\f$,
 * where \f$E_D(\vec{u})\f$ is the discrete Dirichlet energy and
 * \f$A(\vec{u})\f$ is the total area of the parameteric domain.
 * The geometric cotangent weights are used to discretize the Laplace operator
 * needed for the computation of the derivative of the Dirichlet energy.
 *
 * The LSCM requires at least two fixed points. When no point constraints are
 * given, a boundary point is selected automatically. A second fixed boundary
 * point is chosen to be as farthest away from the first fixed point as possible
 * based on geodesic distances on the input surface mesh.
 *
 * - Levy et al. (2002). Least squares conformal maps for automatic texture atlas
 *   generation. ACM Trans. Graphics, 21(3), 362–371.
 * - Desbrun, Meyer, and Alliez (2002). Intrinsic parameterizations of surface meshes.
 *   Computer Graphics Forum, 21(3), 209–218.
 * - Mullen et al. (2008). Spectral conformal parameterization.
 *   Eurographics Symposium on Geometry Processing, 27(5), 1487–1494.
 *
 * \todo Implement area weighting extension as described in Mullen et al. (2008)
 *       to account for irregular surface sampling.
 */
class LeastSquaresConformalSurfaceMapper : public FreeBoundarySurfaceMapper
{
  mirtkObjectMacro(LeastSquaresConformalSurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

private:

  /// Maximum number of iterations
  ///
  /// When the number of iterations is set to 1 a sparse direct solver is used.
  mirtkPublicAttributeMacro(int, NumberOfIterations);

  /// Tolerance for sparse linear solver
  mirtkPublicAttributeMacro(double, Tolerance);

  /// Index of point in set of points with free (i >= 0) or fixed (i < 0) values
  mirtkAttributeMacro(Array<int>, PointIndex);

  /// IDs of surface points with free map values
  mirtkAttributeMacro(Array<int>, FreePoints);

  /// IDs of surface points with fixed map values
  mirtkAttributeMacro(Array<int>, FixedPoints);

  /// Map values at fixed points
  mirtkAttributeMacro(Matrix, FixedValues);

  /// Computed map values at surface points
  mirtkAttributeMacro(vtkSmartPointer<vtkDataArray>, Values);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const LeastSquaresConformalSurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  LeastSquaresConformalSurfaceMapper();

  /// Copy constructor
  LeastSquaresConformalSurfaceMapper(const LeastSquaresConformalSurfaceMapper &);

  /// Assignment operator
  LeastSquaresConformalSurfaceMapper &operator =(const LeastSquaresConformalSurfaceMapper &);

  /// Destructor
  virtual ~LeastSquaresConformalSurfaceMapper();

  // ---------------------------------------------------------------------------
  // Constraints

  /// Add hard point constraint
  ///
  /// \param[in] i Surface point index.
  /// \param[in] u First  constraint value component.
  /// \param[in] v Second constraint value component.
  void AddFixedPoint(int i, double u, double v);

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Weight of undirected edge (i, j)
  ///
  /// \param[in] i First end point.
  /// \param[in] j Second end point.
  ///
  /// \returns Weight of undirected edge (i, j).
  virtual double Weight(int i, int j) const;

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Compute surface map
  virtual void ComputeMap();

  /// Finalize filter execution
  virtual void Finalize();

  // ---------------------------------------------------------------------------
  // Auxiliaries

public:

  /// Number of surface points with free map value
  int NumberOfFreePoints() const;

  /// Number of surface points with fixed map value
  int NumberOfFixedPoints() const;

  /// Whether the map value of the specified surface point is free
  bool IsFreePoint(int i) const;

  /// Whether the map value of the specified surface point is fixed
  bool IsFixedPoint(int i) const;

  /// Get index of i-th point with free map value or -1 if map value of point is fixed
  ///
  /// \param[in] i Surface point index.
  ///
  /// \return Index of point in set of points with free map value.
  int FreePointIndex(int i) const;

  /// Get surface point ID of i-th point with free map value
  ///
  /// \param[in] i Index of point in set of points with free map value.
  ///
  /// \return Surface point index.
  int FreePointId(int i) const;

  /// Get index of i-th point with fixed map value or -1 if map value of point is free
  ///
  /// \param[in] i Surface point index.
  ///
  /// \return Index of point in set of points with fixed map value.
  int FixedPointIndex(int i) const;

  /// Get surface point ID of i-th point with fixed map value
  ///
  /// \param[in] i Index of point in set of points with fixed map value.
  ///
  /// \return Surface point index.
  int FixedPointId(int i) const;

  /// Get component of map value at surface vertex
  ///
  /// \param[in] i Surface point index.
  /// \param[in] j Map value component index.
  ///
  /// \return The j-th component of the map value evaluated at the i-th surface point.
  double GetValue(int i, int j = 0) const;

  /// Get component of map value at i-th free point
  ///
  /// \param[in] i Index of point in set of points with free map value.
  /// \param[in] j Map value component index.
  ///
  /// \return The j-th component of the map value evaluated at the i-th free point.
  double GetFreeValue(int i, int j = 0) const;

  /// Get component of map value at i-th fixed point
  ///
  /// \param[in] i Index of point in set of points with fixed map value.
  /// \param[in] j Map value component index.
  ///
  /// \return The j-th component of the map value evaluated at the i-th fixed point.
  double GetFixedValue(int i, int j = 0) const;

protected:

  /// Set scalar map value at surface vertex
  ///
  /// \param[in] i Surface point index.
  /// \param[in] v Map value.
  void SetValue(int i, double v);

  /// Set component of map value at surface vertex
  ///
  /// \param[in] i Surface point index.
  /// \param[in] j Map component index.
  /// \param[in] v Map component value.
  void SetValue(int i, int j, double v);

  /// Set scalar map value at i-th free point
  ///
  /// \param[in] i Index of point in set of points with free map value.
  /// \param[in] v Map value.
  void SetFreeValue(int i, double v);

  /// Set component of map value at i-th free point
  ///
  /// \param[in] i Index of point in set of points with free map value.
  /// \param[in] j Map component index.
  /// \param[in] v Map component value.
  void SetFreeValue(int i, int j, double v);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
inline int LeastSquaresConformalSurfaceMapper::NumberOfFreePoints() const
{
  return static_cast<int>(_FreePoints.size());
}

// -----------------------------------------------------------------------------
inline int LeastSquaresConformalSurfaceMapper::NumberOfFixedPoints() const
{
  return static_cast<int>(_FixedPoints.size());
}

// -----------------------------------------------------------------------------
inline int LeastSquaresConformalSurfaceMapper::FreePointIndex(int ptId) const
{
  const int i = _PointIndex[ptId];
  return (i < 0 ? -1 : i);
}

// -----------------------------------------------------------------------------
inline int LeastSquaresConformalSurfaceMapper::FreePointId(int i) const
{
  return _FreePoints[i];
}

// -----------------------------------------------------------------------------
inline bool LeastSquaresConformalSurfaceMapper::IsFreePoint(int ptId) const
{
  return FreePointIndex(ptId) != -1;
}

// -----------------------------------------------------------------------------
inline int LeastSquaresConformalSurfaceMapper::FixedPointIndex(int ptId) const
{
  const int i = _PointIndex[ptId];
  return (i < 0 ? (-i) - 1 : -1);
}

// -----------------------------------------------------------------------------
inline int LeastSquaresConformalSurfaceMapper::FixedPointId(int i) const
{
  return _FixedPoints[i];
}

// -----------------------------------------------------------------------------
inline bool LeastSquaresConformalSurfaceMapper::IsFixedPoint(int ptId) const
{
  return FixedPointIndex(ptId) != -1;
}

// -----------------------------------------------------------------------------
inline double LeastSquaresConformalSurfaceMapper::GetValue(int i, int j) const
{
  return _Values->GetComponent(static_cast<vtkIdType>(i), j);
}

// -----------------------------------------------------------------------------
inline double LeastSquaresConformalSurfaceMapper::GetFreeValue(int i, int j) const
{
  return GetValue(FreePointId(i), j);
}

// -----------------------------------------------------------------------------
inline double LeastSquaresConformalSurfaceMapper::GetFixedValue(int i, int j) const
{
  return GetValue(FixedPointId(i), j);
}

// -----------------------------------------------------------------------------
inline void LeastSquaresConformalSurfaceMapper::SetValue(int i, double v)
{
  _Values->SetComponent(static_cast<vtkIdType>(i), 0, v);
}

// -----------------------------------------------------------------------------
inline void LeastSquaresConformalSurfaceMapper::SetValue(int i, int j, double v)
{
  _Values->SetComponent(static_cast<vtkIdType>(i), j, v);
}

// -----------------------------------------------------------------------------
inline void LeastSquaresConformalSurfaceMapper::SetFreeValue(int i, double v)
{
  SetValue(FreePointId(i), v);
}

// -----------------------------------------------------------------------------
inline void LeastSquaresConformalSurfaceMapper::SetFreeValue(int i, int j, double v)
{
  SetValue(FreePointId(i), j, v);
}


} // namespace mirtk

#endif // MIRTK_LeastSquaresConformalSurfaceMapper_H
