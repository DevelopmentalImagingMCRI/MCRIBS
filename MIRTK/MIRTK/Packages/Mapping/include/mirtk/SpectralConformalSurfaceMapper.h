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

#ifndef MIRTK_SpectralConformalSurfaceMapper_H
#define MIRTK_SpectralConformalSurfaceMapper_H

#include "mirtk/FreeBoundarySurfaceMapper.h"

#include "mirtk/Array.h"

#include "vtkSmartPointer.h"
#include "vtkDataArray.h"


namespace mirtk {


/**
 * Least squares conformal map without boundary constraints
 *
 * - Mullen et al. (2008). Spectral conformal parameterization.
 *   Eurographics Symposium on Geometry Processing, 27(5), 1487–1494.
 *
 * \todo Implement area weighting extension as described in Mullen et al. (2008)
 *       to account for irregular surface sampling.
 */
class SpectralConformalSurfaceMapper : public FreeBoundarySurfaceMapper
{
  mirtkObjectMacro(SpectralConformalSurfaceMapper);

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
  void CopyAttributes(const SpectralConformalSurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  SpectralConformalSurfaceMapper();

  /// Copy constructor
  SpectralConformalSurfaceMapper(const SpectralConformalSurfaceMapper &);

  /// Assignment operator
  SpectralConformalSurfaceMapper &operator =(const SpectralConformalSurfaceMapper &);

  /// Destructor
  virtual ~SpectralConformalSurfaceMapper();

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
inline int SpectralConformalSurfaceMapper::NumberOfFreePoints() const
{
  return static_cast<int>(_FreePoints.size());
}

// -----------------------------------------------------------------------------
inline int SpectralConformalSurfaceMapper::NumberOfFixedPoints() const
{
  return static_cast<int>(_FixedPoints.size());
}

// -----------------------------------------------------------------------------
inline int SpectralConformalSurfaceMapper::FreePointIndex(int ptId) const
{
  const int i = _PointIndex[ptId];
  return (i < 0 ? -1 : i);
}

// -----------------------------------------------------------------------------
inline int SpectralConformalSurfaceMapper::FreePointId(int i) const
{
  return _FreePoints[i];
}

// -----------------------------------------------------------------------------
inline bool SpectralConformalSurfaceMapper::IsFreePoint(int ptId) const
{
  return FreePointIndex(ptId) != -1;
}

// -----------------------------------------------------------------------------
inline int SpectralConformalSurfaceMapper::FixedPointIndex(int ptId) const
{
  const int i = _PointIndex[ptId];
  return (i < 0 ? (-i) - 1 : -1);
}

// -----------------------------------------------------------------------------
inline int SpectralConformalSurfaceMapper::FixedPointId(int i) const
{
  return _FixedPoints[i];
}

// -----------------------------------------------------------------------------
inline bool SpectralConformalSurfaceMapper::IsFixedPoint(int ptId) const
{
  return FixedPointIndex(ptId) != -1;
}

// -----------------------------------------------------------------------------
inline double SpectralConformalSurfaceMapper::GetValue(int i, int j) const
{
  return _Values->GetComponent(static_cast<vtkIdType>(i), j);
}

// -----------------------------------------------------------------------------
inline double SpectralConformalSurfaceMapper::GetFreeValue(int i, int j) const
{
  return GetValue(FreePointId(i), j);
}

// -----------------------------------------------------------------------------
inline double SpectralConformalSurfaceMapper::GetFixedValue(int i, int j) const
{
  return GetValue(FixedPointId(i), j);
}

// -----------------------------------------------------------------------------
inline void SpectralConformalSurfaceMapper::SetValue(int i, double v)
{
  _Values->SetComponent(static_cast<vtkIdType>(i), 0, v);
}

// -----------------------------------------------------------------------------
inline void SpectralConformalSurfaceMapper::SetValue(int i, int j, double v)
{
  _Values->SetComponent(static_cast<vtkIdType>(i), j, v);
}

// -----------------------------------------------------------------------------
inline void SpectralConformalSurfaceMapper::SetFreeValue(int i, double v)
{
  SetValue(FreePointId(i), v);
}

// -----------------------------------------------------------------------------
inline void SpectralConformalSurfaceMapper::SetFreeValue(int i, int j, double v)
{
  SetValue(FreePointId(i), j, v);
}


} // namespace mirtk

#endif // MIRTK_SpectralConformalSurfaceMapper_H
