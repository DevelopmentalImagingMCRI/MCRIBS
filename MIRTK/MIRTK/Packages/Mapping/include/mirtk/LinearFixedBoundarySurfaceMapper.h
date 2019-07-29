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

#ifndef MIRTK_LinearFixedBoundarySurfaceMapper_H
#define MIRTK_LinearFixedBoundarySurfaceMapper_H

#include "mirtk/FixedBoundarySurfaceMapper.h"

#include "mirtk/Array.h"

#include "vtkSmartPointer.h"
#include "vtkDataArray.h"


namespace mirtk {


/**
 * Base class of discrete piecewise linear fixed boundary surface mapping methods
 *
 * For a given choice of edge weights, a sparse system of linear equations is solved
 * iteratively to obtain a surface parameterization (Marchandise et al., 2014).
 * Many published surface parameterization/texture mapping methods can be seen 
 * as spring and mass model, where the mesh vertices are the masses and the mesh
 * edges are linear springs. The only difference between the different methods
 * are the spring constants used, which result from a discretization/modeling
 * of different map energies/properties.
 *
 * \todo Implement interactive boundary map update as proposed in
 *       Desbrun et al. (2002). Intrinsic parameterizations of surface meshes.
 *
 * \todo Implement boundary map optimization as proposed in
 *       Desbrun et al. (2002). Intrinsic parameterizations of surface meshes.
 */
class LinearFixedBoundarySurfaceMapper : public FixedBoundarySurfaceMapper
{
  mirtkAbstractMacro(LinearFixedBoundarySurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

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

  /// Computed map values at surface points
  mirtkAttributeMacro(vtkSmartPointer<vtkDataArray>, Values);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const LinearFixedBoundarySurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  LinearFixedBoundarySurfaceMapper();

  /// Copy constructor
  LinearFixedBoundarySurfaceMapper(const LinearFixedBoundarySurfaceMapper &);

  /// Assignment operator
  LinearFixedBoundarySurfaceMapper &operator =(const LinearFixedBoundarySurfaceMapper &);

public:

  /// Destructor
  virtual ~LinearFixedBoundarySurfaceMapper();

  // ---------------------------------------------------------------------------
  // Execution

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Compute map values at interior points
  virtual void ComputeMap() = 0;

  /// Assemble output surface map
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
inline int LinearFixedBoundarySurfaceMapper::NumberOfFreePoints() const
{
  return static_cast<int>(_FreePoints.size());
}

// -----------------------------------------------------------------------------
inline int LinearFixedBoundarySurfaceMapper::NumberOfFixedPoints() const
{
  return static_cast<int>(_FixedPoints.size());
}

// -----------------------------------------------------------------------------
inline int LinearFixedBoundarySurfaceMapper::FreePointIndex(int ptId) const
{
  const int i = _PointIndex[ptId];
  return (i < 0 ? -1 : i);
}

// -----------------------------------------------------------------------------
inline int LinearFixedBoundarySurfaceMapper::FreePointId(int i) const
{
  return _FreePoints[i];
}

// -----------------------------------------------------------------------------
inline bool LinearFixedBoundarySurfaceMapper::IsFreePoint(int ptId) const
{
  return FreePointIndex(ptId) != -1;
}

// -----------------------------------------------------------------------------
inline int LinearFixedBoundarySurfaceMapper::FixedPointIndex(int ptId) const
{
  const int i = _PointIndex[ptId];
  return (i < 0 ? (-i) - 1 : -1);
}

// -----------------------------------------------------------------------------
inline int LinearFixedBoundarySurfaceMapper::FixedPointId(int i) const
{
  return _FixedPoints[i];
}

// -----------------------------------------------------------------------------
inline bool LinearFixedBoundarySurfaceMapper::IsFixedPoint(int ptId) const
{
  return FixedPointIndex(ptId) != -1;
}

// -----------------------------------------------------------------------------
inline double LinearFixedBoundarySurfaceMapper::GetValue(int i, int j) const
{
  return _Values->GetComponent(static_cast<vtkIdType>(i), j);
}

// -----------------------------------------------------------------------------
inline double LinearFixedBoundarySurfaceMapper::GetFreeValue(int i, int j) const
{
  return GetValue(FreePointId(i), j);
}

// -----------------------------------------------------------------------------
inline double LinearFixedBoundarySurfaceMapper::GetFixedValue(int i, int j) const
{
  return GetValue(FixedPointId(i), j);
}

// -----------------------------------------------------------------------------
inline void LinearFixedBoundarySurfaceMapper::SetValue(int i, double v)
{
  _Values->SetComponent(static_cast<vtkIdType>(i), 0, v);
}

// -----------------------------------------------------------------------------
inline void LinearFixedBoundarySurfaceMapper::SetValue(int i, int j, double v)
{
  _Values->SetComponent(static_cast<vtkIdType>(i), j, v);
}

// -----------------------------------------------------------------------------
inline void LinearFixedBoundarySurfaceMapper::SetFreeValue(int i, double v)
{
  SetValue(FreePointId(i), v);
}

// -----------------------------------------------------------------------------
inline void LinearFixedBoundarySurfaceMapper::SetFreeValue(int i, int j, double v)
{
  SetValue(FreePointId(i), j, v);
}


} // namespace mirtk

#endif // MIRTK_LinearFixedBoundarySurfaceMapper_H
