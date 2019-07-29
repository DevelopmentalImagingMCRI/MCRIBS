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

#ifndef MIRTK_BoundaryMapper_H
#define MIRTK_BoundaryMapper_H

#include "mirtk/Object.h"

#include "mirtk/Matrix.h"
#include "mirtk/Memory.h"
#include "mirtk/SurfaceBoundary.h"
#include "mirtk/PiecewiseLinearMap.h"


namespace mirtk {


/**
 * Base class of filters which assign map values to surface boundary points
 *
 * Filters of this type take a piecewice linear surface boundary segments
 * extracted from a surface mesh as input and compute a map value for each point
 * on the surface boundary. Such piecewise linear boundary map can be used as
 * fixed boundary condition for a surface map.
 *
 * \sa SurfaceMapper
 */
class BoundaryMapper : public Object
{
  mirtkAbstractMacro(BoundaryMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Extracted surface mesh boundary
  mirtkPublicAttributeMacro(SharedPtr<SurfaceBoundary>, Boundary);

  /// Map values at boundary points
  ///
  /// Not all boundary points need to have a boundary map value assigned by the
  /// subclass implementation. The map value at these points is Not-a-Number (NaN).
  mirtkAttributeMacro(Matrix, Values);

  /// Piecewise linear boundary map
  mirtkReadOnlyAttributeMacro(SharedPtr<PiecewiseLinearMap>, Output);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const BoundaryMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  BoundaryMapper();

  /// Copy constructor
  BoundaryMapper(const BoundaryMapper &);

  /// Assignment operator
  BoundaryMapper &operator =(const BoundaryMapper &);

public:

  /// Destructor
  virtual ~BoundaryMapper();

  /// Create new copy of this instance
  virtual BoundaryMapper *NewCopy() const = 0;

  // ---------------------------------------------------------------------------
  // Execution

  /// Compute boundary map values
  virtual void Run();

  /// Number of map value components
  virtual int NumberOfComponents() const;

  /// Whether a given boundary point has a mapped value
  ///
  /// \param[in] i Boundary point index.
  ///
  /// \returns Whether this mapper assigned a map value to the specified boundary point.
  bool HasBoundaryValue(int i) const;

  /// Get map value at boundary point
  ///
  /// \param[in] i Boundary point index.
  /// \param[in] j Index of map value component.
  ///
  /// \returns Map value component at surface boundary point.
  double GetBoundaryValue(int i, int j = 0) const;

  /// Get map value at surface point
  ///
  /// \param[in] ptId Surface point index.
  /// \param[in] j    Index of map value component.
  ///
  /// \returns Map value component at surface point.
  double GetSurfaceValue(int ptId, int j = 0) const;

protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Assign map values to boundary points of surface mesh
  virtual void ComputeMap() = 0;

  /// Finalize boundary map
  virtual void Finalize();

  /// Set map value at boundary point
  ///
  /// \param[in] i Boundary point index.
  /// \param[in] v Value of j-th map component.
  ///
  /// \returns Map value component at surface point.
  void SetBoundaryValue(int i, double v);

  /// Set map value at boundary point
  ///
  /// \param[in] i Boundary point index.
  /// \param[in] j Index of map value component.
  /// \param[in] v Value of j-th map component.
  ///
  /// \returns Map value component at surface point.
  void SetBoundaryValue(int i, int j, double v);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Boundary map
// =============================================================================

// -----------------------------------------------------------------------------
inline bool BoundaryMapper::HasBoundaryValue(int i) const
{
  const int dim = _Values.Rows();
  const double *value = _Values.Col(i);
  for (int j = 0; j < dim; ++j) {
    if (IsNaN(value[j])) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
inline void BoundaryMapper::SetBoundaryValue(int i, int j, double v)
{
  _Values.Col(i)[j] = v;
}

// -----------------------------------------------------------------------------
inline void BoundaryMapper::SetBoundaryValue(int i, double v)
{
  SetBoundaryValue(i, 0, v);
}

// -----------------------------------------------------------------------------
inline double BoundaryMapper::GetBoundaryValue(int i, int j) const
{
  return (HasBoundaryValue(i) ? _Values.Col(i)[j] : .0);
}

// -----------------------------------------------------------------------------
inline double BoundaryMapper::GetSurfaceValue(int ptId, int j) const
{
  const int i = _Boundary->Find(ptId);
  return (i < .0 ? .0 : GetBoundaryValue(i, j));
}


} // namespace mirtk

#endif // MIRTK_BoundaryMapper_H
