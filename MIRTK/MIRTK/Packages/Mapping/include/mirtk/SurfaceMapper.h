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

#ifndef MIRTK_SurfaceMapper_H
#define MIRTK_SurfaceMapper_H

#include "mirtk/Object.h"

#include "mirtk/Memory.h"
#include "mirtk/Point.h"

#include "mirtk/EdgeTable.h"
#include "mirtk/SurfaceBoundary.h"
#include "mirtk/Mapping.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"


namespace mirtk {


/**
 * Base class of solvers for the computation of a surface map
 *
 * Solvers of this type compute a map value for each point on the surface
 * of the input mesh embedded in 3D Euclidean space. The properties of the
 * output map depend on the boundary conditions and the specific solver used to
 * compute the surface map. Examples of surface maps are the assignment of 2D
 * texture coordinates and a bijective mapping of a surface mesh to a disk,
 * square, or sphere, respectively.
 */
class SurfaceMapper : public Object
{
  mirtkAbstractMacro(SurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Polygonal representation of surface map domain
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Surface);

  /// Pre-computed edge table of surface mesh
  mirtkPublicAttributeMacro(SharedPtr<mirtk::EdgeTable>, EdgeTable);

  /// Extracted surface mesh boundary
  mirtkPublicAttributeMacro(SharedPtr<SurfaceBoundary>, Boundary);

  /// Output surface map
  ///
  /// \note The output map is uninitialized! Mapping::Initialize must be
  ///       called before this map can be evaluated at map domain points.
  mirtkReadOnlyAttributeMacro(SharedPtr<Mapping>, Output);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const SurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  SurfaceMapper();

  /// Copy constructor
  SurfaceMapper(const SurfaceMapper &);

  /// Assignment operator
  SurfaceMapper &operator =(const SurfaceMapper &);

public:

  /// Destructor
  virtual ~SurfaceMapper();

  // ---------------------------------------------------------------------------
  // Execution

  /// Compute surface map
  void Run();

protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Compute surface map
  virtual void ComputeMap() = 0;

  /// Finalize filter execution
  virtual void Finalize();

  // ---------------------------------------------------------------------------
  // Auxiliaries

protected:

  /// Number of surface points
  int NumberOfPoints() const;

  /// Number of surface points minus the number of boundary points
  int NumberOfInteriorPoints() const;

  /// Get surface point coordinates
  ///
  /// \param[in]  ptId Surface point ID.
  /// \param[out] p    Point coordinates.
  void GetPoint(int ptId, double p[3]) const;

  /// Get surface point coordinates
  ///
  /// \param[in] ptId Surface point ID.
  ///
  /// \returns Point point.
  class Point Point(int ptId) const;

  /// Get IDs of surface points belong to the two triangles sharing an edge
  ///
  /// \param[in]  i Edge start point.
  /// \param[in]  j Edge end point.
  /// \param[out] k Other point of first adjacent triangle.
  /// \param[out] l Other point of second adjacent triangle. When the edge is
  ///               on the surface boundary and therefore has only as single
  ///               adjacent triangle, -1 is returned.
  ///
  /// \returns Number of cells adjacent to the specified edge. When the surface
  ///          mesh is triangulated, the return value is 1 for a boundary edge,
  ///          and 2 for an interior edge.
  int GetEdgeNeighborPoints(int i, int j, int &k, int &l) const;
};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
inline int SurfaceMapper::NumberOfPoints() const
{
  return static_cast<int>(_Surface->GetNumberOfPoints());
}

// -----------------------------------------------------------------------------
inline int SurfaceMapper::NumberOfInteriorPoints() const
{
  return NumberOfPoints() - _Boundary->NumberOfPoints();
}

// -----------------------------------------------------------------------------
inline void SurfaceMapper::GetPoint(int ptId, double p[3]) const
{
  _Surface->GetPoint(static_cast<vtkIdType>(ptId), p);
}

// -----------------------------------------------------------------------------
inline Point SurfaceMapper::Point(int ptId) const
{
  double p[3];
  GetPoint(ptId, p);
  return p;
}


} // namespace mirtk

#endif // MIRTK_SurfaceMapper_H
