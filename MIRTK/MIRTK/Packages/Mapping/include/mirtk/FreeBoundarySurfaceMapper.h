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

#ifndef MIRTK_FreeBoundarySurfaceMapper_H
#define MIRTK_FreeBoundarySurfaceMapper_H

#include "mirtk/SurfaceMapper.h"


namespace mirtk {


/**
 * Base class of solvers computing a surface map with free boundary values
 *
 * Solvers of this type compute a surface parameterization of a non-closed surface
 * mesh, where the boundary points are not assigned a predefined fixed boundary
 * map value. These solvers cannot be used to interpolate arbitrary values
 * given on the boundary at interior points. Instead, the result will always
 * be a 2D parameterization of the surface embedded in 3D Euclidean space.
 */
class FreeBoundarySurfaceMapper : public SurfaceMapper
{
  mirtkAbstractMacro(FreeBoundarySurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Copy attributes of this class from another instance
  void CopyAttributes(const FreeBoundarySurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  FreeBoundarySurfaceMapper();

  /// Copy constructor
  FreeBoundarySurfaceMapper(const FreeBoundarySurfaceMapper &);

  /// Assignment operator
  FreeBoundarySurfaceMapper &operator =(const FreeBoundarySurfaceMapper &);

public:

  /// Destructor
  virtual ~FreeBoundarySurfaceMapper();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  // ---------------------------------------------------------------------------
  // Auxiliaries

public:

  /// Dimension of map codomain
  int NumberOfComponents() const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline int FreeBoundarySurfaceMapper::NumberOfComponents() const
{
  return 2;
}


} // namespace mirtk

#endif // MIRTK_FreeBoundarySurfaceMapper_H
