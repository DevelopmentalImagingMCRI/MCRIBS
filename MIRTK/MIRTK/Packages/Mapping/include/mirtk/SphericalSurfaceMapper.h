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

#ifndef MIRTK_SphericalSurfaceMapper_H
#define MIRTK_SphericalSurfaceMapper_H

#include "mirtk/SurfaceMapper.h"


namespace mirtk {


/**
 * Base class of solvers computing a parameterization of a spherical surface
 *
 * Solvers of this type compute a surface parameterization of a closed surface
 * genus-0 mesh, i.e., a surface that is topologically equivalent to a sphere.
 */
class SphericalSurfaceMapper : public SurfaceMapper
{
  mirtkAbstractMacro(SphericalSurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Copy attributes of this class from another instance
  void CopyAttributes(const SphericalSurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  SphericalSurfaceMapper();

  /// Copy constructor
  SphericalSurfaceMapper(const SphericalSurfaceMapper &);

  /// Assignment operator
  SphericalSurfaceMapper &operator =(const SphericalSurfaceMapper &);

public:

  /// Destructor
  virtual ~SphericalSurfaceMapper();

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
inline int SphericalSurfaceMapper::NumberOfComponents() const
{
  return 2;
}


} // namespace mirtk

#endif // MIRTK_SphericalSurfaceMapper_H
