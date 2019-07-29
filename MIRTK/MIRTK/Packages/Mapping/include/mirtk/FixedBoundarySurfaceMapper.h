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

#ifndef MIRTK_FixedBoundarySurfaceMapper_H
#define MIRTK_FixedBoundarySurfaceMapper_H

#include "mirtk/SurfaceMapper.h"

#include "mirtk/PiecewiseLinearMap.h"


namespace mirtk {


/**
 * Base class of solvers computing a surface map given a fixed boundary map
 *
 * Solvers of this type compute a surface map value for each interior point of
 * the non-closed input surface mesh given a fixed map value for each point on
 * at least one of the boundary segments. Additional either weak or hard
 * constraints at interior points may also be considered.
 *
 * \todo Mesh holes formed by those boundaries for which no Dirichlet boundary
 *       conditions are given as suggest in Marchandise et al. (2014)
 *       Optimal parametrizations for surface remeshing. See Figure 9.
 */
class FixedBoundarySurfaceMapper : public SurfaceMapper
{
  mirtkAbstractMacro(FixedBoundarySurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Surface boundary map
  ///
  /// This discrete mapping must map at least one of the closed surface boundary
  /// segments and possibly map also individual interior surface points. Whether
  /// these interior map values are considered as weak constraints, hard constraints,
  /// or not at all by a surface mapping method depends on the particular subclass
  /// implementation and respective mapping method. Only map values at boundary
  /// points are guaranteed to be enforced by subclasses.
  mirtkPublicAttributeMacro(SharedPtr<PiecewiseLinearMap>, Input);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const FixedBoundarySurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  FixedBoundarySurfaceMapper();

  /// Copy constructor
  FixedBoundarySurfaceMapper(const FixedBoundarySurfaceMapper &);

  /// Assignment operator
  FixedBoundarySurfaceMapper &operator =(const FixedBoundarySurfaceMapper &);

public:

  /// Destructor
  virtual ~FixedBoundarySurfaceMapper();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();


  // ---------------------------------------------------------------------------
  // Auxiliaries

protected:

  /// Number of boundary/surface map components
  int NumberOfComponents() const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline int FixedBoundarySurfaceMapper::NumberOfComponents() const
{
  return _Input->NumberOfComponents();
}


} // namespace mirtk

#endif // MIRTK_FixedBoundarySurfaceMapper_H
