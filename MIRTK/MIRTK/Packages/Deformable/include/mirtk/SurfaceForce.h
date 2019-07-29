/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Andreas Schuh
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

#ifndef MIRTK_SurfaceForce_H
#define MIRTK_SurfaceForce_H

#include "mirtk/ExternalForce.h"


namespace mirtk {


/**
 * Base class for external surface force terms
 *
 * Subclasses implement in particular external forces for deformable surface
 * models such as inflation/balloon forces and intensity edge forces. In case
 * of a tetrahedral mesh, these forces only apply to the corresponding
 * triangular surface of the simplicial complex boundary.
 */
class SurfaceForce : public ExternalForce
{
  mirtkAbstractMacro(SurfaceForce);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  SurfaceForce(const char * = "", double = 1.0);

  /// Copy constructor
  SurfaceForce(const SurfaceForce &);

  /// Assignment operator
  SurfaceForce &operator =(const SurfaceForce &);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const SurfaceForce &);

public:

  /// Destructor
  virtual ~SurfaceForce();

  /// Get distance to closest intersection with specified ray
  ///
  /// \param[in] p Ray starting point.
  /// \param[in] e Ray direction.
  /// \param[in] l Length of ray. When negative, a ray of length abs(l) is
  ///              cast in the direction opposite to \p e.
  ///
  /// \attention This function is not thread-safe because the VTK cell locator
  ///            used by this function is not thread-safe.
  ///
  /// \returns Distance to closest intersection with specified ray,
  ///          clamped to maximum distance equal to length of ray \p l.
  double IntersectWithRay(const double p[3], const double e[3], double l = .0) const;

  /// Get self-distance value at given world position along specified direction
  ///
  /// \param[in] p    Node position.
  /// \param[in] n    Node normal.
  /// \param[in] maxd Maximum distance, length of rays cast in opposite directions.
  ///
  /// \attention This function is not thread-safe because the VTK cell locator
  ///            used by this function is not thread-safe.
  ///
  /// \returns Distance to closest self-intersection along normal direction
  ///          (two rays cast in opposing directions), clamped to maximum
  ///          distance equal to length of rays \p l.
  double SelfDistance(const double p[3], const double n[3], double maxd = .0) const;

};


} // namespace mirtk

#endif // MIRTK_SurfaceForce_H
