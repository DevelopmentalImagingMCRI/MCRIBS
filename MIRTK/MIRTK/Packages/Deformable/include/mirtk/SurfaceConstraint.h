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

#ifndef MIRTK_SurfaceConstraint_H
#define MIRTK_SurfaceConstraint_H

#include "mirtk/InternalForce.h"


namespace mirtk {


/**
 * Base class for a penalty term imposed on a point set surface
 *
 * Subclasses represent in particular internal forces of the boundary surface of
 * a simplicial complex such as curvature and non-self-intersection.
 */
class SurfaceConstraint : public InternalForce
{
  mirtkAbstractMacro(SurfaceConstraint);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Constructor
  SurfaceConstraint(const char * = "", double = 1.0);

  /// Copy constructor
  SurfaceConstraint(const SurfaceConstraint &);

  /// Assignment operator
  SurfaceConstraint &operator =(const SurfaceConstraint &);

public:

  /// Destructor
  virtual ~SurfaceConstraint();

};


} // namespace mirtk

#endif // MIRTK_SurfaceConstraint_H
