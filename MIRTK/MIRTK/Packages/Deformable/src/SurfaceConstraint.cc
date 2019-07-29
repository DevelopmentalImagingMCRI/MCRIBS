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

#include "mirtk/SurfaceConstraint.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
SurfaceConstraint::SurfaceConstraint(const char *name, double weight)
:
  InternalForce(name, weight)
{
  _SurfaceForce = true;
}

// -----------------------------------------------------------------------------
SurfaceConstraint::SurfaceConstraint(const SurfaceConstraint &other)
:
  InternalForce(other)
{
}

// -----------------------------------------------------------------------------
SurfaceConstraint &SurfaceConstraint::operator =(const SurfaceConstraint &other)
{
  InternalForce::operator =(other);
  return *this;
}

// -----------------------------------------------------------------------------
SurfaceConstraint::~SurfaceConstraint()
{
}


} // namespace mirtk
