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

#include "mirtk/UniformBoundarySegmentParameterizer.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void UniformBoundarySegmentParameterizer::CopyAttributes(const UniformBoundarySegmentParameterizer &)
{
}

// -----------------------------------------------------------------------------
UniformBoundarySegmentParameterizer::UniformBoundarySegmentParameterizer()
{
}

// -----------------------------------------------------------------------------
UniformBoundarySegmentParameterizer
::UniformBoundarySegmentParameterizer(const UniformBoundarySegmentParameterizer &other)
:
  BoundarySegmentParameterizer(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
UniformBoundarySegmentParameterizer &UniformBoundarySegmentParameterizer
::operator =(const UniformBoundarySegmentParameterizer &other)
{
  if (this != &other) {
    BoundarySegmentParameterizer::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
UniformBoundarySegmentParameterizer::~UniformBoundarySegmentParameterizer()
{
}

// -----------------------------------------------------------------------------
BoundarySegmentParameterizer *UniformBoundarySegmentParameterizer::NewCopy() const
{
  return new UniformBoundarySegmentParameterizer(*this);
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void UniformBoundarySegmentParameterizer::Parameterize()
{
  double       t  = 0.;
  const double dt = 1. / _Boundary.NumberOfPoints();
  for (auto && v : _Values) {
    v = t, t += dt;
  }
}


} // namespace mirtk
