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

#include "mirtk/ChordLengthBoundarySegmentParameterizer.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void ChordLengthBoundarySegmentParameterizer::CopyAttributes(const ChordLengthBoundarySegmentParameterizer &)
{
}

// -----------------------------------------------------------------------------
ChordLengthBoundarySegmentParameterizer::ChordLengthBoundarySegmentParameterizer()
{
}

// -----------------------------------------------------------------------------
ChordLengthBoundarySegmentParameterizer
::ChordLengthBoundarySegmentParameterizer(const ChordLengthBoundarySegmentParameterizer &other)
:
  BoundarySegmentParameterizer(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
ChordLengthBoundarySegmentParameterizer &ChordLengthBoundarySegmentParameterizer
::operator =(const ChordLengthBoundarySegmentParameterizer &other)
{
  if (this != &other) {
    BoundarySegmentParameterizer::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
ChordLengthBoundarySegmentParameterizer::~ChordLengthBoundarySegmentParameterizer()
{
}

// -----------------------------------------------------------------------------
BoundarySegmentParameterizer *ChordLengthBoundarySegmentParameterizer::NewCopy() const
{
  return new ChordLengthBoundarySegmentParameterizer(*this);
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void ChordLengthBoundarySegmentParameterizer::Parameterize()
{
  const Vector l = _Boundary.EdgeLengths();
  const double L = l.Sum();

  double t = .0;
  for (int i = 0; i < _Boundary.NumberOfPoints(); ++i) {
    _Values[i] = t;
    t += l(i) / L;
  }
}


} // namespace mirtk
