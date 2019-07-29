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

#include "mirtk/ChordLengthSurfaceMapper.h"

#include "mirtk/Math.h"
#include "mirtk/VtkMath.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void ChordLengthSurfaceMapper::CopyAttributes(const ChordLengthSurfaceMapper &other)
{
  _Exponent = other._Exponent;
}

// -----------------------------------------------------------------------------
ChordLengthSurfaceMapper::ChordLengthSurfaceMapper(int p)
:
  _Exponent(p)
{
}

// -----------------------------------------------------------------------------
ChordLengthSurfaceMapper::ChordLengthSurfaceMapper(const ChordLengthSurfaceMapper &other)
:
  SymmetricWeightsSurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
ChordLengthSurfaceMapper &ChordLengthSurfaceMapper
::operator =(const ChordLengthSurfaceMapper &other)
{
  if (this != &other) {
    SymmetricWeightsSurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
ChordLengthSurfaceMapper::~ChordLengthSurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
double ChordLengthSurfaceMapper::Weight(int i, int j) const
{
  double p[3], q[3];
  _Surface->GetPoint(static_cast<vtkIdType>(i), p);
  _Surface->GetPoint(static_cast<vtkIdType>(j), q);
  const double dist2 = vtkMath::Distance2BetweenPoints(p, q);
  if (_Exponent <= 0) return 1.0;
  if (_Exponent == 1) return sqrt(dist2);
  if (_Exponent == 2) return dist2;
  if (_Exponent % 2 == 0) return pow(dist2, _Exponent / 2);
  return pow(sqrt(dist2), _Exponent);
}


} // namespace mirtk
