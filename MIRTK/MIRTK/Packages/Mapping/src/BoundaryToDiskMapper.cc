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

#include "mirtk/BoundaryToDiskMapper.h"

#include "mirtk/Math.h"
#include "mirtk/PointSetUtils.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void BoundaryToDiskMapper::CopyAttributes(const BoundaryToDiskMapper &other)
{
  _Radius = other._Radius;
}

// -----------------------------------------------------------------------------
BoundaryToDiskMapper::BoundaryToDiskMapper()
:
  _Radius(.0)
{
}

// -----------------------------------------------------------------------------
BoundaryToDiskMapper::BoundaryToDiskMapper(const BoundaryToDiskMapper &other)
:
  BoundarySegmentMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
BoundaryToDiskMapper &BoundaryToDiskMapper::operator =(const BoundaryToDiskMapper &other)
{
  if (this != &other) {
    BoundarySegmentMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
BoundaryToDiskMapper::~BoundaryToDiskMapper()
{
}

// -----------------------------------------------------------------------------
BoundaryMapper *BoundaryToDiskMapper::NewCopy() const
{
  return new BoundaryToDiskMapper(*this);
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void BoundaryToDiskMapper::Initialize()
{
  // Initialize base class
  BoundarySegmentMapper::Initialize();

  // Set default radius
  if (_Radius <= .0) {
    _Radius = sqrt(Area(_Boundary->Surface()) / pi);
  }
}

// -----------------------------------------------------------------------------
void BoundaryToDiskMapper
::MapBoundarySegment(int, const Array<int>    &indices,
                          const Array<double> &tvalues,
                          const Array<int>    &/*selection*/)
{
  const size_t npoints = indices.size();
  for (size_t n = 0; n < npoints; ++n) {
    const auto  &i = indices[n];
    const auto  &t = tvalues[n];
    const double a = t * two_pi;
    SetBoundaryValue(i, 0, _Radius * cos(a));
    SetBoundaryValue(i, 1, _Radius * sin(a));
  }
}


} // namespace mirtk
