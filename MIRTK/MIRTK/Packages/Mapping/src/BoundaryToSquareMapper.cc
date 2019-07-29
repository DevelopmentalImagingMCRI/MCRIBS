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

#include "mirtk/BoundaryToSquareMapper.h"

#include "mirtk/Math.h"
#include "mirtk/PointSetUtils.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void BoundaryToSquareMapper::CopyAttributes(const BoundaryToSquareMapper &other)
{
  _SideLength = other._SideLength;
}

// -----------------------------------------------------------------------------
BoundaryToSquareMapper::BoundaryToSquareMapper()
:
  _SideLength(.0)
{
}

// -----------------------------------------------------------------------------
BoundaryToSquareMapper::BoundaryToSquareMapper(const BoundaryToSquareMapper &other)
:
  BoundarySegmentMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
BoundaryToSquareMapper &BoundaryToSquareMapper::operator =(const BoundaryToSquareMapper &other)
{
  if (this != &other) {
    BoundarySegmentMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
BoundaryToSquareMapper::~BoundaryToSquareMapper()
{
}

// -----------------------------------------------------------------------------
BoundaryMapper *BoundaryToSquareMapper::NewCopy() const
{
  return new BoundaryToSquareMapper(*this);
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void BoundaryToSquareMapper::Initialize()
{
  // Initialize base class
  BoundarySegmentMapper::Initialize();

  // Set default side length
  if (_SideLength <= .0) {
    _SideLength = sqrt(Area(_Boundary->Surface()));
  }
}

// -----------------------------------------------------------------------------
void BoundaryToSquareMapper
::MapBoundarySegment(int, const Array<int>    &indices,
                          const Array<double> &tvalues,
                          const Array<int>    &selection)
{
  const size_t npoints = indices.size();
  const double radius  = .5 * _SideLength;

  if (npoints < 4) {
    cerr << this->NameOfClass() << ":MapBoundarySegment: Boundary segment must have at least 4 points!" << endl;
    exit(1);
  }

  // Determine corner points
  size_t cornerIdx [4] = {0};
  for (size_t n = 0, s = 0; n < npoints; ++n) {
    const auto &t0 = tvalues[cornerIdx[s]];
    const auto &t  = tvalues[n];
    if (t - t0 >= .25) {
      ++s;
      if (n > cornerIdx[s-1] + 1 && .25 - tvalues[n-1] + t0 < t - t0 - .25) {
        cornerIdx[s] = n - 1;
      } else {
        cornerIdx[s] = n;
      }
    }
  }

  // Side 1
  double L = tvalues[cornerIdx[1]];
  for (size_t n = cornerIdx[0]; n < cornerIdx[1]; ++n) {
    const auto &i = indices[n];
    const auto &t = tvalues[n];
    SetBoundaryValue(i, 0, -radius + (t / L) * _SideLength);
    SetBoundaryValue(i, 1, +radius);
  }

  // Side 2
  L = (tvalues[cornerIdx[2]] - tvalues[cornerIdx[1]]);
  for (size_t n = cornerIdx[1]; n < cornerIdx[2]; ++n) {
    const auto &i = indices[n];
    const auto &t = tvalues[n];
    SetBoundaryValue(i, 0, +radius);
    SetBoundaryValue(i, 1, +radius - ((t - .25) / L) * _SideLength);
  }

  // Side 3
  L = (tvalues[cornerIdx[3]] - tvalues[cornerIdx[2]]);
  for (size_t n = cornerIdx[2]; n < cornerIdx[3]; ++n) {
    const auto &i = indices[n];
    const auto &t = tvalues[n];
    SetBoundaryValue(i, 0, +radius - ((t - .5) / L) * _SideLength);
    SetBoundaryValue(i, 1, -radius);
  }

  // Side 4
  L = (1.0 - tvalues[cornerIdx[3]]);
  for (size_t n = cornerIdx[3]; n < npoints; ++n) {
    const auto &i = indices[n];
    const auto &t = tvalues[n];
    SetBoundaryValue(i, 0, -radius);
    SetBoundaryValue(i, 1, -radius + ((t - .75) / L) * _SideLength);
  }
}


} // namespace mirtk
