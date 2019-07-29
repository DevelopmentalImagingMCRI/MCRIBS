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

#include "mirtk/BoundaryToPolygonMapper.h"

#include "mirtk/Math.h"
#include "mirtk/PointSetUtils.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void BoundaryToPolygonMapper::CopyAttributes(const BoundaryToPolygonMapper &other)
{
  _Radius = other._Radius;
}

// -----------------------------------------------------------------------------
BoundaryToPolygonMapper::BoundaryToPolygonMapper()
:
  _Radius(.0)
{
}

// -----------------------------------------------------------------------------
BoundaryToPolygonMapper::BoundaryToPolygonMapper(const BoundaryToPolygonMapper &other)
:
  BoundarySegmentMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
BoundaryToPolygonMapper &BoundaryToPolygonMapper::operator =(const BoundaryToPolygonMapper &other)
{
  if (this != &other) {
    BoundarySegmentMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
BoundaryToPolygonMapper::~BoundaryToPolygonMapper()
{
}

// -----------------------------------------------------------------------------
BoundaryMapper *BoundaryToPolygonMapper::NewCopy() const
{
  return new BoundaryToPolygonMapper(*this);
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void BoundaryToPolygonMapper::Initialize()
{
  // Initialize base class
  BoundarySegmentMapper::Initialize();

  // Set default radius
  if (_Radius <= .0) {
    _Radius = sqrt(Area(_Boundary->Surface()) / pi);
  }
}

// -----------------------------------------------------------------------------
void BoundaryToPolygonMapper
::MapBoundarySegment(int, const Array<int>    &indices,
                          const Array<double> &tvalues,
                          const Array<int>    &selection)
{
  const int npoints = static_cast<int>(indices.size());

  const Array<int> &corners = selection;
  if (corners.empty()) {
    // TODO: Auto-select corner points
  }
  if (corners.size() < 3) {
    cerr << this->NameOfType() << "::MapBoundarySegment: Select at least 3 boundary points as corners of the polygon!" << endl;
    exit(1);
  }

  // Map corner points to circumcircle
  for (size_t c = 0; c < corners.size(); ++c) {
    const auto  &n = corners[c];
    const auto  &i = indices[n];
    const auto  &t = tvalues[n];
    const double a = t * two_pi;
    SetBoundaryValue(i, 0, _Radius * sin(a));
    SetBoundaryValue(i, 1, _Radius * cos(a));
  }

  // Map intermediate points to the polygon edges
  double x1, x2, y1, y2, mx, my;
  for (size_t c = 0; c < corners.size(); ++c) {
    const auto &n1 = corners[c];
    const auto &n2 = corners[(c+1) % corners.size()];
    const auto &i1 = indices[n1];
    const auto &i2 = indices[n2];
    const auto &t1 = tvalues[n1];
    const auto &t2 = tvalues[n2];
    x1 = GetBoundaryValue(i1, 0);
    y1 = GetBoundaryValue(i1, 1);
    x2 = GetBoundaryValue(i2, 0);
    y2 = GetBoundaryValue(i2, 1);
    mx = (x2 - x1) / (t2 - t1);
    my = (y2 - y1) / (t2 - t1);
    if (n2 < n1) {
      for (int n = n1 + 1; n < npoints; ++n) {
        const auto &i = indices[n];
        const auto &t = tvalues[n];
        SetBoundaryValue(i, 0, x1 + mx * (t - t1));
        SetBoundaryValue(i, 1, y1 + my * (t - t1));
      }
      for (int n = 0; n < n2; ++n) {
        const auto &i = indices[n];
        const auto &t = tvalues[n];
        SetBoundaryValue(i, 0, x1 + mx * (t - t1));
        SetBoundaryValue(i, 1, y1 + my * (t - t1));
      }
    } else {
      for (int n = n1 + 1; n < n2; ++n) {
        const auto &i = indices[n];
        const auto &t = tvalues[n];
        SetBoundaryValue(i, 0, x1 + mx * (t - t1));
        SetBoundaryValue(i, 1, y1 + my * (t - t1));
      }
    }
  }
}


} // namespace mirtk
