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

#include "mirtk/BoundarySegmentMapper.h"

#include "mirtk/Algorithm.h"
#include "mirtk/ChordLengthBoundarySegmentParameterizer.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void BoundarySegmentMapper::CopyAttributes(const BoundarySegmentMapper &other)
{
  _Parameterizer = SharedPtr<BoundarySegmentParameterizer>(other._Parameterizer->NewCopy());
}

// -----------------------------------------------------------------------------
BoundarySegmentMapper::BoundarySegmentMapper()
:
  _Parameterizer(NewShared<ChordLengthBoundarySegmentParameterizer>())
{
}

// -----------------------------------------------------------------------------
BoundarySegmentMapper::BoundarySegmentMapper(const BoundarySegmentMapper &other)
:
  BoundaryMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
BoundarySegmentMapper &BoundarySegmentMapper::operator =(const BoundarySegmentMapper &other)
{
  if (this != &other) {
    BoundaryMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
BoundarySegmentMapper::~BoundarySegmentMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void BoundarySegmentMapper::Initialize()
{
  // Initialize base class
  BoundaryMapper::Initialize();

  // Ensure that boundary segment parameterizer is set
  if (_Parameterizer == nullptr) {
    _Parameterizer = NewShared<ChordLengthBoundarySegmentParameterizer>();
  }
}

// -----------------------------------------------------------------------------
void BoundarySegmentMapper::ComputeMap()
{
  this->MapLongest();
}

// -----------------------------------------------------------------------------
void BoundarySegmentMapper::MapLongest()
{
  this->MapSegment(_Boundary->FindLongestSegment());
}

// -----------------------------------------------------------------------------
void BoundarySegmentMapper::MapLargest()
{
  this->MapSegment(_Boundary->FindLargestSegment());
}

// -----------------------------------------------------------------------------
void BoundarySegmentMapper::MapSegment(int n)
{
  // Get boundary segment
  const BoundarySegment &segment = _Boundary->Segment(n);

  // Parameterize boundary curve
  _Parameterizer->Boundary(segment);
  _Parameterizer->Run();

  // Sort points by increasing parameter value
  const int        npoints = segment.NumberOfPoints();
  const Array<int> indices = _Boundary->PointIndices(n);
  const Array<int> order   = IncreasingOrder(_Parameterizer->Values());

  Array<int>    i(npoints);
  Array<double> t(npoints);
  Array<int>    selection;
  selection.reserve(segment.NumberOfSelectedPoints());

  for (int j = 0; j < npoints; ++j) {
    const auto &pos = order[j];
    i[j] = indices[pos];
    t[j] = _Parameterizer->Value(pos);
    if (segment.IsSelected(pos)) {
      selection.push_back(j);
    }
  }

  // Assign map values to points of boundary segment
  this->MapBoundarySegment(n, i, t, selection);
}

// -----------------------------------------------------------------------------
void BoundarySegmentMapper::Finalize()
{
  // Finalize base class
  BoundaryMapper::Finalize();
}


} // namespace mirtk
