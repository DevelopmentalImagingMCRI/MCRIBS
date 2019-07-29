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

#include "mirtk/FixedBoundarySurfaceMapper.h"

#include "mirtk/BoundaryToSquareMapper.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void FixedBoundarySurfaceMapper
::CopyAttributes(const FixedBoundarySurfaceMapper &other)
{
  _Input = other._Input;
}

// -----------------------------------------------------------------------------
FixedBoundarySurfaceMapper::FixedBoundarySurfaceMapper()
{
}

// -----------------------------------------------------------------------------
FixedBoundarySurfaceMapper
::FixedBoundarySurfaceMapper(const FixedBoundarySurfaceMapper &other)
:
  SurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
FixedBoundarySurfaceMapper &FixedBoundarySurfaceMapper
::operator =(const FixedBoundarySurfaceMapper &other)
{
  if (this != &other) {
    SurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
FixedBoundarySurfaceMapper::~FixedBoundarySurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void FixedBoundarySurfaceMapper::Initialize()
{
  // Initialize base class
  SurfaceMapper::Initialize();

  // Input surface must have a boundary
  if (_Boundary->NumberOfSegments() == 0) {
    cerr << this->NameOfType() << "::Initialize: Surface mesh must have at least one boundary segment!" << endl;
    exit(1);
  }

  // Check boundary map
  if (!_Input) {
    cerr << this->NameOfType() << "::Initialize: No boundary map is set!" << endl;
    exit(1);
  }
  if (_Input->NumberOfPoints() == 0) {
    cerr << this->NameOfType() << "::Initialize: Invalid boundary map!" << endl;
    exit(1);
  }
}


} // namespace mirtk
