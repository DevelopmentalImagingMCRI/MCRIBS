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

#include "mirtk/FreeBoundarySurfaceMapper.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void FreeBoundarySurfaceMapper
::CopyAttributes(const FreeBoundarySurfaceMapper &other)
{
}

// -----------------------------------------------------------------------------
FreeBoundarySurfaceMapper::FreeBoundarySurfaceMapper()
{
}

// -----------------------------------------------------------------------------
FreeBoundarySurfaceMapper
::FreeBoundarySurfaceMapper(const FreeBoundarySurfaceMapper &other)
:
  SurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
FreeBoundarySurfaceMapper &FreeBoundarySurfaceMapper
::operator =(const FreeBoundarySurfaceMapper &other)
{
  if (this != &other) {
    SurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
FreeBoundarySurfaceMapper::~FreeBoundarySurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void FreeBoundarySurfaceMapper::Initialize()
{
  // Initialize base class
  SurfaceMapper::Initialize();

  // Input surface must have a boundary
  if (_Boundary->NumberOfSegments() == 0) {
    cerr << this->NameOfType() << "::Initialize: Surface mesh must have at least one boundary segment!" << endl;
    exit(1);
  }
}


} // namespace mirtk
