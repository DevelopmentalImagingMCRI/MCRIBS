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

#include "mirtk/SphericalSurfaceMapper.h"

#include "mirtk/PointSetUtils.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void SphericalSurfaceMapper::CopyAttributes(const SphericalSurfaceMapper &other)
{
}

// -----------------------------------------------------------------------------
SphericalSurfaceMapper::SphericalSurfaceMapper()
{
}

// -----------------------------------------------------------------------------
SphericalSurfaceMapper::SphericalSurfaceMapper(const SphericalSurfaceMapper &other)
:
  SurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
SphericalSurfaceMapper &SphericalSurfaceMapper::operator =(const SphericalSurfaceMapper &other)
{
  if (this != &other) {
    SurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
SphericalSurfaceMapper::~SphericalSurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void SphericalSurfaceMapper::Initialize()
{
  // Initialize base class
  SurfaceMapper::Initialize();

  // Input surface must not have a boundary
  if (_Boundary->NumberOfPoints() > 0) {
    cerr << this->NameOfType() << "::Initialize: Surface mesh must be closed without boundary!" << endl;
    exit(1);
  }

  // Check genus of surface mesh
  if (Genus(_Surface, *_EdgeTable) != 0.) {
    cerr << this->NameOfType() << "::Initialize: Surface mesh topology must be equivalent to a sphere!" << endl;
    exit(1);
  }
}


} // namespace mirtk
