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

#include "mirtk/SurfaceMapper.h"

#include "mirtk/Assert.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void SurfaceMapper::CopyAttributes(const SurfaceMapper &other)
{
  _Surface   = other._Surface;
  _Boundary  = other._Boundary;
  _EdgeTable = other._EdgeTable;

  if (other._Output) {
    _Output = SharedPtr<Mapping>(other._Output->NewCopy());
  } else {
    _Output = nullptr;
  }
}

// -----------------------------------------------------------------------------
SurfaceMapper::SurfaceMapper()
{
}

// -----------------------------------------------------------------------------
SurfaceMapper::SurfaceMapper(const SurfaceMapper &other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
SurfaceMapper &SurfaceMapper::operator =(const SurfaceMapper &other)
{
  if (this != &other) {
    Object::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
SurfaceMapper::~SurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void SurfaceMapper::Run()
{
  this->Initialize();
  this->ComputeMap();
  this->Finalize();
}

// -----------------------------------------------------------------------------
void SurfaceMapper::Initialize()
{
  // Free previous map
  _Output = nullptr;

  // Check input
  if (!_Surface) {
    cerr << this->NameOfType() << "::Initialize: Missing input surface" << endl;
    exit(1);
  }
  if (_Surface->GetNumberOfPolys() == 0 && _Surface->GetNumberOfLines() == 0) {
    cerr << this->NameOfType() << "::Initialize: Input point set must be a surface mesh or polygon" << endl;
    exit(1);
  }

  // Build links
  _Surface->BuildLinks();

  // Determine adjacencies
  if (!_EdgeTable) {
    _EdgeTable = SharedPtr<class EdgeTable>(new class EdgeTable(_Surface));
  }

  // Extract boundaries
  if (!_Boundary) {
    _Boundary = SharedPtr<SurfaceBoundary>(new SurfaceBoundary(_Surface, _EdgeTable));
  }
}

// -----------------------------------------------------------------------------
void SurfaceMapper::Finalize()
{
  // Check that subclass produced output map
  mirtkAssert(_Output != nullptr, "output surface map is not nullptr");
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
int SurfaceMapper::GetEdgeNeighborPoints(int i, int j, int &k, int &l) const
{
  k = l = -1;
  unsigned short ncells1, ncells2, ncells = 0;
  vtkIdType      *cells1, *cells2, npts, *pts, ptId;
  _Surface->GetPointCells(static_cast<vtkIdType>(i), ncells1, cells1);
  _Surface->GetPointCells(static_cast<vtkIdType>(j), ncells2, cells2);
  for (unsigned short idx1 = 0; idx1 < ncells1; ++idx1) {
    for (unsigned short idx2 = 0; idx2 < ncells2; ++idx2) {
      if (cells1[idx1] == cells2[idx2]) {
        ++ncells;
        if (ncells < 3) {
          _Surface->GetCellPoints(cells1[idx1], npts, pts);
          if (npts == 3) {
            if      (pts[0] != i && pts[0] != j) ptId = pts[0];
            else if (pts[1] != i && pts[1] != j) ptId = pts[1];
            else if (pts[2] != i && pts[2] != j) ptId = pts[2];
            switch (ncells) {
              case 1: k = ptId; break;
              case 2: l = ptId; break;
            }
          }
        }
      }
    }
  }
  return ncells;
}


} // namespace mirtk
