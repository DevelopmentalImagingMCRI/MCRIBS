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
#include "mirtk/Vtk.h"

#include "vtkIdList.h"


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
  vtkIdType ptId, ncells = 0;
  vtkNew<vtkIdList> cellIds1, cellIds2, ptIds;
  _Surface->GetPointCells(static_cast<vtkIdType>(i), cellIds1.GetPointer());
  _Surface->GetPointCells(static_cast<vtkIdType>(j), cellIds2.GetPointer());
  for (vtkIdType idx1 = 0; idx1 < cellIds1->GetNumberOfIds(); ++idx1)
  for (vtkIdType idx2 = 0; idx2 < cellIds2->GetNumberOfIds(); ++idx2) {
    if (cellIds1->GetId(idx1) == cellIds2->GetId(idx2)) {
      ++ncells;
      if (ncells < 3) {
        GetCellPoints(_Surface, cellIds1->GetId(idx1), ptIds.GetPointer());
        if (ptIds->GetNumberOfIds() == 3) {
          if      (ptIds->GetId(0) != i && ptIds->GetId(0) != j) ptId = ptIds->GetId(0);
          else if (ptIds->GetId(1) != i && ptIds->GetId(1) != j) ptId = ptIds->GetId(1);
          else if (ptIds->GetId(2) != i && ptIds->GetId(2) != j) ptId = ptIds->GetId(2);
          switch (ncells) {
            case 1: k = ptId; break;
            case 2: l = ptId; break;
          }
        }
      }
    }
  }
  return ncells;
}


} // namespace mirtk
