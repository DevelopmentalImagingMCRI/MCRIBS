/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
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

#include "mirtk/TetrahedralMeshMapper.h"

#include "mirtk/Vtk.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/PiecewiseLinearMap.h"

#include "vtkSmartPointer.h"
#include "vtkPointSet.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void TetrahedralMeshMapper::CopyAttributes(const TetrahedralMeshMapper &other)
{
  _InputMask = other._InputMask;
  if (other._Volume && other._Coords && other._BoundaryMask) {
    _Coords.TakeReference(other._Coords->NewInstance());
    _Coords->DeepCopy(other._Coords);
    _BoundaryMask.TakeReference(other._BoundaryMask->NewInstance());
    _BoundaryMask->DeepCopy(other._BoundaryMask);
    _Volume.TakeReference(other._Volume->NewInstance());
    _Volume->ShallowCopy(other._Volume);
    _Volume->GetPointData()->Initialize();
    _Volume->GetPointData()->AddArray(_Coords);
    _Volume->GetPointData()->AddArray(_BoundaryMask);
  } else {
    _Volume = NULL;
    _Coords = NULL;
    _BoundaryMask = NULL;
  }
  _NumberOfPoints         = other._NumberOfPoints;
  _NumberOfBoundaryPoints = other._NumberOfBoundaryPoints;
  _NumberOfInteriorPoints = other._NumberOfInteriorPoints;
}

// -----------------------------------------------------------------------------
TetrahedralMeshMapper::TetrahedralMeshMapper()
{
}

// -----------------------------------------------------------------------------
TetrahedralMeshMapper::TetrahedralMeshMapper(const TetrahedralMeshMapper &other)
:
  VolumeMapper(other),
  _NumberOfPoints(0),
  _NumberOfBoundaryPoints(0),
  _NumberOfInteriorPoints(0)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
TetrahedralMeshMapper &TetrahedralMeshMapper::operator =(const TetrahedralMeshMapper &other)
{
  if (this != &other) {
    VolumeMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
TetrahedralMeshMapper::~TetrahedralMeshMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void TetrahedralMeshMapper::Initialize()
{
  // Initialize base class
  VolumeMapper::Initialize();

  // Tetrahedralize interior of input point set
  int map_index, mask_index = -1;
  vtkSmartPointer<vtkPointSet> input;
  input.TakeReference(_InputSet->NewInstance());
  input->ShallowCopy(_InputSet);
  input->GetCellData ()->Initialize();
  input->GetPointData()->Initialize();
  map_index = input->GetPointData()->AddArray(_InputMap);
  if (_InputMask) mask_index = input->GetPointData()->AddArray(_InputMask);
  _Volume = Tetrahedralize(input);
  _Coords = _Volume->GetPointData()->GetArray(map_index);
  _Coords->SetName("VolumetricMap");

  // Extract surface of volume mesh
  this->InitializeBoundary(_Volume, _Coords);

  // Initialize boundary mask
  if (mask_index != -1) {
    _BoundaryMask = _Volume->GetPointData()->GetArray(mask_index);
  } else {
    _BoundaryMask = NewVtkDataArray(VTK_UNSIGNED_CHAR);
    _BoundaryMask->SetNumberOfComponents(1);
    _BoundaryMask->SetNumberOfTuples(_NumberOfPoints);
    _BoundaryMask->FillComponent(0, .0);
  }
  _BoundaryMask->SetName("BoundaryMask");
  vtkDataArray *origPtIds = _Boundary->GetPointData()->GetArray("vtkOriginalPointIds");
  for (vtkIdType ptId = 0, origPtId; ptId < _Boundary->GetNumberOfPoints(); ++ptId) {
    origPtId = static_cast<vtkIdType>(origPtIds->GetComponent(ptId, 0));
    _BoundaryMask->SetComponent(origPtId, 0, 1.0);
  }
  _NumberOfPoints         = static_cast<int>(_Volume->GetNumberOfPoints());
  _NumberOfBoundaryPoints = 0;
  for (vtkIdType ptId = 0; ptId < _BoundaryMask->GetNumberOfTuples(); ++ptId) {
    if (IsBoundaryPoint(ptId)) ++_NumberOfBoundaryPoints;
  }
  _NumberOfInteriorPoints = _NumberOfPoints - _NumberOfBoundaryPoints;
}

// -----------------------------------------------------------------------------
void TetrahedralMeshMapper::Finalize()
{
  // Create output map
  SharedPtr<PiecewiseLinearMap> map = NewShared<PiecewiseLinearMap>();
  map->Domain(_Volume);
  map->Values(_Coords);
  _Output = map;

  // Finalize base class
  VolumeMapper::Finalize();
}


} // namespace mirtk
