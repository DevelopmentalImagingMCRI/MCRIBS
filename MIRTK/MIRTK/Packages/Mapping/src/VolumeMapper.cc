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

#include "mirtk/VolumeMapper.h"

#include "mirtk/Memory.h"
#include "mirtk/PointSetUtils.h"

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
void VolumeMapper::CopyAttributes(const VolumeMapper &other)
{
  _InputSet    = other._InputSet;
  _InputMap    = other._InputMap;
  _Boundary    = other._Boundary;
  _BoundaryMap = other._BoundaryMap;

  if (other._Output) {
    _Output = SharedPtr<Mapping>(other._Output->NewCopy());
  } else {
    _Output = nullptr;
  }
}

// -----------------------------------------------------------------------------
VolumeMapper::VolumeMapper()
{
}

// -----------------------------------------------------------------------------
VolumeMapper::VolumeMapper(const VolumeMapper &other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
VolumeMapper &VolumeMapper::operator =(const VolumeMapper &other)
{
  if (this != &other) {
    Object::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
VolumeMapper::~VolumeMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void VolumeMapper::Run()
{
  this->Initialize();
  this->Solve();
  this->Finalize();
}

// -----------------------------------------------------------------------------
void VolumeMapper::Initialize()
{
  // Free previous output map
  _Output = nullptr;

  // Check input
  if (!_InputSet) {
    cerr << this->NameOfType() << "::Initialize: Missing input point set" << endl;
    exit(1);
  }
  if (_InputSet->GetNumberOfCells() == 0) {
    cerr << this->NameOfType() << "::Initialize: Input has no cells" << endl;
    exit(1);
  }
  if (!_InputMap) {
    cerr << this->NameOfType() << "::Initialize: Missing input map" << endl;
    exit(1);
  }
  if (_InputMap->GetNumberOfTuples() != _InputSet->GetNumberOfPoints()) {
    cerr << this->NameOfType() << "::Initialize: Invalid input map" << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
void VolumeMapper::InitializeBoundary(vtkPointSet *input, vtkDataArray *map)
{
  vtkSmartPointer<vtkPointSet> volume;
  volume.TakeReference(input->NewInstance());
  volume->ShallowCopy(input);
  volume->GetCellData ()->Initialize();
  volume->GetPointData()->Initialize();
  const int i  = volume->GetPointData()->AddArray(map);
  _Boundary    = DataSetSurface(volume, true);
  _BoundaryMap = _Boundary->GetPointData()->GetArray(i);
  _BoundaryMap->SetName("BoundaryMap");
}

// -----------------------------------------------------------------------------
void VolumeMapper::Finalize()
{
}


} // namespace mirtk
