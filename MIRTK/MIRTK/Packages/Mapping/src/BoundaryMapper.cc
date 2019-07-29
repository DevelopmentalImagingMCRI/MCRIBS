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

#include "mirtk/BoundaryMapper.h"

#include "mirtk/Algorithm.h"
#include "mirtk/Memory.h"
#include "mirtk/UnorderedMap.h"

#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void BoundaryMapper::CopyAttributes(const BoundaryMapper &other)
{
  _Boundary = other._Boundary;
  _Values   = other._Values;

  if (other._Output) {
    PiecewiseLinearMap *output;
    output = reinterpret_cast<PiecewiseLinearMap *>(other._Output->NewCopy());
    _Output = SharedPtr<PiecewiseLinearMap>(output);
  } else {
    _Output = nullptr;
  }
}

// -----------------------------------------------------------------------------
BoundaryMapper::BoundaryMapper()
{
}

// -----------------------------------------------------------------------------
BoundaryMapper::BoundaryMapper(const BoundaryMapper &other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
BoundaryMapper &BoundaryMapper::operator =(const BoundaryMapper &other)
{
  if (this != &other) {
    Object::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
BoundaryMapper::~BoundaryMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void BoundaryMapper::Run()
{
  this->Initialize();
  this->ComputeMap();
  this->Finalize();
}

// -----------------------------------------------------------------------------
int BoundaryMapper::NumberOfComponents() const
{
  return 2;
}

// -----------------------------------------------------------------------------
void BoundaryMapper::Initialize()
{
  // Free previous boundary map
  _Output = nullptr;

  // Check surface boundary
  if (!_Boundary) {
    cerr << this->NameOfType() << "::Initialize: No surface boundary set!" << endl;
    exit(1);
  }

  // Initialize boundary values
  const int num = _Boundary->NumberOfPoints();
  const int dim = this->NumberOfComponents();

  if (num < 1) {
    cerr << this->NameOfType() << "::Initialize: Surface is closed and has no boundary points to map!" << endl;
    exit(1);
  }
  if (dim < 1) {
    cerr << this->NameOfType() << "::Initialize: Subclass must return positive number of map components!" << endl;
    exit(1);
  }

  _Values.Resize(dim, num);
  _Values = mirtk::nan;
}

// -----------------------------------------------------------------------------
void BoundaryMapper::Finalize()
{
  const int num = _Boundary->NumberOfPoints();
  const int dim = this->NumberOfComponents();

  vtkSmartPointer<vtkPoints>    points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> verts  = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkCellArray> lines  = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkDataArray> values = vtkSmartPointer<vtkDoubleArray>::New();

  // Get only those boundary points with valid map value
  values->SetName("BoundaryMap");
  values->SetNumberOfComponents(dim);
  points->Allocate(static_cast<vtkIdType>(num));
  values->Allocate(static_cast<vtkIdType>(dim * num));

  double p[3];
  UnorderedMap<int, vtkIdType> ptIds;
  ptIds.reserve(num);

  for (int i = 0; i < num;  ++i) {
    if (HasBoundaryValue(i)) {
      _Boundary->GetPoint(i, p);
      ptIds[i] = points->InsertNextPoint(p);
      values->InsertNextTuple(_Values.Col(i));
    }
  }
  if (points->GetNumberOfPoints() == 0) {
    cerr << this->NameOfType() << "::Finalize: No boundary map values have been assigned!" << endl;
    exit(1);
  }

  points->Squeeze();
  values->Squeeze();

  // Determine topology of boundary map domain
  int       curPt, prePt, nxtPt;
  vtkIdType curId, preId, line[2];

  verts->Allocate(verts->EstimateSize(static_cast<vtkIdType>(ptIds.size()), 1));
  lines->Allocate(lines->EstimateSize(static_cast<vtkIdType>(ptIds.size()), 2));

  for (int n = 0; n < _Boundary->NumberOfSegments(); ++n) {
    prePt = _Boundary->PointIndex(n, -1);
    preId = (ptIds.find(prePt) == ptIds.end() ? -1 : -2);
    for (int i = 0; i < _Boundary->NumberOfPoints(n); ++i) {
      curPt = _Boundary->PointIndex(n, i);
      auto curIt = ptIds.find(curPt);
      if (curIt != ptIds.end()) {
        curId = curIt->second;
        nxtPt = _Boundary->PointIndex(n, i+1);
        auto nxtIt = ptIds.find(nxtPt);
        if (nxtIt != ptIds.end()) {
          line[0] = curId;
          line[1] = nxtIt->second;
          lines->InsertNextCell(2, line);
        } else if (preId == -1) {
          verts->InsertNextCell(1, &curId);
        }
      } else {
        curId = -1;
      }
      prePt = curPt;
      preId = curId;
    }
  }

  // Assemble piecewise linear output map
  vtkSmartPointer<vtkPolyData> domain = vtkSmartPointer<vtkPolyData>::New();
  domain->SetPoints(points);
  if (verts->GetNumberOfCells() > 0) {
    verts ->Squeeze();
    domain->SetVerts(verts);
  }
  if (lines->GetNumberOfCells() > 0) {
    lines ->Squeeze();
    domain->SetLines(lines);
  }

  _Output = NewShared<PiecewiseLinearMap>();
  _Output->Domain(domain);
  _Output->Values(values);
}


} // namespace mirtk
