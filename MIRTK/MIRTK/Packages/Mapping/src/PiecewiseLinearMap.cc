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

#include "mirtk/PiecewiseLinearMap.h"

#include "mirtk/Vtk.h"
#include "mirtk/Path.h"
#include "mirtk/PointSetIO.h"
#include "mirtk/PointSetUtils.h"

#include "vtkPoints.h"
#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkGenericCell.h"
#include "vtkIdList.h"
#include "vtkCellLocator.h"

#include "vtkXMLImageDataWriter.h"
#include "vtkXMLImageDataReader.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

const double PiecewiseLinearMap::_Tolerance2 = 1e-9;

// -----------------------------------------------------------------------------
void PiecewiseLinearMap::CopyAttributes(const PiecewiseLinearMap &other)
{
  _Domain      = other._Domain;
  _Values      = other._Values;
  _MaxCellSize = other._MaxCellSize;

  if (other._Locator) {
    _Locator = vtkSmartPointer<vtkCellLocator>::New();
    _Locator->SetDataSet(_Domain);
    _Locator->BuildLocator();
  } else {
    _Locator = nullptr;
  }
}

// -----------------------------------------------------------------------------
PiecewiseLinearMap::PiecewiseLinearMap()
:
  _MaxCellSize(0)
{
}

// -----------------------------------------------------------------------------
PiecewiseLinearMap::PiecewiseLinearMap(const PiecewiseLinearMap &other)
:
  Mapping(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
PiecewiseLinearMap &PiecewiseLinearMap::operator =(const PiecewiseLinearMap &other)
{
  if (this != &other) {
    Mapping::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
void PiecewiseLinearMap::Initialize()
{
  // Check input mesh
  if (!_Domain) {
    cerr << this->NameOfType() << "::Initialize: No domain mesh is set" << endl;
    exit(1);
  }
  if (_Domain->GetNumberOfCells() == 0) {
    cerr << this->NameOfType() << "::Initialize: Domain mesh has no cells" << endl;
    exit(1);
  }
  if (!_Values) {
    cerr << this->NameOfType() << "::Initialize: No discrete map values at domain mesh points is set!" << endl;
    exit(1);
  }
  // Determine maximum number of cell points
  _MaxCellSize = _Domain->GetMaxCellSize();

  // Build cell locator
  _Locator = vtkSmartPointer<vtkCellLocator>::New();
  _Locator->SetDataSet(_Domain);
  _Locator->BuildLocator();
}

// -----------------------------------------------------------------------------
Mapping *PiecewiseLinearMap::NewCopy() const
{
  return new PiecewiseLinearMap(*this);
}

// -----------------------------------------------------------------------------
PiecewiseLinearMap::~PiecewiseLinearMap()
{
}

// =============================================================================
// Input domain
// =============================================================================

// -----------------------------------------------------------------------------
void PiecewiseLinearMap::BoundingBox(double &x1, double &y1, double &z1,
                                     double &x2, double &y2, double &z2) const
{
  if (_Domain) {
    double bounds[6];
    _Domain->GetBounds(bounds);
    x1 = bounds[0], x2 = bounds[1];
    y1 = bounds[2], y2 = bounds[3];
    z1 = bounds[4], z2 = bounds[5];
  } else {
    x1 = x2 = y1 = y2 = z1 = z2 = .0;
  }
}

// =============================================================================
// Map codomain
// =============================================================================

// -----------------------------------------------------------------------------
int PiecewiseLinearMap::NumberOfComponents() const
{
  return static_cast<int>(_Values->GetNumberOfComponents());
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkDataSet> PiecewiseLinearMap::Codomain() const
{
  const int dim = this->NumberOfComponents();
  vtkPolyData *domain = vtkPolyData::SafeDownCast(_Domain);
  if (domain == nullptr || dim < 2 || dim > 3) return nullptr;
  double p[3] = {.0};
  vtkSmartPointer<vtkPoints> points;
  points.TakeReference(domain->GetPoints()->NewInstance());
  points->SetNumberOfPoints(domain->GetNumberOfPoints());
  for (vtkIdType ptId = 0; ptId < domain->GetNumberOfPoints(); ++ptId) {
    _Values->GetTuple(ptId, p);
    points->SetPoint(ptId, p);
  }
  vtkSmartPointer<vtkPolyData> codomain;
  codomain.TakeReference(domain->NewInstance());
  codomain->ShallowCopy(domain);
  codomain->SetPoints(points);
  codomain->GetPointData()->Initialize();
  codomain->GetCellData()->Initialize();
  return codomain;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
bool PiecewiseLinearMap::Evaluate(double *v, double x, double y, double z) const
{
  vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
  double p[3] = {x, y, z}, pcoords[3];
  double *weight = new double[_MaxCellSize];
  if (_Locator->FindCell(p, _Tolerance2, cell, pcoords, weight) == -1) {
      for (int j = 0; j < _Values->GetNumberOfComponents(); ++j) {
        v[j] = _OutsideValue;
      }
      delete[] weight;
      return false;
  }
  for (int j = 0; j < _Values->GetNumberOfComponents(); ++j) {
    v[j] = .0;
  }
  vtkIdList * const ptIds = cell->GetPointIds();
  for (vtkIdType i = 0; i < ptIds->GetNumberOfIds(); ++i) {
    for (int j = 0; j < _Values->GetNumberOfComponents(); ++j) {
      v[j] += weight[i] * _Values->GetComponent(ptIds->GetId(i), j);
    }
  }
  delete[] weight;
  return true;
}

// -----------------------------------------------------------------------------
double PiecewiseLinearMap::Evaluate(double x, double y, double z, int l) const
{
  vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
  double p[3] = {x, y, z}, pcoords[3];
  double *weight = new double[_MaxCellSize];
  if (_Locator->FindCell(p, _Tolerance2, cell, pcoords, weight) == -1) {
      delete[] weight;
      return _OutsideValue;
  }
  double value = .0;
  vtkIdList * const ptIds = cell->GetPointIds();
  for (vtkIdType i = 0; i < ptIds->GetNumberOfIds(); ++i) {
    value += weight[i] * _Values->GetComponent(ptIds->GetId(i), l);
  }
  delete[] weight;
  return value;
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
bool PiecewiseLinearMap::Read(const char *fname)
{
  const string ext = Extension(fname);
  _Domain = nullptr;
  if (ext != ".vti") {
    const bool exit_on_failure = false;
    _Domain = ReadPointSet(fname, exit_on_failure);
  }
  if (_Domain == nullptr || _Domain->GetNumberOfPoints() == 0) {
    vtkNew<vtkXMLImageDataReader> reader;
    reader->SetFileName(fname);
    reader->Update();
    _Domain = reader->GetOutput();
  }
  if (_Domain->GetNumberOfPoints() == 0 ||
      _Domain->GetNumberOfCells()  == 0) return false;
  if (_Domain->GetPointData()->GetNumberOfArrays() == 1) {
    _Values = _Domain->GetPointData()->GetArray(0);
    _Domain->GetPointData()->Initialize();
  } else {
    _Values = nullptr;
  }
  this->Initialize();
  return true;
}

// -----------------------------------------------------------------------------
bool PiecewiseLinearMap::Write(const char *fname) const
{
  vtkSmartPointer<vtkDataSet> output;
  output.TakeReference(_Domain->NewInstance());
  output->ShallowCopy(_Domain);
  output->GetCellData ()->Initialize();
  output->GetPointData()->Initialize();
  if (_Values->GetNumberOfComponents() == 1) {
    output->GetPointData()->SetScalars(_Values);
  } else if (_Values->GetNumberOfComponents() == 2 ||
             _Values->GetNumberOfComponents() == 3) {
    output->GetPointData()->SetTCoords(_Values);
  } else {
    output->GetPointData()->AddArray(_Values);
  }
  vtkPointSet *pointset = vtkPointSet::SafeDownCast(output);
  if (pointset) {
    return WritePointSet(fname, pointset);
  }
  vtkNew<vtkXMLImageDataWriter> writer;
  writer->SetFileName(fname);
  SetVTKInput(writer, output);
  return writer->Write() == 0 ? true : false;
}


} // namespace mirtk
