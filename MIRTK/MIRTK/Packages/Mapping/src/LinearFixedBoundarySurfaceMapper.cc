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

#include "mirtk/LinearFixedBoundarySurfaceMapper.h"

#include "mirtk/Memory.h"
#include "mirtk/PiecewiseLinearMap.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void LinearFixedBoundarySurfaceMapper
::CopyAttributes(const LinearFixedBoundarySurfaceMapper &other)
{
  _NumberOfIterations = other._NumberOfIterations;
  _Tolerance          = other._Tolerance;
  _PointIndex         = other._PointIndex;
  _FreePoints         = other._FreePoints;
  _FixedPoints        = other._FixedPoints;

  if (other._Values) {
    _Values.TakeReference(other._Values->NewInstance());
    _Values->DeepCopy(other._Values);
  } else {
    _Values = nullptr;
  }
}

// -----------------------------------------------------------------------------
LinearFixedBoundarySurfaceMapper::LinearFixedBoundarySurfaceMapper()
:
  _NumberOfIterations(-1),
  _Tolerance(-1.0)
{
}

// -----------------------------------------------------------------------------
LinearFixedBoundarySurfaceMapper
::LinearFixedBoundarySurfaceMapper(const LinearFixedBoundarySurfaceMapper &other)
:
  FixedBoundarySurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
LinearFixedBoundarySurfaceMapper &LinearFixedBoundarySurfaceMapper
::operator =(const LinearFixedBoundarySurfaceMapper &other)
{
  if (this != &other) {
    FixedBoundarySurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
LinearFixedBoundarySurfaceMapper::~LinearFixedBoundarySurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void LinearFixedBoundarySurfaceMapper::Initialize()
{
  // Initialize base class
  FixedBoundarySurfaceMapper::Initialize();

  // Check input
  if (_Surface->GetNumberOfPolys() == 0) {
    cerr << this->NameOfType() << "::Initialize: Input point set must be a surface mesh" << endl;
    exit(1);
  }

  // Initialize map values and
  // determine sets of points with free or fixed map values, respectively
  const int num = this->NumberOfPoints();
  const int dim = this->NumberOfComponents();

  #if MIRTK_USE_FLOAT_BY_DEFAULT
    _Values = vtkSmartPointer<vtkFloatArray>::New();
  #else
    _Values = vtkSmartPointer<vtkDoubleArray>::New();
  #endif
  _Values->SetName("SurfaceMap");
  _Values->SetNumberOfComponents(dim);
  _Values->SetNumberOfTuples(static_cast<vtkIdType>(num));

  _FreePoints .reserve(num);
  _FixedPoints.reserve(num);
  _PointIndex .resize (num);

  double p[3], *v = reinterpret_cast<double *>(_Values->GetVoidPointer(0));
  for (int i = 0; i < num; ++i, v += dim) {
    _Surface->GetPoint(static_cast<vtkIdType>(i), p);
    if (_Input->Evaluate(v, p)) {
      _PointIndex[i] = -(static_cast<int>(_FixedPoints.size()) + 1);
      _FixedPoints.push_back(i);
    } else {
      _PointIndex[i] = static_cast<int>(_FreePoints.size());
      _FreePoints.push_back(i);
    }
  }

  _FixedPoints.shrink_to_fit();
  _FreePoints .shrink_to_fit();
}

// -----------------------------------------------------------------------------
void LinearFixedBoundarySurfaceMapper::Finalize()
{
  // Assemble surface map
  SharedPtr<PiecewiseLinearMap> map = NewShared<PiecewiseLinearMap>();
  vtkSmartPointer<vtkPolyData> domain;
  domain.TakeReference(_Surface->NewInstance());
  domain->ShallowCopy(_Surface);
  domain->GetPointData()->Initialize();
  domain->GetCellData()->Initialize();
  map->Domain(domain);
  map->Values(_Values);
  _Output = map;

  // Finalize base class
  FixedBoundarySurfaceMapper::Finalize();
}


} // namespace mirtk
