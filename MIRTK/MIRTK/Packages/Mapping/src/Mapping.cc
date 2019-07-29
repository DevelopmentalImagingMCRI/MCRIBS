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

#include "mirtk/Mapping.h"

#include "mirtk/Config.h" // WINDOWS
#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Cfstream.h"
#include "mirtk/BaseImage.h"
#include "mirtk/VoxelFunction.h"
#include "mirtk/PointSetUtils.h"

#include "vtkSmartPointer.h"
#include "vtkPointSet.h"
#include "vtkPolyData.h"
#include "vtkImageData.h"


#include "mirtk/PiecewiseLinearMap.h"
#include "mirtk/MeshlessHarmonicMap.h"
#include "mirtk/MeshlessBiharmonicMap.h"


namespace mirtk {

// =============================================================================
// Factory method
// =============================================================================

// -----------------------------------------------------------------------------
Mapping *Mapping::New(const char *fname)
{
  UniquePtr<Mapping> map;

  const size_t max_name_len = 32;
  char         map_type_name[max_name_len];

  Cifstream is(fname);
  is.ReadAsChar(map_type_name, max_name_len);

  if (strncmp(map_type_name, MeshlessHarmonicMap::NameOfType(), max_name_len) == 0) {
    map.reset(new MeshlessHarmonicMap());
  } else if (strncmp(map_type_name, MeshlessBiharmonicMap::NameOfType(), max_name_len) == 0) {
    map.reset(new MeshlessBiharmonicMap());
  } else {
    // Note that a picewise linear map is stored using a VTK file format and
    // therefore the map file contains no "mirtk::PiecewiseLinearMap" header.
    map.reset(new PiecewiseLinearMap());
  }

  map->Read(fname);
  return map.release();
}

// =============================================================================
// Auxiliary functors
// =============================================================================

namespace MappingUtils {


// -----------------------------------------------------------------------------
/// Evaluate map at lattice points
class EvaluateMap : public VoxelFunction
{
  const Mapping   *_Map;
  vtkImageData    *_Domain;
  const BaseImage *_Output;
  const int        _NumberOfVoxels;
  const int        _l1, _l2;

public:

  EvaluateMap(const Mapping   *map,
              vtkImageData    *domain,
              const BaseImage *output,
              int l1, int l2)
  :
    _Map(map),
    _Domain(domain),
    _Output(output),
    _NumberOfVoxels(output->NumberOfSpatialVoxels()),
    _l1(l1), _l2(l2)
  {}

  template <class T>
  void operator ()(int i, int j, int k, int, T *v) const
  {
    if (!_Domain || _Domain->GetScalarComponentAsFloat(i, j, k, 0) != .0) {
      double x = i, y = j, z = k;
      _Output->ImageToWorld(x, y, z);
      double *f = new double[_Map->NumberOfComponents()];
      _Map->Evaluate(f, x, y, z);
      for (int l = _l1; l < _l2; ++l, v += _NumberOfVoxels) {
        *v = f[l];
      }
      delete[] f;
    } else {
      for (int l = _l1; l < _l2; ++l, v += _NumberOfVoxels) {
        *v = numeric_limits<T>::quiet_NaN();
      }
    }
  }
};


} // namespace MappingUtils
using namespace MappingUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void Mapping::CopyAttributes(const Mapping &other)
{
  _OutsideValue = other._OutsideValue;
}

// -----------------------------------------------------------------------------
Mapping::Mapping()
:
  _OutsideValue(mirtk::nan)
{
}

// -----------------------------------------------------------------------------
Mapping::Mapping(const Mapping &other)
:
  Object(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
Mapping &Mapping::operator =(const Mapping &other)
{
  if (this != &other) {
    Object::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
void Mapping::Initialize()
{
}

// -----------------------------------------------------------------------------
Mapping::~Mapping()
{
}

// =============================================================================
// Map domain
// =============================================================================

// -----------------------------------------------------------------------------
ImageAttributes Mapping::Attributes(int nx, int ny, int nz) const
{
  if (this->NumberOfArguments() >= 2) {
    if (ny <= 0) ny = nx;
  } else {
    ny = 0;
  }

  if (this->NumberOfArguments() >= 3) {
    if (nz <= 0) nz = nx;
  } else {
    nz = 0;
  }

  double x1, y1, z1, x2, y2, z2;
  this->BoundingBox(x1, y1, z1, x2, y2, z2);
  const double lx = x2 - x1;
  const double ly = y2 - y1;
  const double lz = z2 - z1;

  ImageAttributes lattice;
  lattice._xorigin = x1 + .5 * lx;
  lattice._yorigin = y1 + .5 * ly;
  lattice._zorigin = z1 + .5 * lz;
  lattice._x       = nx;
  lattice._y       = ny;
  lattice._z       = nz;
  lattice._dx      = (lattice._x == 0 ? .0 : lx / nx);
  lattice._dy      = (lattice._y == 0 ? .0 : ly / ny);
  lattice._dz      = (lattice._z == 0 ? .0 : lz / nz);
  if (lattice._x == 0) lattice._x = 1;
  if (lattice._y == 0) lattice._y = 1;
  if (lattice._z == 0) lattice._z = 1;

  return lattice;
}

// -----------------------------------------------------------------------------
ImageAttributes Mapping::Attributes(double dx, double dy, double dz) const
{
  double x1, y1, z1, x2, y2, z2;
  this->BoundingBox(x1, y1, z1, x2, y2, z2);
  const double lx = x2 - x1;
  const double ly = y2 - y1;
  const double lz = z2 - z1;

  if (dx <= .0) {
    dx = sqrt(lx*lx + ly*ly + lz*lz) / 256;
  }

  if (this->NumberOfArguments() >= 2) {
    if (dy <= .0) dy = dx;
  } else {
    dy = 0;
  }

  if (this->NumberOfArguments() >= 3) {
    if (dz <= .0) dz = dx;
  } else {
    dz = 0;
  }

  ImageAttributes lattice;
  lattice._xorigin = x1 + .5 * lx;
  lattice._yorigin = y1 + .5 * ly;
  lattice._zorigin = z1 + .5 * lz;
  lattice._x       = iceil(lx / dx);
  lattice._y       = iceil(ly / dy);
  lattice._z       = iceil(lz / dz);
  lattice._dx      = (lattice._x == 0 ? .0 : dx);
  lattice._dy      = (lattice._y == 0 ? .0 : dy);
  lattice._dz      = (lattice._z == 0 ? .0 : dz);
  if (lattice._x == 0) lattice._x = 1;
  if (lattice._y == 0) lattice._y = 1;
  if (lattice._z == 0) lattice._z = 1;

  return lattice;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void Mapping::Evaluate(GenericImage<float> &f, int l, vtkSmartPointer<vtkPointSet> m) const
{
  ImageAttributes lattice = f.Attributes();
  lattice._dt = .0;

  if (l >= NumberOfComponents() || l + lattice._t > NumberOfComponents()) {
    cerr << this->NameOfType() << "::Evaluate: Component index out of range" << endl;
    exit(1);
  }

  vtkSmartPointer<vtkImageData> mask;
  if (m) {
    mask = NewVtkMask(lattice._x, lattice._y, lattice._z);
    ImageStencilToMask(ImageStencil(mask, WorldToImage(m, &f)), mask);
  }

  ParallelForEachVoxel(EvaluateMap(this, mask, &f, l, l + lattice._t), lattice, f);
}

// -----------------------------------------------------------------------------
void Mapping::Evaluate(GenericImage<double> &f, int l, vtkSmartPointer<vtkPointSet> m) const
{
  ImageAttributes lattice = f.Attributes();
  lattice._dt = .0;

  if (l >= NumberOfComponents() || l + lattice._t > NumberOfComponents()) {
    cerr << this->NameOfType() << "::Evaluate: Component index out of range" << endl;
    exit(1);
  }

  vtkSmartPointer<vtkImageData> mask;
  if (m) {
    mask = NewVtkMask(lattice._x, lattice._y, lattice._z);
    ImageStencilToMask(ImageStencil(mask, WorldToImage(m, &f)), mask);
  }

  ParallelForEachVoxel(EvaluateMap(this, mask, &f, l, l + lattice._t), lattice, f);
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
bool Mapping::Read(const char *fname)
{
  Cifstream is(fname);
  const size_t max_name_len = 32;
  char         map_type_name[max_name_len];
  is.ReadAsChar(map_type_name, max_name_len);
  if (strncmp(map_type_name, this->NameOfClass(), max_name_len) != 0) {
    return false;
  }
  this->ReadMap(is);
  this->Initialize();
  return true;
}

// -----------------------------------------------------------------------------
bool Mapping::Write(const char *fname) const
{
  Cofstream os(fname);
  const size_t max_name_len = 32;
  char         map_type_name[max_name_len] = {0};
  #ifdef WINDOWS
    strncpy_s(map_type_name, max_name_len, this->NameOfClass(), max_name_len);
  #else
    strncpy(map_type_name, this->NameOfClass(), max_name_len);
  #endif
  os.WriteAsChar(map_type_name, max_name_len);
  this->WriteMap(os);
  return true;
}

// -----------------------------------------------------------------------------
void Mapping::ReadMap(Cifstream &)
{
  cerr << this->NameOfClass() << "::ReadMap not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void Mapping::WriteMap(Cofstream &) const
{
  cerr << this->NameOfClass() << "::WriteMap not implemented" << endl;
  exit(1);
}


} // namespace mirtk
