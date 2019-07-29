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

#include "mirtk/MeshlessMap.h"

#include "mirtk/Point.h"


namespace mirtk {

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void MeshlessMap::CopyAttributes(const MeshlessMap &other)
{
  _SourcePoints = other._SourcePoints;
  _Coefficients = other._Coefficients;
}

// -----------------------------------------------------------------------------
MeshlessMap::MeshlessMap()
{
}

// -----------------------------------------------------------------------------
MeshlessMap::MeshlessMap(const MeshlessMap &other)
:
  Mapping(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
MeshlessMap &MeshlessMap::operator =(const MeshlessMap &other)
{
  if (this != &other) {
    Mapping::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
void MeshlessMap::Initialize()
{
  // Initialize base class
  Mapping::Initialize();

  // Check parameters
  if (_SourcePoints.Size() == 0) {
    cerr << this->NameOfType() << "::Initialize: Set of source points is empty" << endl;
    exit(1);
  }
  if (_Coefficients.Rows() == 0 || _Coefficients.Cols() == 0) {
    cerr << this->NameOfType() << "::Initialize: Coefficients not set" << endl;
    exit(1);
  }
  if (_Coefficients.Rows() != _SourcePoints.Size()) {
    cerr << this->NameOfType() << "::Initialize: Number of coefficients must be equal number of source points" << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
MeshlessMap::~MeshlessMap()
{
}

// =============================================================================
// Input domain
// =============================================================================

// -----------------------------------------------------------------------------
void MeshlessMap::BoundingBox(double &x1, double &y1, double &z1,
                              double &x2, double &y2, double &z2) const
{
  if (_SourcePoints.Size() > 0) {
    const Point &p0 = _SourcePoints(0);
    x1 = x2 = p0._x;
    y1 = y2 = p0._y;
    z1 = z1 = p0._z;
    for (int i = 1; i < _SourcePoints.Size(); ++i) {
      const Point &p = _SourcePoints(i);
      if (p._x < x1) x1 = p._x;
      if (p._x > x2) x2 = p._x;
      if (p._y < y1) y1 = p._y;
      if (p._y > y2) y2 = p._y;
      if (p._z < z1) z1 = p._z;
      if (p._z > z2) z2 = p._z;
    }
  } else {
    x1 = x2 = y1 = y2 = z1 = z2 = .0;
  }
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void MeshlessMap::ReadMap(Cifstream &is)
{
  _SourcePoints.Clear();
  is >> _SourcePoints >> _Coefficients;
}

// -----------------------------------------------------------------------------
void MeshlessMap::WriteMap(Cofstream &os) const
{
  os << _SourcePoints << _Coefficients;
}


} // namespace mirtk
