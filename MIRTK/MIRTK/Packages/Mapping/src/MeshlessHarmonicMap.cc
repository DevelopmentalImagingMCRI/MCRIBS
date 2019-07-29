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

#include "mirtk/MeshlessHarmonicMap.h"

#include "mirtk/Point.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
MeshlessHarmonicMap::MeshlessHarmonicMap()
{
}

// -----------------------------------------------------------------------------
MeshlessHarmonicMap::MeshlessHarmonicMap(const MeshlessHarmonicMap &other)
:
  MeshlessMap(other)
{
}

// -----------------------------------------------------------------------------
MeshlessHarmonicMap &MeshlessHarmonicMap::operator =(const MeshlessHarmonicMap &other)
{
  if (this != &other) {
    MeshlessMap::operator =(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
Mapping *MeshlessHarmonicMap::NewCopy() const
{
  return new MeshlessHarmonicMap(*this);
}

// -----------------------------------------------------------------------------
MeshlessHarmonicMap::~MeshlessHarmonicMap()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
bool MeshlessHarmonicMap::Evaluate(double *v, double x, double y, double z) const
{
  Point  p(x, y, z);
  double d, h;

  for (int j = 0; j < _Coefficients.Cols(); ++j) {
    v[j] = .0;
  }

  for (int i = 0; i < _SourcePoints.Size(); ++i) {
    d = p.Distance(_SourcePoints(i));
    if (d < 1e-12) {
      for (int j = 0; j < _Coefficients.Cols(); ++j) {
        v[j] = _OutsideValue;
      }
      return false;
    }
    h = H(d);
    for (int j = 0; j < _Coefficients.Cols(); ++j) {
      v[j] += h * _Coefficients(i, j);
    }
  }

  return true;
}

// -----------------------------------------------------------------------------
double MeshlessHarmonicMap::Evaluate(double x, double y, double z, int l) const
{
  Point p(x, y, z);
  double    d, v = .0;

  for (int i = 0; i < _SourcePoints.Size(); ++i) {
    d = p.Distance(_SourcePoints(i));
    if (d < 1e-12) return _OutsideValue;
    v += H(d) * _Coefficients(i, l);
  }

  return v;
}


} // namespace mirtk
