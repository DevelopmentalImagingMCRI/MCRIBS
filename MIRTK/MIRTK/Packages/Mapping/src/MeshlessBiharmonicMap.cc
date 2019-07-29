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

#include "mirtk/MeshlessBiharmonicMap.h"

#include "mirtk/Point.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
MeshlessBiharmonicMap::MeshlessBiharmonicMap()
{
}

// -----------------------------------------------------------------------------
MeshlessBiharmonicMap::MeshlessBiharmonicMap(const MeshlessBiharmonicMap &other)
:
  MeshlessHarmonicMap(other)
{
}

// -----------------------------------------------------------------------------
MeshlessBiharmonicMap &MeshlessBiharmonicMap::operator =(const MeshlessBiharmonicMap &other)
{
  if (this != &other) {
    MeshlessHarmonicMap::operator =(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
void MeshlessBiharmonicMap::Initialize()
{
  // Initialize base class (skip MeshlessHarmonicMap::Initialize)
  MeshlessMap::Initialize();

  // Check parameters
  if (_SourcePoints.Size() == 0) {
    cerr << this->NameOfType() << "::Initialize: Set of source points is empty" << endl;
    exit(1);
  }
  if (_Coefficients.Rows() == 0 || _Coefficients.Cols() == 0) {
    cerr << this->NameOfType() << "::Initialize: Coefficients not set" << endl;
    exit(1);
  }
  if (_Coefficients.Rows() != 2 * _SourcePoints.Size()) {
    cerr << this->NameOfType() << "::Initialize: Number of coefficients must be twice the number of source points" << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
Mapping *MeshlessBiharmonicMap::NewCopy() const
{
  return new MeshlessBiharmonicMap(*this);
}

// -----------------------------------------------------------------------------
MeshlessBiharmonicMap::~MeshlessBiharmonicMap()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
bool MeshlessBiharmonicMap::Evaluate(double *v, double x, double y, double z) const
{
  const int n = _SourcePoints.Size();

  Point p(x, y, z);
  double    d, h, b;

  for (int j = 0; j < _Coefficients.Cols(); ++j) {
    v[j] = .0;
  }

  for (int i = 0; i < n; ++i) {
    d = p.Distance(_SourcePoints(i));
    if (d < 1e-12) {
      for (int j = 0; j < _Coefficients.Cols(); ++j) {
        v[j] = _OutsideValue;
      }
      return false;
    }
    h = H(d);
    b = B(d);
    for (int j = 0; j < _Coefficients.Cols(); ++j) {
      v[j] += h * _Coefficients(i,     j);
      v[j] += b * _Coefficients(i + n, j);
    }
  }

  return true;
}

// -----------------------------------------------------------------------------
double MeshlessBiharmonicMap::Evaluate(double x, double y, double z, int l) const
{
  const int n = _SourcePoints.Size();

  Point p(x, y, z);
  double    d, v = .0;

  for (int i = 0; i < n; ++i) {
    d = p.Distance(_SourcePoints(i));
    if (d < 1e-12) return _OutsideValue;
    v += H(d) * _Coefficients(i,     l);
    v += B(d) * _Coefficients(i + n, l);
  }

  return v;
}


} // namespace mirtk
