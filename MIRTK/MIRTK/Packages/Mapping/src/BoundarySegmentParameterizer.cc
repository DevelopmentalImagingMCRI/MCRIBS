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

#include "mirtk/BoundarySegmentParameterizer.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void BoundarySegmentParameterizer::CopyAttributes(const BoundarySegmentParameterizer &other)
{
  _Boundary = other._Boundary;
  _Values   = other._Values;
}

// -----------------------------------------------------------------------------
BoundarySegmentParameterizer::BoundarySegmentParameterizer()
{
}

// -----------------------------------------------------------------------------
BoundarySegmentParameterizer
::BoundarySegmentParameterizer(const BoundarySegmentParameterizer &other)
:
  Object(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
BoundarySegmentParameterizer &BoundarySegmentParameterizer
::operator =(const BoundarySegmentParameterizer &other)
{
  if (this != &other) {
    Object::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
BoundarySegmentParameterizer::~BoundarySegmentParameterizer()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void BoundarySegmentParameterizer::Run()
{
  this->Initialize();
  this->Parameterize();
  this->Finalize();
}

// -----------------------------------------------------------------------------
void BoundarySegmentParameterizer::Initialize()
{
  // Check boundary segment
  const int npoints = _Boundary.NumberOfPoints();
  if (npoints < 2) {
    cerr << this->NameOfType() << "::Initialize: Boundary segment must have at least two points!" << endl;
    exit(1);
  }

  // Initialize parameter values
  _Values.clear();
  _Values.resize(npoints, -1.0);
}

// -----------------------------------------------------------------------------
void BoundarySegmentParameterizer::Finalize()
{
#ifndef NDEBUG
  for (auto t : _Values) {
    if (t < 0. || t >= 1.) {
      cerr << this->NameOfType() << "::Finalize: Curve parameter value must be in [0, 1), but t=" << t << endl;
      exit(1);
    }
  }
#endif

  // Reparameterize such that first selected point has t=0
  if (_Boundary.NumberOfSelectedPoints() > 0) {
    int    i0 = _Boundary.SelectedPointIndex(0);
    double t0 = _Values[i0];
    if (t0 != .0) {
      for (int i = 0; i < _Boundary.NumberOfPoints(); ++i) {
        _Values[i] -= t0;
        if (_Values[i] < .0) _Values[i] += 1.0;
      }
      _Values[i0] = .0;
    }
  }

  // Revert orientation of curve if points where selected in reverse order
  // or when only two points where selected and the second point is further
  // away from the first point at t=0 then if curve is being reversed
  if (_Boundary.NumberOfSelectedPoints() > 1) {
    double t1 = _Values[_Boundary.SelectedPointIndex(1)];
    double t2 = 1.0 - t1;
    if (_Boundary.NumberOfSelectedPoints() > 2) {
      t2 = _Values[_Boundary.SelectedPointIndex(2)];
    }
    if (t1 > t2) {
      for (int n = 0; n < _Boundary.NumberOfPoints(); ++n) {
        if (_Values[n] != .0) {
          _Values[n] = 1.0 - _Values[n];
        }
      }
    }
  }
}


} // namespace mirtk
