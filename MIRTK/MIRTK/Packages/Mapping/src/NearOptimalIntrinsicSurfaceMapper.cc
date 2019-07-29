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

#include "mirtk/NearOptimalIntrinsicSurfaceMapper.h"

#include "mirtk/IntrinsicSurfaceMapper.h"

#include "vtkPointData.h"
#include "vtkCellData.h"

#if MIRTK_USE_FLOAT_BY_DEFAULT
  #include "vtkFloatArray.h"
#else
  #include "vtkDoubleArray.h"
#endif


namespace mirtk {


// Global flags (cf. mirtk/Options.h)
MIRTK_Common_EXPORT extern int verbose;


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void NearOptimalIntrinsicSurfaceMapper
::CopyAttributes(const NearOptimalIntrinsicSurfaceMapper &other)
{
}

// -----------------------------------------------------------------------------
NearOptimalIntrinsicSurfaceMapper::NearOptimalIntrinsicSurfaceMapper()
{
}

// -----------------------------------------------------------------------------
NearOptimalIntrinsicSurfaceMapper
::NearOptimalIntrinsicSurfaceMapper(const NearOptimalIntrinsicSurfaceMapper &other)
:
  IntrinsicSurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
NearOptimalIntrinsicSurfaceMapper &NearOptimalIntrinsicSurfaceMapper
::operator =(const NearOptimalIntrinsicSurfaceMapper &other)
{
  if (this != &other) {
    IntrinsicSurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
NearOptimalIntrinsicSurfaceMapper::~NearOptimalIntrinsicSurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void NearOptimalIntrinsicSurfaceMapper::ComputeMap()
{
  const int       m = NumberOfComponents();
  const vtkIdType n = NumberOfPoints();

  vtkSmartPointer<vtkDataArray> u0; // Authalic  map, i.e., lambda=0
  vtkSmartPointer<vtkDataArray> u1; // Conformal map, i.e., lambda=1

  // Compute discrete authalic map
  if (verbose > 0) {
    cout << "\n  Computing discrete authalic map...";
  }
  _Lambda = 0., IntrinsicSurfaceMapper::ComputeMap();
  u0.TakeReference(_Values->NewInstance());
  u0->DeepCopy(_Values);
  if (verbose > 0) {
    cout << "  Computing discrete authalic map... done\n";
  }

  // Compute discrete conformal map
  if (verbose > 0) {
    cout << "\n  Computing discrete conformal map...";
  }
  _Lambda = 1., IntrinsicSurfaceMapper::ComputeMap();
  u1 = _Values;
  if (verbose > 0) {
    cout << "  Computing discrete conformal map... done\n\n";
  }

  // Determine optimal lambda
  _Lambda = this->ComputeLambda(u0, u1);
  const double mu = 1. - _Lambda;
  if (verbose > 0) {
    cout << "  Optimal lambda value         = " << _Lambda << "\n";
    cout.flush();
  }

  // Set near-optimal surface map values
  for (vtkIdType i = 0; i < n; ++i)
  for (int       j = 0; j < m; ++j) {
    _Values->SetComponent(i, j, mu * u0->GetComponent(i, j) + _Lambda * u1->GetComponent(i, j));
  }
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
double NearOptimalIntrinsicSurfaceMapper::Scale(vtkPolyData *surface) const
{
  double bounds[6];
  surface->GetBounds(bounds);
  return 3. / ((bounds[1] - bounds[0]) +
               (bounds[3] - bounds[2]) +
               (bounds[5] - bounds[4]));
}

// -----------------------------------------------------------------------------
double NearOptimalIntrinsicSurfaceMapper::Scale(vtkDataArray *u0, vtkDataArray *u1) const
{
  double range[2], u_min, u_max, v_min, v_max;

  u0->GetRange(range, 0);
  u_min = range[0];
  u_max = range[1];

  u0->GetRange(range, 1);
  v_min = range[0];
  v_max = range[1];

  if (u1 != nullptr) {
    u1->GetRange(range, 0);
    u_min = min(u_min, range[0]);
    u_max = max(u_max, range[1]);

    u1->GetRange(range, 1);
    v_min = min(v_min, range[0]);
    v_max = max(v_max, range[1]);
  }

  return 2. / ((u_max - u_min) + (v_max - v_min));
}


} // namespace mirtk
