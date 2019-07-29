/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Andreas Schuh
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

#include "mirtk/SurfaceForce.h"

#include "mirtk/Math.h"

#include "vtkAbstractCellLocator.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void SurfaceForce::CopyAttributes(const SurfaceForce &other)
{
}

// -----------------------------------------------------------------------------
SurfaceForce::SurfaceForce(const char *name, double weight)
:
  ExternalForce(name, weight)
{
  _SurfaceForce = true;
}

// -----------------------------------------------------------------------------
SurfaceForce::SurfaceForce(const SurfaceForce &other)
:
  ExternalForce(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
SurfaceForce &SurfaceForce::operator =(const SurfaceForce &other)
{
  if (this != &other) {
    ExternalForce::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
SurfaceForce::~SurfaceForce()
{
}

// =============================================================================
// Auxiliary functions
// =============================================================================

// -----------------------------------------------------------------------------
double SurfaceForce::IntersectWithRay(const double p[3], const double e[3], double l) const
{
  if (l == .0) l = _PointSet->Diameter();

  const double tol  = 1e-3; // unused by VTK locators (at least up to VTK 6.3)
  const double minl = copysign(2.0 * tol, l);

  double p1[3], p2[3], x[3], pcoords[3], t;
  int    subId;

  vtkAbstractCellLocator *locator = _PointSet->SurfaceCellLocator();

  p1[0] = p[0] + minl * e[0];
  p1[1] = p[1] + minl * e[1];
  p1[2] = p[2] + minl * e[2];
  p2[0] = p[0] + l * e[0];
  p2[1] = p[1] + l * e[1];
  p2[2] = p[2] + l * e[2];

  if (locator->IntersectWithLine(p1, p2, tol, t, x, pcoords, subId)) {
    return abs(minl) + t * abs(l - minl);
  } else {
    return abs(l);
  }
}

// -----------------------------------------------------------------------------
double SurfaceForce::SelfDistance(const double p[3], const double n[3], double maxd) const
{
  if (maxd == .0) maxd = _PointSet->Diameter();
  return min(IntersectWithRay(p, n, +maxd), IntersectWithRay(p, n, -maxd));
}


} // namespace mirtk
