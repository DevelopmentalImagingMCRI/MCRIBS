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

#include "mirtk/ShapePreservingSurfaceMapper.h"

#include "mirtk/Math.h"
#include "mirtk/Matrix.h"
#include "mirtk/Vector3D.h"
#include "mirtk/Triangle.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void ShapePreservingSurfaceMapper
::CopyAttributes(const ShapePreservingSurfaceMapper &other)
{
}

// -----------------------------------------------------------------------------
ShapePreservingSurfaceMapper::ShapePreservingSurfaceMapper()
{
}

// -----------------------------------------------------------------------------
ShapePreservingSurfaceMapper::ShapePreservingSurfaceMapper(const ShapePreservingSurfaceMapper &other)
:
  NonSymmetricWeightsSurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
ShapePreservingSurfaceMapper &ShapePreservingSurfaceMapper
::operator =(const ShapePreservingSurfaceMapper &other)
{
  if (this != &other) {
    NonSymmetricWeightsSurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
ShapePreservingSurfaceMapper::~ShapePreservingSurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void ShapePreservingSurfaceMapper
::Weights(int i, const int *j, double *w_i, int d_i) const
{
  if (d_i < 3) {
    cerr << this->NameOfType() << "::Weights: Every point must have at least three neighbors!" << endl;
    exit(1);
  }

  // Compute side lengths and arcs
  Vector l(d_i);
  Vector a(d_i);
  {
    const class Point q = Point(i);
    class Point       p0, p1, p2;
    Vector3D<double>  e0, e1, e2;
    p0 = Point(j[0]);
    e0 = Vector3D<double>(p0 - q);
    e0.Normalize();
    p1 = p0, e1 = e0;
    for (int k = 1; k < d_i; ++k) {
      p2 = Point(j[k]);
      e2 = Vector3D<double>(p2 - q);
      e2.Normalize();
      l(k-1) = p1.Distance(p2);
      a(k-1) = acos(clamp(e1.DotProduct(e2), -1., 1.));
      p1 = p2, e1 = e2;
    }
    l(d_i-1) = p1.Distance(p0);
    a(d_i-1) = acos(clamp(e1.DotProduct(e0), -1., 1.));
  }

  // Step i) Project subgraph into plane
  Matrix p(2, d_i);
  {
    double alpha = 0.;
    a *= two_pi / a.Sum();
    for (int k = 0; k < d_i; ++k) {
      p.Col(k)[0] = l(k) * cos(alpha);
      p.Col(k)[1] = l(k) * sin(alpha);
      alpha += a(k);
    }
  }

  // Step ii) Average barycentric coordinates
  const double center[2] = {0., 0.};
  double alpha, area, area1, area2, area3;
  for (int k = 0; k < d_i; ++k) {
    w_i[k] = 0.;
  }
  for (int c1 = 0, c2, c3; c1 < d_i; ++c1) {
    alpha = a(c1);
    c3 = (c1 + 1) % d_i;
    while (alpha < pi) {
      alpha += a(c3);
      c3 = (c3 + 1) % d_i;
    }
    c2 = (c3 == 0 ? d_i : c3) - 1;
    area  = Triangle::Area2D(p.Col(c1), p.Col(c2), p.Col(c3));
    area1 = Triangle::Area2D(center,    p.Col(c2), p.Col(c3));
    area2 = Triangle::Area2D(p.Col(c1), center,    p.Col(c3));
    area3 = Triangle::Area2D(p.Col(c1), p.Col(c2), center);
    w_i[c1] += area1 / area;
    w_i[c2] += area2 / area;
    w_i[c3] += area3 / area;
  }
  for (int k = 0; k < d_i; ++k) {
    w_i[k] /= d_i;
  }
}


} // namespace mirtk
