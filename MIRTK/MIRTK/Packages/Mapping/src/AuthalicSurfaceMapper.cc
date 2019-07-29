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

#include "mirtk/AuthalicSurfaceMapper.h"

#include "mirtk/Triangle.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void AuthalicSurfaceMapper::CopyAttributes(const AuthalicSurfaceMapper &)
{
}

// -----------------------------------------------------------------------------
AuthalicSurfaceMapper::AuthalicSurfaceMapper()
{
}

// -----------------------------------------------------------------------------
AuthalicSurfaceMapper::AuthalicSurfaceMapper(const AuthalicSurfaceMapper &other)
:
  NonSymmetricWeightsSurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
AuthalicSurfaceMapper &AuthalicSurfaceMapper::operator =(const AuthalicSurfaceMapper &other)
{
  if (this != &other) {
    NonSymmetricWeightsSurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
AuthalicSurfaceMapper::~AuthalicSurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
double AuthalicSurfaceMapper::Weight(int i, int j) const
{
  double a[3], b[3], c[3], w;

  int k, l;
  if (GetEdgeNeighborPoints(i, j, k, l) > 2 || k < 0) {
    cerr << this->NameOfType() << "::Weight: Surface mesh must be triangulated!" << endl;
    exit(1);
  }

  _Surface->GetPoint(static_cast<vtkIdType>(i), a);
  _Surface->GetPoint(static_cast<vtkIdType>(j), b);
  _Surface->GetPoint(static_cast<vtkIdType>(k), c);
  w = Triangle::Cotangent(a, b, c);
  if (l >= 0) {
    _Surface->GetPoint(static_cast<vtkIdType>(l), c);
    w += Triangle::Cotangent(a, b, c);
  }

  return w / vtkMath::Distance2BetweenPoints(a, b);
}


} // namespace mirtk
