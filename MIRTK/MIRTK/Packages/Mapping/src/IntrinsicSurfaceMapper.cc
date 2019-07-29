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

#include "mirtk/IntrinsicSurfaceMapper.h"

#include "mirtk/Math.h"
#include "mirtk/Triangle.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void IntrinsicSurfaceMapper::CopyAttributes(const IntrinsicSurfaceMapper &other)
{
  _Lambda = other._Lambda;
}

// -----------------------------------------------------------------------------
IntrinsicSurfaceMapper::IntrinsicSurfaceMapper(double lambda)
:
  _Lambda(clamp(lambda, 0., 1.))
{
}

// -----------------------------------------------------------------------------
IntrinsicSurfaceMapper::IntrinsicSurfaceMapper(const IntrinsicSurfaceMapper &other)
:
  NonSymmetricWeightsSurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
IntrinsicSurfaceMapper &IntrinsicSurfaceMapper::operator =(const IntrinsicSurfaceMapper &other)
{
  if (this != &other) {
    NonSymmetricWeightsSurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
IntrinsicSurfaceMapper::~IntrinsicSurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
double IntrinsicSurfaceMapper::Weight(int i, int j) const
{
  double a[3], b[3], c[3], w_conformal = 0., w_authalic = 0.;

  int k, l;
  if (GetEdgeNeighborPoints(i, j, k, l) > 2 || k < 0) {
    cerr << this->NameOfType() << "::Weight: Surface mesh must be triangulated!" << endl;
    exit(1);
  }

  _Surface->GetPoint(static_cast<vtkIdType>(i), a);
  _Surface->GetPoint(static_cast<vtkIdType>(j), b);

  double mu = (1. - _Lambda);
  if (mu != 0.) mu /= vtkMath::Distance2BetweenPoints(a, b);

  _Surface->GetPoint(static_cast<vtkIdType>(k), c);
  if (_Lambda != 0.) w_conformal = Triangle::Cotangent(a, c, b);
  if (mu      != 0.) w_authalic  = Triangle::Cotangent(a, b, c);

  if (l >= 0) {
    _Surface->GetPoint(static_cast<vtkIdType>(l), c);
    if (_Lambda != 0.) w_conformal += Triangle::Cotangent(a, c, b);
    if (mu      != 0.) w_authalic  += Triangle::Cotangent(a, b, c);
  }

  return _Lambda * w_conformal + mu * w_authalic;
}


} // namespace mirtk
