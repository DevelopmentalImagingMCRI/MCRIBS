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

#include "mirtk/HarmonicSurfaceMapper.h"

#include "mirtk/VtkMath.h"
#include "mirtk/Triangle.h"

#include "vtkIdList.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void HarmonicSurfaceMapper::CopyAttributes(const HarmonicSurfaceMapper &other)
{
}

// -----------------------------------------------------------------------------
HarmonicSurfaceMapper::HarmonicSurfaceMapper()
{
}

// -----------------------------------------------------------------------------
HarmonicSurfaceMapper::HarmonicSurfaceMapper(const HarmonicSurfaceMapper &other)
:
  SymmetricWeightsSurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
HarmonicSurfaceMapper &HarmonicSurfaceMapper
::operator =(const HarmonicSurfaceMapper &other)
{
  if (this != &other) {
    SymmetricWeightsSurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
HarmonicSurfaceMapper::~HarmonicSurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
double HarmonicSurfaceMapper::Weight(int i, int j) const
{
  double a[3], b[3], c[3], w;

  int k, l;
  if (GetEdgeNeighborPoints(i, j, k, l) > 2 || k < 0) {
    cerr << this->NameOfType() << "::Weight: Surface mesh must be triangulated!" << endl;
    exit(1);
  }

  // The value computed by the Cotangent function found in Meyer et al. (2002)
  // is identical to the spring constants in Eck et al. (1995) up to a constant
  // factor of 1/4 missing from Eck et al.'s formula.
  _Surface->GetPoint(static_cast<vtkIdType>(i), a);
  _Surface->GetPoint(static_cast<vtkIdType>(k), b);
  _Surface->GetPoint(static_cast<vtkIdType>(j), c);
  w = Triangle::Cotangent(a, b, c);
  if (l >= 0) {
    _Surface->GetPoint(static_cast<vtkIdType>(l), b);
    w += Triangle::Cotangent(a, b, c);
  }

  // Factor 1/2 canceled by multiplying both sides of the linear equations with 2.
  return w;
}


} // namespace mirtk
