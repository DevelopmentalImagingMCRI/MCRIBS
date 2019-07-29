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

#include "mirtk/MeanValueSurfaceMapper.h"

#include "mirtk/Math.h"
#include "mirtk/VtkMath.h"
#include "mirtk/PointSetUtils.h"

#include "vtkSmartPointer.h"
#include "vtkIdList.h"
#include "vtkPolyData.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void MeanValueSurfaceMapper::CopyAttributes(const MeanValueSurfaceMapper &other)
{
}

// -----------------------------------------------------------------------------
MeanValueSurfaceMapper::MeanValueSurfaceMapper()
{
}

// -----------------------------------------------------------------------------
MeanValueSurfaceMapper::MeanValueSurfaceMapper(const MeanValueSurfaceMapper &other)
:
  NonSymmetricWeightsSurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
MeanValueSurfaceMapper &MeanValueSurfaceMapper::operator =(const MeanValueSurfaceMapper &other)
{
  if (this != &other) {
    NonSymmetricWeightsSurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
MeanValueSurfaceMapper::~MeanValueSurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
// The following uses the equality tan(alpha/2) = (1 - cos(alpha)) / sin(alpha)
// to compute the required tangens weights using only inner and outer vector products
double MeanValueSurfaceMapper::Weight(int i, int j) const
{
  double a[3], b[3], c[3], ab[3], ac[3], n[3], d_ab, d_ac, n_norm, w;

  int k, l;
  if (GetEdgeNeighborPoints(i, j, k, l) > 2 || k < 0) {
    cerr << this->NameOfType() << "::Weight: Surface mesh must be triangulated!" << endl;
    exit(1);
  }

  _Surface->GetPoint(static_cast<vtkIdType>(i), a);
  _Surface->GetPoint(static_cast<vtkIdType>(j), b);
  vtkMath::Subtract(b, a, ab);
  d_ab = vtkMath::Norm(ab);

  _Surface->GetPoint(static_cast<vtkIdType>(k), c);
  vtkMath::Subtract(c, a, ac);
  vtkMath::Cross(ab, ac, n);
  d_ac   = vtkMath::Norm(ac);
  n_norm = vtkMath::Norm(n);

  w = (d_ab * d_ac - vtkMath::Dot(ab, ac)) / n_norm;

  if (l >= 0) {
    _Surface->GetPoint(static_cast<vtkIdType>(l), c);
    vtkMath::Subtract(c, a, ac);
    vtkMath::Cross(ab, ac, n);
    d_ac   = vtkMath::Norm(ac);
    n_norm = vtkMath::Norm(n);
    w += (d_ab * d_ac - vtkMath::Dot(ab, ac)) / n_norm;
  }

  return w / d_ab;
}


} // namespace mirtk
