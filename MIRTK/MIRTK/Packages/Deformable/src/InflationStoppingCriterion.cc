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

#include "mirtk/InflationStoppingCriterion.h"

#include "mirtk/Parallel.h"
#include "mirtk/DeformableSurfaceModel.h"
#include "mirtk/RegisteredPointSet.h"
#include "mirtk/EdgeTable.h"

#include "vtkPoints.h"
#include "vtkPlane.h"
#include "vtkDataArray.h"


namespace mirtk {


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace InflationStoppingCriterionUtils {


// -----------------------------------------------------------------------------
/// Evaluate RMS of distance of adjacent nodes to the tangent plane
struct Evaluate
{
  typedef RegisteredPointSet::NodeNeighbors NodeNeighbors;

  vtkPoints           *_Points;
  vtkDataArray        *_Normals;
  const NodeNeighbors *_Neighbors;
  double               _SSE;
  int                  _Num;

  Evaluate() : _SSE(.0), _Num(0) {}

  Evaluate(const Evaluate &other, split)
  :
    _Points(other._Points),
    _Normals(other._Normals),
    _Neighbors(other._Neighbors),
    _SSE(.0), _Num(0)
  {}

  void join(const Evaluate &other)
  {
    _SSE += other._SSE;
    _Num += other._Num;
  }

  void operator ()(const blocked_range<int> &ptIds)
  {
    int        numNbrPts;
    const int *nbrPtIds;
    double     p[3], n[3], e[3], dp, dist2;
    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      _Points ->GetPoint(ptId, p);
      _Normals->GetTuple(ptId, n);
      _Neighbors->GetConnectedPoints(ptId, numNbrPts, nbrPtIds);
      for (int i = 0; i < numNbrPts; ++i) {
        _Points->GetPoint(nbrPtIds[i], e);
        e[0] -= p[0], e[1] -= p[1], e[2] -= p[2];
        dist2 = e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
        if (dist2 > .0) {
          dp = n[0]*e[0] + n[1]*e[1] + n[2]*e[2];
          _SSE += dp * dp / dist2;
          ++_Num;
        }
      }
    }
  }
};


} // namespace InflationStoppingCriterionUtils
using namespace InflationStoppingCriterionUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
InflationStoppingCriterion::InflationStoppingCriterion(const ObjectiveFunction *f)
:
  StoppingCriterion(f),
  _Threshold(.015),
  _AverageDistance(numeric_limits<double>::infinity())
{
}

// -----------------------------------------------------------------------------
void InflationStoppingCriterion::CopyAttributes(const InflationStoppingCriterion &other)
{
  _Threshold       = other._Threshold;
  _AverageDistance = other._AverageDistance;
}

// -----------------------------------------------------------------------------
InflationStoppingCriterion::InflationStoppingCriterion(const InflationStoppingCriterion &other)
:
  StoppingCriterion(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
InflationStoppingCriterion &InflationStoppingCriterion::operator =(const InflationStoppingCriterion &other)
{
  if (this != &other) {
    StoppingCriterion::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
StoppingCriterion *InflationStoppingCriterion::New() const
{
  return new InflationStoppingCriterion(*this);
}

// -----------------------------------------------------------------------------
InflationStoppingCriterion::~InflationStoppingCriterion()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
bool InflationStoppingCriterion::Fulfilled(int, double, const double *)
{
  const DeformableSurfaceModel *model;
  model = dynamic_cast<const DeformableSurfaceModel *>(_Function);
  if (!model) return false;

  if (model->NumberOfPoints() == 0) return true;

  InflationStoppingCriterionUtils::Evaluate eval;
  eval._Points    = model->PointSet().SurfacePoints();
  eval._Normals   = model->PointSet().SurfaceNormals();
  eval._Neighbors = model->PointSet().SurfaceNeighbors();
  parallel_reduce(blocked_range<int>(0, model->NumberOfPoints()), eval);
  _AverageDistance = (eval._Num > 0 ? sqrt(eval._SSE / eval._Num) : .0);
  return _AverageDistance < _Threshold;
}

// =============================================================================
// Logging
// =============================================================================

// -----------------------------------------------------------------------------
void InflationStoppingCriterion::Print(ostream &os) const
{
  const ios::fmtflags fmt = os.flags();
  os << "inflation RMS error = " << fixed << _AverageDistance;
  os.flags(fmt);
}


} // namespace mirtk
