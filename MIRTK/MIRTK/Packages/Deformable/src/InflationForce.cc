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

#include "mirtk/InflationForce.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/ObjectFactory.h"
#include "mirtk/VtkMath.h"

#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(InflationForce);


// =============================================================================
// Auxiliaries
// =============================================================================

namespace InflationForceUtils {


// -----------------------------------------------------------------------------
/// Evaluate penalty
struct Evaluate
{
  vtkPoints       *_Points;
  const EdgeTable *_EdgeTable;
  double           _Sum;

  Evaluate() : _Sum(.0) {}

  Evaluate(const Evaluate &other, split)
  :
    _Points(other._Points),
    _EdgeTable(other._EdgeTable),
    _Sum(.0)
  {}

  void join(const Evaluate &other)
  {
    _Sum += other._Sum;
  }

  // TODO: Implement energy term evaluation whose gradient actually corresponds
  //       to the force computed by EvaluateGradient
  void operator ()(const blocked_range<int> &ptIds)
  {
    double     p1[3], p2[3], s;
    int        numAdjPts;
    const int *adjPtIds;

    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPtIds);
      if (numAdjPts > 0) {
        s = .0;
        _Points->GetPoint(ptId, p1);
        for (int i = 0; i < numAdjPts; ++i) {
          _Points->GetPoint(adjPtIds[i], p2);
          s += vtkMath::Distance2BetweenPoints(p1, p2);
        }
        _Sum += s / numAdjPts;
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate gradient of penalty term (i.e., negative inflation force)
struct EvaluateGradient
{
  typedef InflationForce::GradientType GradientType;

  vtkPoints       *_Points;
  vtkDataArray    *_Status;
  const EdgeTable *_EdgeTable;
  GradientType    *_Gradient;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    const int   *adjPtIds;
    int          numAdjPts;
    double       p1[3], p2[3], e[3];
    GradientType gradient;

    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      if (_Status && _Status->GetComponent(ptId, 0) == .0) continue;
      _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPtIds);
      if (numAdjPts > 0) {
        gradient = .0;
        _Points->GetPoint(ptId, p1);
        for (int i = 0; i < numAdjPts; ++i) {
          _Points->GetPoint(adjPtIds[i], p2);
          vtkMath::Subtract(p1, p2, e);
          gradient._x += e[0];
          gradient._y += e[1];
          gradient._z += e[2];
        }
        gradient /= numAdjPts;
        _Gradient[ptId] = gradient;
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Determine magnitude of average normal component of inflation force
struct SumMagnitudeOfNormalComponents
{
  typedef InflationForce::GradientType GradientType;

  const GradientType *_Gradient;
  vtkDataArray       *_Normals;
  double              _Sum;

  SumMagnitudeOfNormalComponents() : _Sum(.0) {}

  SumMagnitudeOfNormalComponents(const SumMagnitudeOfNormalComponents &other, split)
  :
    _Gradient(other._Gradient),
    _Normals(other._Normals),
    _Sum(.0)
  {}

  void join(const SumMagnitudeOfNormalComponents &other)
  {
    _Sum += other._Sum;
  }

  void operator ()(const blocked_range<int> &ptIds)
  {
    double n[3];
    const GradientType *g = _Gradient + ptIds.begin();
    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId, ++g) {
      _Normals->GetTuple(ptId, n);
      _Sum += g->_x * n[0] + g->_y * n[1] + g->_z * n[2];
    }
  }
};

// -----------------------------------------------------------------------------
/// Subtract average normal component of inflation force
struct SubtractAverageNormalComponent
{
  typedef InflationForce::GradientType GradientType;

  vtkDataArray *_Normals;
  double        _Magnitude;
  GradientType *_Gradient;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    double n[3];
    GradientType *g = _Gradient + ptIds.begin();
    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId, ++g) {
      _Normals->GetTuple(ptId, n);
      g->_x -= _Magnitude * n[0];
      g->_y -= _Magnitude * n[1];
      g->_z -= _Magnitude * n[2];
    }
  }
};


} // namespace InflationForceUtils

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
InflationForce::InflationForce(const char *name, double weight)
:
  SurfaceConstraint(name, weight)
{
  _ParameterPrefix.push_back("Surface inflation ");
  _ParameterPrefix.push_back("Mesh inflation ");
  _ParameterPrefix.push_back("Surface mesh inflation ");
}

// -----------------------------------------------------------------------------
InflationForce::InflationForce(const InflationForce &other)
:
  SurfaceConstraint(other)
{
}

// -----------------------------------------------------------------------------
InflationForce &InflationForce::operator =(const InflationForce &other)
{
  SurfaceConstraint::operator =(other);
  return *this;
}

// -----------------------------------------------------------------------------
InflationForce::~InflationForce()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double InflationForce::Evaluate()
{
  if (_NumberOfPoints == 0) return .0;
  MIRTK_START_TIMING();
  InflationForceUtils::Evaluate eval;
  eval._Points    = _PointSet->SurfacePoints();
  eval._EdgeTable = _PointSet->SurfaceEdges();
  parallel_reduce(blocked_range<int>(0, _NumberOfPoints), eval);
  MIRTK_DEBUG_TIMING(3, "evaluation of inflation energy");
  return eval._Sum / _NumberOfPoints;
}

// -----------------------------------------------------------------------------
void InflationForce::EvaluateGradient(double *gradient, double step, double weight)
{
  if (_NumberOfPoints == 0) return;

  MIRTK_START_TIMING();
  memset(_Gradient, 0, _NumberOfPoints * sizeof(GradientType));

  InflationForceUtils::EvaluateGradient eval;
  eval._Points    = _PointSet->SurfacePoints();
  eval._Status    = _PointSet->SurfaceStatus();
  eval._EdgeTable = _PointSet->SurfaceEdges();
  eval._Gradient  = _Gradient;
  parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);

  InflationForceUtils::SumMagnitudeOfNormalComponents mag;
  mag._Gradient = _Gradient;
  mag._Normals  = _PointSet->SurfaceNormals();
  parallel_reduce(blocked_range<int>(0, _NumberOfPoints), mag);

  InflationForceUtils::SubtractAverageNormalComponent sub;
  sub._Normals   = _PointSet->SurfaceNormals();
  sub._Magnitude = mag._Sum / _NumberOfPoints;
  sub._Gradient  = _Gradient;
  parallel_for(blocked_range<int>(0, _NumberOfPoints), sub);

  const double dist_scale = sqrt(_PointSet->InputSurfaceArea() / _PointSet->SurfaceArea());
  weight = 2.0 * dist_scale * weight / _NumberOfPoints;

  SurfaceConstraint::EvaluateGradient(gradient, step, weight);
  MIRTK_DEBUG_TIMING(3, "evaluation of inflation force");
}


} // namespace mirtk
