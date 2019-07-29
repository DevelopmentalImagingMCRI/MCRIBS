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

#include "mirtk/SpringForce.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/VtkMath.h"

#include "vtkPoints.h"
#include "vtkDataArray.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(SpringForce);


// =============================================================================
// Auxiliaries
// =============================================================================

namespace SpringForceUtils {


// -----------------------------------------------------------------------------
/// Evaluate spring force term
struct Evaluate
{
  vtkPoints       *_Points;
  const EdgeTable *_EdgeTable;
  double           _SSE;

  Evaluate() : _SSE(.0) {}

  Evaluate(const Evaluate &other, split)
  :
    _Points(other._Points),
    _EdgeTable(other._EdgeTable),
    _SSE(.0)
  {}

  void join(const Evaluate &other)
  {
    _SSE += other._SSE;
  }

  void operator ()(const blocked_range<int> &ptIds)
  {
    int        numAdjPts;
    const int *adjPtIds;
    double     c[3], p[3], sum;

    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      _Points->GetPoint(ptId, c);
      _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPtIds);
      if (numAdjPts > 0) {
        sum = .0;
        for (int i = 0; i < numAdjPts; ++i) {
          _Points->GetPoint(adjPtIds[i], p);
          sum += vtkMath::Distance2BetweenPoints(c, p);
        }
        _SSE += sum / numAdjPts;
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate spring force term
struct EvaluateWithWeightedComponents
{
  vtkPoints       *_Points;
  vtkDataArray    *_Normals;
  const EdgeTable *_EdgeTable;
  double           _InwardNormalWeight;
  double           _OutwardNormalWeight;
  double           _TangentialWeight;
  double           _SSE;

  EvaluateWithWeightedComponents() : _SSE(.0) {}

  EvaluateWithWeightedComponents(const EvaluateWithWeightedComponents &other, split)
  :
    _Points(other._Points),
    _Normals(other._Normals),
    _EdgeTable(other._EdgeTable),
    _InwardNormalWeight(other._InwardNormalWeight),
    _OutwardNormalWeight(other._OutwardNormalWeight),
    _TangentialWeight(other._TangentialWeight),
    _SSE(.0)
  {}

  void join(const EvaluateWithWeightedComponents &other)
  {
    _SSE += other._SSE;
  }

  void operator ()(const blocked_range<int> &ptIds)
  {
    int        numAdjPts;
    const int *adjPtIds;
    double     c[3], p[3], n[3], nc, sum, penalty, wsum;

    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      _Points ->GetPoint(ptId, c);
      _Normals->GetTuple(ptId, n);
      _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPtIds);
      if (numAdjPts > 0) {
        sum = .0;
        for (int i = 0; i < numAdjPts; ++i) {
          _Points->GetPoint(adjPtIds[i], p);
          vtkMath::Subtract(c, p, p);
          nc = vtkMath::Dot(p, n);
          p[0] -= nc * n[0];
          p[1] -= nc * n[1];
          p[2] -= nc * n[2];
          penalty = _TangentialWeight * vtkMath::Dot(p, p);
          if (nc < 0.) {
            wsum = _TangentialWeight + _OutwardNormalWeight;
            if (wsum > 0.) {
              penalty += _OutwardNormalWeight * nc * nc;
              penalty /= wsum;
            } else {
              penalty = 0.;
            }
          } else {
            wsum = _TangentialWeight + _InwardNormalWeight;
            if (wsum > 0.) {
              penalty += _InwardNormalWeight * nc * nc;
              penalty /= wsum;
            } else {
              penalty = 0.;
            }
          }
          sum += penalty;
        }
        _SSE += sum / numAdjPts;
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate gradient of spring force term
struct EvaluateGradient
{
  typedef SpringForce::GradientType GradientType;

  vtkPoints       *_Points;
  vtkDataArray    *_Status;
  const EdgeTable *_EdgeTable;
  GradientType    *_Gradient;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    int        numAdjPts;
    const int *adjPtIds;
    double     p[3], c[3];

    GradientType *g = _Gradient + ptIds.begin();
    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId, ++g) {
      if (_Status && _Status->GetComponent(ptId, 0) == .0) continue;
      _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPtIds);
      if (numAdjPts > 0) {
        _Points->GetPoint(ptId, c);
        for (int i = 0; i < numAdjPts; ++i) {
          _Points->GetPoint(adjPtIds[i], p);
          g->_x += c[0] - p[0];
          g->_y += c[1] - p[1];
          g->_z += c[2] - p[2];
        }
        (*g) /= numAdjPts;
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Weight normal and tangential component
struct WeightComponents
{
  typedef SpringForce::GradientType GradientType;

  vtkDataArray *_Normals;
  GradientType *_Gradient;
  double        _InwardNormalWeight;
  double        _OutwardNormalWeight;
  double        _TangentialWeight;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    double n[3], nc, wsum;
    GradientType *g = _Gradient + ptIds.begin();
    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId, ++g) {
      if (g->_x || g->_y || g->_z) {
        _Normals->GetTuple(ptId, n);
        nc = g->_x * n[0] + g->_y * n[1] + g->_z * n[2];
        vtkMath::MultiplyScalar(n, nc);
        if (nc < 0.) {
          wsum = _TangentialWeight + _OutwardNormalWeight;
          if (wsum > 0.) {
            g->_x = _TangentialWeight * (g->_x - n[0]) + _OutwardNormalWeight * n[0];
            g->_y = _TangentialWeight * (g->_y - n[1]) + _OutwardNormalWeight * n[1];
            g->_z = _TangentialWeight * (g->_z - n[2]) + _OutwardNormalWeight * n[2];
            (*g) /= wsum;
          } else {
            g->_x = g->_y = g->_z = 0.;
          }
        } else {
          wsum = _TangentialWeight + _InwardNormalWeight;
          if (wsum > 0.) {
            g->_x = _TangentialWeight * (g->_x - n[0]) + _InwardNormalWeight * n[0];
            g->_y = _TangentialWeight * (g->_y - n[1]) + _InwardNormalWeight * n[1];
            g->_z = _TangentialWeight * (g->_z - n[2]) + _InwardNormalWeight * n[2];
            (*g) /= wsum;
          } else {
            g->_x = g->_y = g->_z = 0.;
          }
        }
      }
    }
  }
};


} // namespace SpringForceUtils

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
SpringForce::SpringForce(const char *name, double weight)
:
  SurfaceConstraint(name, weight),
  _InwardNormalWeight(.5),
  _OutwardNormalWeight(.5),
  _TangentialWeight(.5)
{
  _ParameterPrefix.push_back("Spring force ");
  _ParameterPrefix.push_back("Surface smoothness ");
  _ParameterPrefix.push_back("Surface bending ");
  _ParameterPrefix.push_back("Mesh smoothness ");
  _ParameterPrefix.push_back("Mesh bending ");
  _ParameterPrefix.push_back("Surface mesh smoothness ");
  _ParameterPrefix.push_back("Surface mesh bending ");
}

// -----------------------------------------------------------------------------
void SpringForce::CopyAttributes(const SpringForce &other)
{
  _InwardNormalWeight  = other._InwardNormalWeight;
  _OutwardNormalWeight = other._OutwardNormalWeight;
  _TangentialWeight    = other._TangentialWeight;
}

// -----------------------------------------------------------------------------
SpringForce::SpringForce(const SpringForce &other)
:
  SurfaceConstraint(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
SpringForce &SpringForce::operator =(const SpringForce &other)
{
  if (this != &other) {
    SurfaceConstraint::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
SpringForce::~SpringForce()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double SpringForce::Evaluate()
{
  if (_NumberOfPoints == 0 || (fequal(_InwardNormalWeight,  0.) &&
                               fequal(_OutwardNormalWeight, 0.) &&
                               fequal(_TangentialWeight,    0.))) return .0;

  double sse;
  if (fequal(_InwardNormalWeight,  _TangentialWeight) &&
      fequal(_OutwardNormalWeight, _TangentialWeight)) {
    SpringForceUtils::Evaluate eval;
    eval._Points    = Points();
    eval._EdgeTable = Edges();
    parallel_reduce(blocked_range<int>(0, _NumberOfPoints), eval);
    sse = eval._SSE;
  } else {
    SpringForceUtils::EvaluateWithWeightedComponents eval;
    eval._Points              = Points();
    eval._Normals             = Normals();
    eval._EdgeTable           = Edges();
    eval._InwardNormalWeight  = _InwardNormalWeight;
    eval._OutwardNormalWeight = _OutwardNormalWeight;
    eval._TangentialWeight    = _TangentialWeight;
    parallel_reduce(blocked_range<int>(0, _NumberOfPoints), eval);
    sse = eval._SSE;
  }

  double area_scale = _PointSet->InputSurfaceArea() / _PointSet->SurfaceArea();
  if (IsNaN(area_scale)) area_scale = 1.0;

  return area_scale * sse / _NumberOfPoints;
}

// -----------------------------------------------------------------------------
void SpringForce::EvaluateGradient(double *gradient, double step, double weight)
{
  if (_NumberOfPoints == 0 || (fequal(_InwardNormalWeight,  0.) &&
                               fequal(_OutwardNormalWeight, 0.) &&
                               fequal(_TangentialWeight,    0.))) return;

  memset(_Gradient, 0, _NumberOfPoints * sizeof(GradientType));

  SpringForceUtils::EvaluateGradient eval;
  eval._Points    = Points();
  eval._Status    = Status();
  eval._EdgeTable = Edges();
  eval._Gradient  = _Gradient;
  parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);

  if (fequal(_InwardNormalWeight,  _TangentialWeight) &&
      fequal(_OutwardNormalWeight, _TangentialWeight)) {
    SpringForceUtils::WeightComponents mul;
    mul._Normals             = Normals();
    mul._InwardNormalWeight  = _InwardNormalWeight;
    mul._OutwardNormalWeight = _OutwardNormalWeight;
    mul._TangentialWeight    = _TangentialWeight;
    mul._Gradient            = _Gradient;
    parallel_for(blocked_range<int>(0, _NumberOfPoints), mul);
  }

  double area_scale = _PointSet->InputSurfaceArea() / _PointSet->SurfaceArea();
  if (IsNaN(area_scale)) area_scale = 1.0;

  SurfaceConstraint::EvaluateGradient(gradient, step, 2.0 * area_scale * weight / _NumberOfPoints);
}


} // namespace mirtk
