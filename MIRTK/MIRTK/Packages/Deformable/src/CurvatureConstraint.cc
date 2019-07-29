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

#include "mirtk/CurvatureConstraint.h"

#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/EdgeTable.h"

#include "mirtk/VtkMath.h"

#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkDataArray.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(CurvatureConstraint);


// =============================================================================
// Auxiliaries
// =============================================================================

namespace CurvatureConstraintUtils {


// -----------------------------------------------------------------------------
/// Compute centroids of adjacent nodes
struct ComputeCentroids
{
  vtkPoints       *_Points;
  const EdgeTable *_EdgeTable;
  vtkPoints       *_Centroids;

  void operator ()(const blocked_range<int> &re) const
  {
    double     c[3], p[3];
    const int *adjPtIds;
    int        numAdjPts;

    for (int ptId = re.begin(); ptId != re.end(); ++ptId) {
      _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPtIds);
      if (numAdjPts > 0) {
        c[0] = c[1] = c[2] = .0;
        for (int i = 0; i < numAdjPts; ++i) {
          _Points->GetPoint(adjPtIds[i], p);
          c[0] += p[0], c[1] += p[1], c[2] += p[2];
        }
        c[0] /= numAdjPts, c[1] /= numAdjPts, c[2] /= numAdjPts;
      } else {
        _Points->GetPoint(ptId, c);
      }
      _Centroids->SetPoint(ptId, c);
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate bending penalty
struct Evaluate
{
  vtkPoints       *_Points;
  const EdgeTable *_EdgeTable;
  vtkPoints       *_Centroids;
  double           _Sum;

  Evaluate() : _Sum(.0) {}

  Evaluate(const Evaluate &other, split)
  :
    _Points(other._Points),
    _EdgeTable(other._EdgeTable),
    _Centroids(other._Centroids),
    _Sum(.0)
  {}

  void join(const Evaluate &other)
  {
    _Sum += other._Sum;
  }

  void operator ()(const blocked_range<int> &re)
  {
    double p[3], c[3];
    for (int ptId = re.begin(); ptId != re.end(); ++ptId) {
      _Points->GetPoint(ptId, p);
      _Centroids->GetPoint(ptId, c);
      _Sum += vtkMath::Distance2BetweenPoints(c, p);
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate gradient of bending penalty (i.e., negative internal bending force)
struct EvaluateGradient
{
  typedef CurvatureConstraint::GradientType Force;

  vtkPoints       *_Points;
  vtkPoints       *_Centroids;
  vtkDataArray    *_Status;
  const EdgeTable *_EdgeTable;
  Force           *_Gradient;

  void operator ()(const blocked_range<int> &re) const
  {
    int        numAdjPts;
    const int *adjPtIds;
    double     p[3], c[3], w;

    for (int ptId = re.begin(); ptId != re.end(); ++ptId) {
      if (_Status && _Status->GetComponent(ptId, 0) == .0) continue;
      // Derivative of sum terms of adjacent points
      _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPtIds);
      for (int i = 0; i < numAdjPts; ++i) {
        _Points->GetPoint(adjPtIds[i], p);
        _Centroids->GetPoint(adjPtIds[i], c);
        w = 1.0 / _EdgeTable->NumberOfAdjacentPoints(adjPtIds[i]);
        _Gradient[ptId] += w * Force(c[0] - p[0], c[1] - p[1], c[2] - p[2]);
      }
      // Derivative of sum term of this point
      _Points->GetPoint(ptId, p);
      _Centroids->GetPoint(ptId, c);
      _Gradient[ptId] -= Force(c[0] - p[0], c[1] - p[1], c[2] - p[2]);
    }
  }
};


} // namespace CurvatureConstraintUtils

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
CurvatureConstraint::CurvatureConstraint(const char *name, double weight)
:
  SurfaceConstraint(name, weight)
{
  _ParameterPrefix.push_back("Surface curvature ");
  _ParameterPrefix.push_back("Surface bending ");
  _ParameterPrefix.push_back("Mesh curvature ");
  _ParameterPrefix.push_back("Mesh bending ");
  _ParameterPrefix.push_back("Surface mesh curvature ");
  _ParameterPrefix.push_back("Surface mesh bending ");
}

// -----------------------------------------------------------------------------
void CurvatureConstraint::CopyAttributes(const CurvatureConstraint &other)
{
  if (other._Centroids) {
    _Centroids.TakeReference(other._Centroids->NewInstance());
    _Centroids->DeepCopy(other._Centroids);
  } else {
    _Centroids = nullptr;
  }
}

// -----------------------------------------------------------------------------
CurvatureConstraint::CurvatureConstraint(const CurvatureConstraint &other)
:
  SurfaceConstraint(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
CurvatureConstraint &CurvatureConstraint::operator =(const CurvatureConstraint &other)
{
  if (this != &other) {
    SurfaceConstraint::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
CurvatureConstraint::~CurvatureConstraint()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void CurvatureConstraint::Initialize()
{
  // Initialize base class
  SurfaceConstraint::Initialize();

  // Initialize this class
  CurvatureConstraint::Init();
}

// -----------------------------------------------------------------------------
void CurvatureConstraint::Reinitialize()
{
  // Reinitialize base class
  SurfaceConstraint::Reinitialize();

  // Reinitialize this class
  CurvatureConstraint::Init();
}

// -----------------------------------------------------------------------------
void CurvatureConstraint::Init()
{
  if (_Centroids == nullptr) {
    _Centroids = vtkSmartPointer<vtkPoints>::New();
  }
  _Centroids->SetNumberOfPoints(_NumberOfPoints);
}

// -----------------------------------------------------------------------------
void CurvatureConstraint::Update(bool gradient)
{
  // Update base class
  SurfaceConstraint::Update(gradient);

  // Update centroids
  MIRTK_START_TIMING();
  CurvatureConstraintUtils::ComputeCentroids eval;
  eval._Points    = _PointSet->SurfacePoints();
  eval._EdgeTable = _PointSet->SurfaceEdges();
  eval._Centroids = _Centroids;
  parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);
  MIRTK_DEBUG_TIMING(3, "update of centroids");
}

// -----------------------------------------------------------------------------
double CurvatureConstraint::Evaluate()
{
  if (_NumberOfPoints == 0) return .0;
  MIRTK_START_TIMING();
  CurvatureConstraintUtils::Evaluate eval;
  eval._Points    = _PointSet->SurfacePoints();
  eval._EdgeTable = _PointSet->SurfaceEdges();
  eval._Centroids = _Centroids;
  parallel_reduce(blocked_range<int>(0, _NumberOfPoints), eval);
  MIRTK_DEBUG_TIMING(3, "evaluation of curvature penalty");
  return eval._Sum / _NumberOfPoints;
}

// -----------------------------------------------------------------------------
void CurvatureConstraint::EvaluateGradient(double *gradient, double step, double weight)
{
  if (_NumberOfPoints == 0) return;

  MIRTK_START_TIMING();
  memset(_Gradient, 0, _NumberOfPoints * sizeof(GradientType));

  CurvatureConstraintUtils::EvaluateGradient eval;
  eval._Points    = _PointSet->SurfacePoints();
  eval._Status    = _PointSet->SurfaceStatus();
  eval._EdgeTable = _PointSet->SurfaceEdges();
  eval._Centroids = _Centroids;
  eval._Gradient  = _Gradient;
  parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);

  SurfaceConstraint::EvaluateGradient(gradient, step, 2.0 * weight / _NumberOfPoints);
  MIRTK_DEBUG_TIMING(3, "evaluation of curvature force");
}


} // namespace mirtk
