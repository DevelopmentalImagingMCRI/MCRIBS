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

#include "mirtk/MaximumCurvatureConstraint.h"

#include "mirtk/Math.h"
#include "mirtk/MeshSmoothing.h"
#include "mirtk/SurfaceCurvature.h"

#include "vtkPointData.h"
#include "vtkDataArray.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(MaximumCurvatureConstraint);


// =============================================================================
// Auxiliaries
// =============================================================================

namespace MaximumCurvatureConstraintUtils {


// -----------------------------------------------------------------------------
/// Evaluate constraint penalty
struct EvaluatePenalty
{
  vtkDataArray *_Curvature;
  double        _Penalty;

  EvaluatePenalty()
  :
    _Penalty(0.)
  {}

  EvaluatePenalty(const EvaluatePenalty &other, split)
  :
    _Curvature(other._Curvature),
    _Penalty(0.)
  {}

  void join(const EvaluatePenalty &other)
  {
    _Penalty += other._Penalty;
  }

  void operator ()(const blocked_range<int> &ptIds)
  {
    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      _Penalty += abs(_Curvature->GetComponent(ptId, 0));
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate negative constraint force (i.e., gradient of constraint term)
///
/// Spring force designed to minimize maximum curvature with magnitude
/// proportional to the absolute value of the maximum curvature.
struct EvaluateForce
{
  typedef MaximumCurvatureConstraint::GradientType GradientType;

  vtkPoints       *_Points;
  vtkDataArray    *_Status;
  const EdgeTable *_EdgeTable;
  vtkDataArray    *_Curvature;
  double           _Threshold;
  GradientType    *_Gradient;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    double     m, k;
    int        numAdjPts;
    const int *adjPtIds;
    Point      p, q;
    Vector3    f;

    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      if (_Status && _Status->GetComponent(ptId, 0) == .0) continue;
      _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPtIds);
      if (numAdjPts > 0) {
        // Magnitude of spring force
        k = _Curvature->GetComponent(ptId, 0);
        m = SShapedMembershipFunction(abs(k), 0., _Threshold);
        // Compute curvature weighted spring force
        _Points->GetPoint(ptId, p);
        f = 0.;
        for (int i = 0; i < numAdjPts; ++i) {
          _Points->GetPoint(adjPtIds[i], q);
          f += q - p;
        }
        f.Normalize();
        _Gradient[ptId] = -m * GradientType(f);
      }
    }
  }
};


} // namespace MaximumCurvatureConstraintUtils

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void MaximumCurvatureConstraint::CopyAttributes(const MaximumCurvatureConstraint &other)
{
  _Threshold = other._Threshold;
}

// -----------------------------------------------------------------------------
MaximumCurvatureConstraint::MaximumCurvatureConstraint(const char *name, double weight)
:
  SurfaceConstraint(name, weight),
  _Threshold(1.)
{
  _ParameterPrefix.push_back("Maximum curvature ");
}

// -----------------------------------------------------------------------------
MaximumCurvatureConstraint
::MaximumCurvatureConstraint(const MaximumCurvatureConstraint &other)
:
  SurfaceConstraint(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
MaximumCurvatureConstraint &MaximumCurvatureConstraint
::operator =(const MaximumCurvatureConstraint &other)
{
  if (this != &other) {
    SurfaceConstraint::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
MaximumCurvatureConstraint::~MaximumCurvatureConstraint()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void MaximumCurvatureConstraint::Initialize()
{
  // Initialize base class
  SurfaceConstraint::Initialize();

  // Add global (i.e., shared) point data array of computed surface curvatures
  const bool global = true;
  AddPointData(SurfaceCurvature::MAXIMUM, 1, VTK_FLOAT, global);
}

// -----------------------------------------------------------------------------
void MaximumCurvatureConstraint::Update(bool gradient)
{
  // Update base class
  SurfaceConstraint::Update(gradient);

  if (_NumberOfPoints == 0) return;

  // Update mean curvature
  vtkPolyData  * const surface   = DeformedSurface();
  vtkDataArray * const curvature = PointData(SurfaceCurvature::MAXIMUM);
  if (curvature->GetMTime() < surface->GetMTime()) {

    SurfaceCurvature curv(SurfaceCurvature::Maximum);
    curv.Input(surface);
    curv.VtkCurvaturesOn();
    curv.Run();

    MeshSmoothing smoother;
    smoother.Input(curv.Output());
    smoother.SmoothPointsOff();
    smoother.SmoothArray(SurfaceCurvature::MAXIMUM);
    smoother.NumberOfIterations(2);
    smoother.Run();

    vtkPointData * const smoothPD = smoother.Output()->GetPointData();
    curvature->DeepCopy(smoothPD->GetArray(SurfaceCurvature::MAXIMUM));
    curvature->Modified();
  }
}

// -----------------------------------------------------------------------------
double MaximumCurvatureConstraint::Evaluate()
{
  if (_NumberOfPoints == 0) return 0.;
  MaximumCurvatureConstraintUtils::EvaluatePenalty eval;
  eval._Curvature = PointData(SurfaceCurvature::MAXIMUM);
  parallel_reduce(blocked_range<int>(0, _NumberOfPoints), eval);
  return eval._Penalty / _NumberOfPoints;
}

// -----------------------------------------------------------------------------
void MaximumCurvatureConstraint
::EvaluateGradient(double *gradient, double step, double weight)
{
  if (_NumberOfPoints == 0) return;

  memset(_Gradient, 0, _NumberOfPoints * sizeof(GradientType));

  MaximumCurvatureConstraintUtils::EvaluateForce eval;
  eval._Points    = Points();
  eval._Status    = Status();
  eval._EdgeTable = Edges();
  eval._Curvature = PointData(SurfaceCurvature::MAXIMUM);
  eval._Threshold = _Threshold;
  eval._Gradient  = _Gradient;
  parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);

  SurfaceConstraint::EvaluateGradient(gradient, step, weight / _NumberOfPoints);
}


} // namespace mirtk
