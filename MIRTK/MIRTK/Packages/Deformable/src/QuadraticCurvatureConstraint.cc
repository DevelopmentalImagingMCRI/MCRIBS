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

#include "mirtk/QuadraticCurvatureConstraint.h"

#include "mirtk/Array.h"
#include "mirtk/Memory.h"
#include "mirtk/Vector.h"
#include "mirtk/Matrix.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/VtkMath.h"

#include "vtkPoints.h"
#include "vtkDataArray.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(QuadraticCurvatureConstraint);


// =============================================================================
// Auxiliaries
// =============================================================================

namespace QuadraticCurvatureConstraintUtils {


// -----------------------------------------------------------------------------
/// Evaluate error of quadratic fit in normal direction
struct ComputeErrorOfQuadraticFit
{
  typedef RegisteredPointSet::NodeNeighbors NodeNeighbors;

  vtkPoints           *_Points;
  vtkDataArray        *_Normals;
  const NodeNeighbors *_Neighbors;
  vtkDataArray        *_ExternalMagnitude;
  vtkDataArray        *_Residuals;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    const double delta_sigma2 = .25 * .25;

    int       numNbrPts;
    const int *nbrPtIds;
    double     c[3], p[3], n[3], e[3], b, m, delta, w, wsum;

    Vector h; // Distances of neighbors to tangent plane
    Matrix r; // Squared radial distances of neighbors to central node

    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      b = 0.;
      _Neighbors->GetConnectedPoints(ptId, numNbrPts, nbrPtIds);
      if (numNbrPts > 0) {
        _Points   ->GetPoint(ptId, c);
        _Normals  ->GetTuple(ptId, n);
        h.Initialize(numNbrPts);
        r.Initialize(numNbrPts, 2);
        for (int i = 0; i < numNbrPts; ++i) {
          _Points->GetPoint(nbrPtIds[i], p);
          vtkMath::Subtract(p, c, e);
          h(i)    = vtkMath::Dot(e, n);
          r(i, 0) = vtkMath::Dot(e, e) - h(i) * h(i);
          r(i, 1) = 1.;
        }
        r.PseudoInvert();
        if (_ExternalMagnitude) {
          wsum = 0.;
          m = abs(_ExternalMagnitude->GetComponent(ptId, 0));
          for (int i = 0; i < numNbrPts; ++i) {
            delta = (abs(_ExternalMagnitude->GetComponent(nbrPtIds[i], 0)) - m) / (m + 1e-6);
            wsum += (w = exp(-.5 * delta * delta / delta_sigma2));
            b += w * r(1, i) * h(i);
          }
          if (wsum > 0.) b /= wsum;
        } else {
          for (int i = 0; i < numNbrPts; ++i) {
            b += r(1, i) * h(i);
          }
        }
      }
      _Residuals->SetComponent(ptId, 0, b);
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate energy of quadratic curvature term
struct Evaluate
{
  vtkDataArray *_Residuals;
  double        _Sum;

  Evaluate() : _Sum(.0) {}

  Evaluate(const Evaluate &other, split)
  :
    _Residuals(other._Residuals),
    _Sum(.0)
  {}

  void join(const Evaluate &other)
  {
    _Sum += other._Sum;
  }

  void operator ()(const blocked_range<int> &ptIds)
  {
    double b;
    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      b = _Residuals->GetComponent(ptId, 0);
      _Sum += b * b;
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate gradient of quadratic curvature term (i.e., negative force)
struct EvaluateGradient
{
  typedef QuadraticCurvatureConstraint::GradientType GradientType;

  vtkDataArray *_Status;
  vtkDataArray *_Normals;
  vtkDataArray *_Residuals;
  GradientType *_Gradient;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    double b, n[3];
    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      if (_Status && _Status->GetComponent(ptId, 0) == .0) continue;
      _Normals->GetTuple(ptId, n);
      b = _Residuals->GetComponent(ptId, 0);
      _Gradient[ptId] = -b * GradientType(n);
    }
  }
};


} // namespace QuadraticCurvatureConstraintUtils
using namespace QuadraticCurvatureConstraintUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
QuadraticCurvatureConstraint
::QuadraticCurvatureConstraint(const char *name, double weight)
:
  SurfaceConstraint(name, weight)
{
  // QuadraticCurvatureConstraint specific prefixes
  _ParameterPrefix.push_back("Quadratic surface curvature ");
  _ParameterPrefix.push_back("Quadratic surface mesh curvature ");
  _ParameterPrefix.push_back("Quadratic mesh curvature ");
  // Alternative CurvatureConstraint prefixes
  _ParameterPrefix.push_back("Surface curvature ");
  _ParameterPrefix.push_back("Surface bending ");
  _ParameterPrefix.push_back("Mesh curvature ");
  _ParameterPrefix.push_back("Mesh bending ");
  _ParameterPrefix.push_back("Surface mesh curvature ");
  _ParameterPrefix.push_back("Surface mesh bending ");
}

// -----------------------------------------------------------------------------
void QuadraticCurvatureConstraint
::CopyAttributes(const QuadraticCurvatureConstraint &other)
{
}

// -----------------------------------------------------------------------------
QuadraticCurvatureConstraint
::QuadraticCurvatureConstraint(const QuadraticCurvatureConstraint &other)
:
  SurfaceConstraint(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
QuadraticCurvatureConstraint &QuadraticCurvatureConstraint
::operator =(const QuadraticCurvatureConstraint &other)
{
  if (this != &other) {
    SurfaceConstraint::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
QuadraticCurvatureConstraint::~QuadraticCurvatureConstraint()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void QuadraticCurvatureConstraint::Initialize()
{
  // Initialize base class
  SurfaceConstraint::Initialize();

  // Add point data array for residuals of quadratic fit
  AddPointData("Residuals", 1, VTK_FLOAT);
}

// -----------------------------------------------------------------------------
void QuadraticCurvatureConstraint::Update(bool gradient)
{
  if (_NumberOfPoints == 0) return;

  // Update base class
  SurfaceConstraint::Update(gradient);

  // Compute normal coefficients of quadratic fit
  MIRTK_START_TIMING();
  ComputeErrorOfQuadraticFit fit;
  fit._Points            = Points();
  fit._Normals           = Normals();
  fit._Neighbors         = Neighbors();
  fit._ExternalMagnitude = ExternalMagnitude();
  fit._Residuals         = PointData("Residuals");
  parallel_for(blocked_range<int>(0, _NumberOfPoints), fit);
  MIRTK_DEBUG_TIMING(3, "quadratic curvature fitting");
}

// -----------------------------------------------------------------------------
double QuadraticCurvatureConstraint::Evaluate()
{
  if (_NumberOfPoints == 0) return .0;
  MIRTK_START_TIMING();
  QuadraticCurvatureConstraintUtils::Evaluate eval;
  eval._Residuals = PointData("Residuals");
  parallel_reduce(blocked_range<int>(0, _NumberOfPoints), eval);
  MIRTK_DEBUG_TIMING(3, "evaluation of quadratic curvature penalty");
  return eval._Sum / _NumberOfPoints;
}

// -----------------------------------------------------------------------------
void QuadraticCurvatureConstraint
::EvaluateGradient(double *gradient, double step, double weight)
{
  if (_NumberOfPoints == 0) return;

  MIRTK_START_TIMING();
  memset(_Gradient, 0, _NumberOfPoints * sizeof(GradientType));

  QuadraticCurvatureConstraintUtils::EvaluateGradient eval;
  eval._Status    = Status();
  eval._Normals   = Normals();
  eval._Residuals = PointData("Residuals");
  eval._Gradient  = _Gradient;
  parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);

  SurfaceConstraint::EvaluateGradient(gradient, step, 2. * weight / _NumberOfPoints);
  MIRTK_DEBUG_TIMING(3, "evaluation of quadratic curvature force");
}


} // namespace mirtk
