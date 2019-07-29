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

#include "mirtk/GaussCurvatureConstraint.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/MeshSmoothing.h"
#include "mirtk/SurfaceCurvature.h"

#include "mirtk/VtkMath.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(GaussCurvatureConstraint);


// =============================================================================
// Auxiliaries
// =============================================================================

namespace GaussCurvatureConstraintUtils {


// -----------------------------------------------------------------------------
/// Evaluate constraint penalty
struct Evaluate
{
  vtkDataArray *_GaussCurvature;
  double        _Penalty;

  Evaluate()
  :
    _Penalty(0.)
  {}

  Evaluate(const Evaluate &other, split)
  :
    _GaussCurvature(other._GaussCurvature),
    _Penalty(0.)
  {}

  void join(const Evaluate &other)
  {
    _Penalty += other._Penalty;
  }

  void operator ()(const blocked_range<int> &ptIds)
  {
    double K;
    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      K = _GaussCurvature->GetComponent(ptId, 0);
      _Penalty += abs(K);
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate negative constraint force (i.e., gradient of constraint term)
struct EvaluateGradient
{
  typedef GaussCurvatureConstraint::Action       Action;
  typedef GaussCurvatureConstraint::GradientType GradientType;

  vtkPoints       *_Points;
  vtkDataArray    *_Status;
  const EdgeTable *_EdgeTable;
  vtkDataArray    *_Normals;
  vtkDataArray    *_GaussCurvature;
  vtkDataArray    *_MeanCurvature;
  vtkDataArray    *_ExternalMagnitude;
  GradientType    *_Gradient;
  double           _MinGaussCurvature;
  double           _MaxGaussCurvature;
  Action           _NegativeGaussCurvatureAction;
  Action           _PositiveGaussCurvatureAction;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    double     c[3], p[3], d[3], f[3], f1[3], f2[3], n[3], m, K, H, dp, dp1, dp2, sign;
    int        numAdjPts;
    const int *adjPtIds;

    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      if (_Status && _Status->GetComponent(ptId, 0) == .0) continue;
      _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPtIds);
      if (numAdjPts > 0) {
        // Magnitude of spring force based on Gauss curvature
        K = _GaussCurvature->GetComponent(ptId, 0);
        m = SShapedMembershipFunction(abs(K), _MinGaussCurvature, _MaxGaussCurvature);
        // When mean curvature given, smooth proportional to it
        if (_MeanCurvature) {
          H = _MeanCurvature->GetComponent(ptId, 0);
          // The following settings are for white to pial brain surface inflation,
          // where vallys of sulci should not be smoothed out by this force
          if (H < 0.) m *= 1. - SShapedMembershipFunction(-H, 0., .5);
          else        m *= SShapedMembershipFunction(H, 0., 1.);
        }
        // Compute curvature weighted spring force
        _Points->GetPoint(ptId, c);
        f[0] = f[1] = f[2] = 0.;
        if (K < 0.) {
          switch (_NegativeGaussCurvatureAction)
          {
            case GaussCurvatureConstraint::NoAction: {
              continue;
            } break;
            case GaussCurvatureConstraint::DefaultAction: {
              if (_ExternalMagnitude) {
                sign = sgn(_ExternalMagnitude->GetComponent(ptId, 0));
              } else {
                sign = 0.;
              }
            } break;
            case GaussCurvatureConstraint::Deflate: {
              sign = -1.;
            } break;
            case GaussCurvatureConstraint::Inflate: {
              sign = +1.;
            } break;
          }
          _Normals->GetTuple(ptId, n);
          if (sign == 0.) {
            f1[0] = f1[1] = f1[2] = dp1 = 0.;
            f2[0] = f2[1] = f2[2] = dp2 = 0.;
            for (int i = 0; i < numAdjPts; ++i) {
              _Points->GetPoint(adjPtIds[i], p);
              vtkMath::Subtract(p, c, d);
              dp = vtkMath::Dot(d, n);
              if (dp > 0.) {
                dp1 += dp;
                vtkMath::Add(f1, d, f1);
              } else if (dp < 0.) {
                dp2 += -dp;
                vtkMath::Add(f2, d, f2);
              }
            }
            if (dp1 >= dp2) {
              f[0] = f2[0], f[1] = f2[1], f[2] = f2[2];
            } else {
              f[0] = f1[0], f[1] = f1[1], f[2] = f1[2];
            }
          } else {
            for (int i = 0; i < numAdjPts; ++i) {
              _Points->GetPoint(adjPtIds[i], p);
              vtkMath::Subtract(p, c, d);
              if (vtkMath::Dot(d, n) * sign > 0.) {
                vtkMath::Add(f, d, f);
              }
            }
          }
        } else if (K > 0.) {
          if (_PositiveGaussCurvatureAction == GaussCurvatureConstraint::NoAction) continue;
          for (int i = 0; i < numAdjPts; ++i) {
            _Points->GetPoint(adjPtIds[i], p);
            vtkMath::Subtract(p, c, d);
            vtkMath::Add(f, d, f);
          }
        }
        vtkMath::Normalize(f);
        _Gradient[ptId] = -m * GradientType(f);
      }
    }
  }
};


} // namespace GaussCurvatureConstraintUtils

// =============================================================================
// Enum <-> string conversion
// =============================================================================

// -----------------------------------------------------------------------------
template <>
bool FromString(const char *str, enum GaussCurvatureConstraint::Action &value)
{
  const string lstr = ToLower(str);
  if (lstr == "default") {
    value = GaussCurvatureConstraint::DefaultAction;
  } else if (lstr == "none") {
    value = GaussCurvatureConstraint::NoAction;
  } else if (lstr == "inflate") {
    value = GaussCurvatureConstraint::Inflate;
  } else if (lstr == "deflate") {
    value = GaussCurvatureConstraint::Deflate;
  } else {
    return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
template <>
string ToString(const enum GaussCurvatureConstraint::Action &value, int w, char c, bool left)
{
  const char *str;
  switch (value) {
    case GaussCurvatureConstraint::DefaultAction: {
      str = "Default";
    } break;
    case GaussCurvatureConstraint::NoAction: {
      str = "None";
    } break;
    case GaussCurvatureConstraint::Inflate: {
      str = "Inflate";
    } break;
    case GaussCurvatureConstraint::Deflate: {
      str = "Deflate";
    } break;
  }
  return ToString(str, w, c, left);
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void GaussCurvatureConstraint::CopyAttributes(const GaussCurvatureConstraint &other)
{
  _MinGaussCurvature            = other._MinGaussCurvature;
  _MaxGaussCurvature            = other._MinGaussCurvature;
  _NegativeGaussCurvatureAction = other._NegativeGaussCurvatureAction;
  _PositiveGaussCurvatureAction = other._PositiveGaussCurvatureAction;
  _UseMeanCurvature             = other._UseMeanCurvature;
}

// -----------------------------------------------------------------------------
GaussCurvatureConstraint::GaussCurvatureConstraint(const char *name, double weight)
:
  SurfaceConstraint(name, weight),
  _MinGaussCurvature(0.),
  _MaxGaussCurvature(1.),
  _NegativeGaussCurvatureAction(DefaultAction),
  _PositiveGaussCurvatureAction(DefaultAction),
  _UseMeanCurvature(false)
{
  _ParameterPrefix.push_back("Gauss curvature ");
  _ParameterPrefix.push_back("Gaussian curvature ");
}

// -----------------------------------------------------------------------------
GaussCurvatureConstraint
::GaussCurvatureConstraint(const GaussCurvatureConstraint &other)
:
  SurfaceConstraint(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
GaussCurvatureConstraint &GaussCurvatureConstraint
::operator =(const GaussCurvatureConstraint &other)
{
  if (this != &other) {
    SurfaceConstraint::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
GaussCurvatureConstraint::~GaussCurvatureConstraint()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void GaussCurvatureConstraint::Initialize()
{
  // Initialize base class
  SurfaceConstraint::Initialize();

  // Add global (i.e., shared) point data array of computed surface curvatures
  const bool global = true;
  AddPointData(SurfaceCurvature::GAUSS, 1, VTK_FLOAT, global);
  if (_UseMeanCurvature) {
    AddPointData(SurfaceCurvature::MEAN,  1, VTK_FLOAT, global);
  }
}

// -----------------------------------------------------------------------------
void GaussCurvatureConstraint::Update(bool gradient)
{
  // Update base class
  SurfaceConstraint::Update(gradient);

  if (_NumberOfPoints == 0) return;

  // Update Gauss and/or mean curvature
  vtkPolyData  * const surface         = DeformedSurface();
  vtkDataArray * const gauss_curvature = PointData(SurfaceCurvature::GAUSS);
  vtkDataArray * const mean_curvature  = (_UseMeanCurvature ? PointData(SurfaceCurvature::MEAN) : nullptr);

  int curv_types = 0;
  if (gauss_curvature->GetMTime() < surface->GetMTime()) {
    curv_types |= SurfaceCurvature::Gauss;
  }
  if (mean_curvature && mean_curvature->GetMTime() < surface->GetMTime()) {
    curv_types |= SurfaceCurvature::Mean;
  }
  if (curv_types != 0) {
    SurfaceCurvature curv(curv_types);
    curv.Input(surface);
    curv.EdgeTable(SharedEdgeTable());
    curv.VtkCurvaturesOn();
    curv.Run();

    MeshSmoothing smoother;
    smoother.Input(curv.Output());
    smoother.EdgeTable(SharedEdgeTable());
    smoother.SmoothPointsOff();
    if ((curv_types & SurfaceCurvature::Gauss) != 0) {
      smoother.SmoothArray(SurfaceCurvature::GAUSS);
    }
    if ((curv_types & SurfaceCurvature::Mean) != 0) {
      smoother.SmoothArray(SurfaceCurvature::MEAN);
    }
    smoother.NumberOfIterations(2);
    smoother.Run();

    vtkPointData * const outputPD = smoother.Output()->GetPointData();
    if ((curv_types & SurfaceCurvature::Gauss) != 0) {
      gauss_curvature->CopyComponent(0, outputPD->GetArray(SurfaceCurvature::GAUSS), 0);
      gauss_curvature->Modified();
    }
    if ((curv_types & SurfaceCurvature::Mean) != 0) {
      mean_curvature->CopyComponent(0, outputPD->GetArray(SurfaceCurvature::MEAN), 0);
      mean_curvature->Modified();
    }
  }
}

// -----------------------------------------------------------------------------
double GaussCurvatureConstraint::Evaluate()
{
  if (_NumberOfPoints == 0) return 0.;
  GaussCurvatureConstraintUtils::Evaluate eval;
  eval._GaussCurvature = PointData(SurfaceCurvature::GAUSS);
  parallel_reduce(blocked_range<int>(0, _NumberOfPoints), eval);
  return eval._Penalty / _NumberOfPoints;
}

// -----------------------------------------------------------------------------
void GaussCurvatureConstraint
::EvaluateGradient(double *gradient, double step, double weight)
{
  if (_NumberOfPoints == 0) return;

  memset(_Gradient, 0, _NumberOfPoints * sizeof(GradientType));

  GaussCurvatureConstraintUtils::EvaluateGradient eval;
  eval._Points            = Points();
  eval._Status            = Status();
  eval._EdgeTable         = Edges();
  eval._Normals           = Normals();
  eval._GaussCurvature    = PointData(SurfaceCurvature::GAUSS);
  eval._MeanCurvature     = (_UseMeanCurvature ? PointData(SurfaceCurvature::MEAN) : nullptr);
  eval._ExternalMagnitude = ExternalMagnitude();
  eval._Gradient          = _Gradient;
  eval._MinGaussCurvature = _MinGaussCurvature;
  eval._MaxGaussCurvature = _MaxGaussCurvature;
  eval._NegativeGaussCurvatureAction = _NegativeGaussCurvatureAction;
  eval._PositiveGaussCurvatureAction = _PositiveGaussCurvatureAction;
  parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);

  SurfaceConstraint::EvaluateGradient(gradient, step, weight / _NumberOfPoints);
}


} // namespace mirtk
