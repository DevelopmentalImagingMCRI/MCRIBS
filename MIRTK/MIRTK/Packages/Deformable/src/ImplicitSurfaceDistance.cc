/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
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

#include "mirtk/ImplicitSurfaceDistance.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/DataStatistics.h"

#include "vtkPoints.h"
#include "vtkDataArray.h"

using namespace mirtk::data::statistic;


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(ImplicitSurfaceDistance);


// =============================================================================
// Auxiliary functions
// =============================================================================

namespace ImplicitSurfaceDistanceUtils {


// -----------------------------------------------------------------------------
/// Compute magnitude of implicit surface distance force
struct ComputeMagnitude
{
  vtkDataArray *_Status;
  vtkDataArray *_Distances;
  vtkDataArray *_Magnitude;
  double        _Threshold;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    double d, m;

    const vtkIdType begin = static_cast<vtkIdType>(ptIds.begin());
    const vtkIdType end   = static_cast<vtkIdType>(ptIds.end  ());

    for (vtkIdType ptId = begin; ptId != end; ++ptId) {
      if (_Status && _Status->GetComponent(ptId, 0) == 0.) {
        _Magnitude->SetComponent(ptId, 0, 0.);
      } else {
        d = _Distances->GetComponent(ptId, 0);
        m = SShapedMembershipFunction(abs(d), 0., _Threshold);
        _Magnitude->SetComponent(ptId, 0, -copysign(m, d));
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Compute magnitude of implicit surface distance force
struct EvaluateMagnitude
{
  typedef ImplicitSurfaceForce::ImageFunction MagnitudeFunction;

  MagnitudeFunction *_Function;
  vtkPoints         *_Points;
  vtkDataArray      *_Status;
  vtkDataArray      *_Distances;
  vtkDataArray      *_Magnitude;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    double p[3], d, m;

    const vtkIdType begin = static_cast<vtkIdType>(ptIds.begin());
    const vtkIdType end   = static_cast<vtkIdType>(ptIds.end  ());

    for (vtkIdType ptId = begin; ptId != end; ++ptId) {
      if (_Status && _Status->GetComponent(ptId, 0) == 0.) {
        _Magnitude->SetComponent(ptId, 0, 0.);
      } else {
        _Points->GetPoint(ptId, p);
        _Function->WorldToImage(p[0], p[1], p[2]);
        m = _Function->Evaluate(p[0], p[1], p[2]);
        d = _Distances->GetComponent(ptId, 0);
        _Magnitude->SetComponent(ptId, 0, -copysign(m, d));
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Compute gradient of implicit surface distance force term (i.e., negative force)
struct ComputeGradient
{
  typedef ImplicitSurfaceDistance::GradientType GradientType;

  vtkDataArray *_Normals;
  vtkDataArray *_Magnitude;
  GradientType *_Gradient;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    double m, n[3];

    const vtkIdType begin = static_cast<vtkIdType>(ptIds.begin());
    const vtkIdType end   = static_cast<vtkIdType>(ptIds.end  ());

    for (vtkIdType ptId = begin; ptId != end; ++ptId) {
      _Normals->GetTuple(ptId, n);
      m = _Magnitude->GetComponent(ptId, 0);
      _Gradient[ptId] = -m * GradientType(n);
    }
  }
};


} // namespace ImplicitSurfaceDistanceUtils
using namespace ImplicitSurfaceDistanceUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void ImplicitSurfaceDistance::CopyAttributes(const ImplicitSurfaceDistance &other)
{
  _MagnitudeImage     = other._MagnitudeImage;
  _NormalizeMagnitude = other._NormalizeMagnitude;
  _InvertMagnitude    = other._InvertMagnitude;
  _MinThreshold       = other._MinThreshold;
  _MaxThreshold       = other._MaxThreshold;
}

// -----------------------------------------------------------------------------
ImplicitSurfaceDistance::ImplicitSurfaceDistance(const char *name, double weight)
:
  ImplicitSurfaceForce(name, weight),
  _MagnitudeImage(nullptr),
  _NormalizeMagnitude(false),
  _InvertMagnitude(false),
  _MinThreshold(1.),
  _MaxThreshold(0.)
{
  _MaxDistance = 5.;
}

// -----------------------------------------------------------------------------
ImplicitSurfaceDistance::ImplicitSurfaceDistance(const ImplicitSurfaceDistance &other)
:
  ImplicitSurfaceForce(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
ImplicitSurfaceDistance &ImplicitSurfaceDistance::operator =(const ImplicitSurfaceDistance &other)
{
  if (this != &other) {
    ImplicitSurfaceForce::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
ImplicitSurfaceDistance::~ImplicitSurfaceDistance()
{
}

// -----------------------------------------------------------------------------
void ImplicitSurfaceDistance::Initialize()
{
  // Initialize base class
  ImplicitSurfaceForce::Initialize();

  // Initialize point data array of implicit surface distances
  InitializeDistances();

  // Add point data array of force magnitude
  AddPointData("Magnitude", 1, VTK_FLOAT);
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void ImplicitSurfaceDistance::Update(bool gradient)
{
  // Update base class
  ImplicitSurfaceForce::Update(gradient);

  // Calculate distances of surface points to the implicit target surface
  UpdateDistances();

  // Calculate force magnitude at surface points
  UpdateMagnitude();
}

// -----------------------------------------------------------------------------
void ImplicitSurfaceDistance::UpdateMagnitude()
{
  vtkDataArray * const status    = Status();
  vtkDataArray * const distances = Distances();
  vtkDataArray * const magnitude = PointData("Magnitude");

  // Evaluate/compute magnitude function at active surface points
  if (_MagnitudeImage) {

    ImageFunction func;
    func.Input(_MagnitudeImage);
    func.Initialize();
    EvaluateMagnitude eval;
    eval._Points    = Points();
    eval._Status    = status;
    eval._Function  = &func;
    eval._Distances = distances;
    eval._Magnitude = magnitude;
    parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);

  } else {

    double threshold = _MaxThreshold;
    if (threshold <= 0.) {
      if (_DistanceMeasure == DM_Minimum) {
        threshold = _MinThreshold;
      } else {
        threshold = AbsPercentile::Calculate(95, distances);
      }
    }

    ComputeMagnitude calc;
    calc._Status    = status;
    calc._Distances = distances;
    calc._Threshold = max(_MinThreshold, threshold);
    calc._Magnitude = magnitude;
    parallel_for(blocked_range<int>(0, _NumberOfPoints), calc);

  }

  // Invert magnitude
  if (_InvertMagnitude) {
    double m;
    double range[2] = { inf, 0. };
    for (vtkIdType ptId = 0; ptId < magnitude->GetNumberOfTuples(); ++ptId) {
      m = abs(magnitude->GetComponent(ptId, 0));
      if (m != 0.) {
        range[0] = min(range[0], m);
        range[1] = max(range[1], m);
      }
    }
    double extent = range[1] - range[0];
    for (vtkIdType ptId = 0; ptId < magnitude->GetNumberOfTuples(); ++ptId) {
      m = magnitude->GetComponent(ptId, 0);
      if (m != 0.) {
        m = copysign(extent - (abs(m) - range[0]) + range[0], m);
        magnitude->SetComponent(ptId, 0, m);
      }
    }
  }

  // Divide by maximum absolute magnitude value
  if (_NormalizeMagnitude) {
    double range[2];
    magnitude->GetRange(range, 0);
    double norm = max(abs(range[0]), abs(range[1]));
    if (norm > 0.) {
      norm = 1. / norm;
      for (vtkIdType ptId = 0; ptId < magnitude->GetNumberOfTuples(); ++ptId) {
        magnitude->SetComponent(ptId, 0, norm * magnitude->GetComponent(ptId, 0));
      }
    }
  }
}

// -----------------------------------------------------------------------------
double ImplicitSurfaceDistance::Evaluate()
{
  if (_NumberOfPoints == 0) return 0.;
  double d, sum = 0.;
  vtkDataArray *distances = this->Distances();
  for (int ptId = 0; ptId < _NumberOfPoints; ++ptId) {
    d = distances->GetComponent(ptId, 0);
    sum += abs(d);
  }
  return sum / _NumberOfPoints;
}

// -----------------------------------------------------------------------------
void ImplicitSurfaceDistance::EvaluateGradient(double *gradient, double step, double weight)
{
  if (_NumberOfPoints == 0) return;

  memset(_Gradient, 0, _NumberOfPoints * sizeof(GradientType));

  ComputeGradient eval;
  eval._Normals   = Normals();
  eval._Magnitude = PointData("Magnitude");
  eval._Gradient  = _Gradient;
  parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);

  ImplicitSurfaceForce::EvaluateGradient(gradient, step, weight / _NumberOfPoints);
}


} // namespace mirtk
