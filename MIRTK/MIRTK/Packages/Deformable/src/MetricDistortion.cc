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

#include "mirtk/MetricDistortion.h"

#include "mirtk/Array.h"
#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/ObjectFactory.h"
#include "mirtk/VtkMath.h"

#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(MetricDistortion);


// =============================================================================
// Auxiliaries
// =============================================================================

namespace MetricDistortionUtils {

// -----------------------------------------------------------------------------
typedef RegisteredPointSet::NodeNeighbors NodeNeighbors;
typedef MetricDistortion::DistancesArray  DistancesArray;

// -----------------------------------------------------------------------------
/// Compute initial distances between neighboring nodes
struct ComputeInitialDistances
{
  vtkPoints             *_InitialPoints;
  const NodeNeighbors   *_Neighbors;
  Array<DistancesArray> *_Distances;
  int                    _Radius;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    int        numNbrPts;
    const int *nbrPtIds;
    double     c[3], p[3];

    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      _InitialPoints->GetPoint(ptId, c);
      _Neighbors->GetConnectedPoints(ptId, numNbrPts, nbrPtIds, _Radius);
      DistancesArray &dists = (*_Distances)[ptId];
      dists.resize(numNbrPts);
      for (int i = 0; i < numNbrPts; ++i) {
        _InitialPoints->GetPoint(nbrPtIds[i], p);
        dists[i]._Distance0 = sqrt(vtkMath::Distance2BetweenPoints(c, p));
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Compute current distances between neighboring nodes
struct ComputeDistances
{
  vtkPoints             *_Points;
  const NodeNeighbors   *_Neighbors;
  Array<DistancesArray> *_Distances;
  int                    _Radius;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    int        numNbrPts;
    const int *nbrPtIds;
    double     c[3], p[3];

    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      _Points->GetPoint(ptId, c);
      _Neighbors->GetConnectedPoints(ptId, numNbrPts, nbrPtIds, _Radius);
      DistancesArray &dists = (*_Distances)[ptId];
      for (int i = 0; i < numNbrPts; ++i) {
        _Points->GetPoint(nbrPtIds[i], p);
        dists[i]._Distance = sqrt(vtkMath::Distance2BetweenPoints(c, p));
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate metric distortion
struct Evaluate
{
  vtkPoints                   *_Points;
  const NodeNeighbors         *_Neighbors;
  const Array<DistancesArray> *_Distances;
  int                          _Radius;
  double                       _Scale;
  double                       _Sum;

  Evaluate() : _Sum(.0) {}

  Evaluate(const Evaluate &other, split)
  :
    _Points   (other._Points),
    _Neighbors(other._Neighbors),
    _Distances(other._Distances),
    _Radius   (other._Radius),
    _Scale    (other._Scale),
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
    int        numNbrPts;
    const int *nbrPtIds;
    double     delta, sum;

    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      const DistancesArray &dists = (*_Distances)[ptId];
      _Neighbors->GetConnectedPoints(ptId, numNbrPts, nbrPtIds, _Radius);
      if (numNbrPts > 0) {
        sum = .0;
        for (int i = 0; i < numNbrPts; ++i) {
          delta = dists[i]._Distance - _Scale * dists[i]._Distance0;
          sum += delta * delta;
        }
        _Sum += (sum / numNbrPts);
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate gradient of metric distortion (i.e., negative force)
struct EvaluateGradient
{
  typedef MetricDistortion::GradientType GradientType;

  vtkPoints                   *_Points;
  vtkDataArray                *_Normals;
  vtkDataArray                *_Status;
  const NodeNeighbors         *_Neighbors;
  const Array<DistancesArray> *_Distances;
  GradientType                *_Gradient;
  int                          _Radius;
  double                       _Scale;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    int        numNbrPts;
    const int *nbrPtIds;
    double     p1[3], p2[3], n[3], e[3], f[3], s;

    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      if (_Status && _Status->GetComponent(ptId, 0) == .0) continue;
      _Points->GetPoint(ptId, p1);
      const DistancesArray &dists = (*_Distances)[ptId];
      _Neighbors->GetConnectedPoints(ptId, numNbrPts, nbrPtIds, _Radius);
      if (numNbrPts == 0) continue;
      f[0] = f[1] = f[2] = .0;
      for (int i = 0; i < numNbrPts; ++i) {
        _Points->GetPoint(nbrPtIds[i], p2);
        s = (dists[i]._Distance - _Scale * dists[i]._Distance0) / dists[i]._Distance;
        vtkMath::Subtract(p2, p1, e);
        vtkMath::MultiplyScalar(e, s);
        vtkMath::Add(f, e, f);
      }
      vtkMath::MultiplyScalar(f, 1.0 / numNbrPts);
      // Take out normal component (cf. FreeSurfer's mrisComputeDistanceTerm)
      _Normals->GetTuple(ptId, n);
      s = -vtkMath::Dot(n, f);
      vtkMath::MultiplyScalar(n, s);
      vtkMath::Add(f, n, f);
      _Gradient[ptId] = GradientType(-f[0], -f[1], -f[2]);
    }
  }
};


} // namespace MetricDistortionUtils

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
MetricDistortion::MetricDistortion(const char *name, double weight)
:
  SurfaceConstraint(name, weight),
  _Radius(-1), // use RegisteredPointSet::NeighborhoodRadius
  _InitialArea(.0)
{
  _ParameterPrefix.push_back("Distortion ");
  _ParameterPrefix.push_back("Metric distortion ");
  _ParameterPrefix.push_back("Surface distortion ");
  _ParameterPrefix.push_back("Mesh distortion ");
  _ParameterPrefix.push_back("Surface mesh distortion ");
}

// -----------------------------------------------------------------------------
void MetricDistortion::CopyAttributes(const MetricDistortion &other)
{
  _Radius      = other._Radius;
  _InitialArea = other._InitialArea;
  _Distances   = other._Distances;
}

// -----------------------------------------------------------------------------
MetricDistortion::MetricDistortion(const MetricDistortion &other)
:
  SurfaceConstraint(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
MetricDistortion &MetricDistortion::operator =(const MetricDistortion &other)
{
  if (this != &other) {
    SurfaceConstraint::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
MetricDistortion::~MetricDistortion()
{
}

// =============================================================================
// Configuration
// =============================================================================

// -----------------------------------------------------------------------------
bool MetricDistortion::SetWithoutPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Radius") == 0) {
    return FromString(value, _Radius);
  }
  return SurfaceConstraint::SetWithoutPrefix(param, value);
}

// -----------------------------------------------------------------------------
ParameterList MetricDistortion::Parameter() const
{
  ParameterList params = SurfaceConstraint::Parameter();
  InsertWithPrefix(params, "Radius", _Radius);
  return params;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void MetricDistortion::Initialize()
{
  // Initialize base class
  SurfaceConstraint::Initialize();

  // Initialize this class
  MetricDistortion::Init();
}

// -----------------------------------------------------------------------------
void MetricDistortion::Reinitialize()
{
  // Reinitialize base class
  SurfaceConstraint::Reinitialize();

  // Reinitialize this class
  MetricDistortion::Init();
}

// -----------------------------------------------------------------------------
void MetricDistortion::Init()
{
  if (_NumberOfPoints > 0) {
    _Distances.resize(_NumberOfPoints);
    MetricDistortionUtils::ComputeInitialDistances eval;
    vtkSmartPointer<vtkPoints> points = GetInitialPoints();
    eval._InitialPoints = points;
    eval._Neighbors     = _PointSet->SurfaceNeighbors(_Radius);
    eval._Distances     = &_Distances;
    eval._Radius        = _Radius;
    parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);
    vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();
    surface->ShallowCopy(_PointSet->InputSurface());
    surface->SetPoints(points);
    _InitialArea = Area(surface);
  } else {
    _Distances.clear();
    _InitialArea = .0;
  }
}

// -----------------------------------------------------------------------------
void MetricDistortion::Update(bool gradient)
{
  SurfaceConstraint::Update(gradient);
  MetricDistortionUtils::ComputeDistances eval;
  eval._Points    = _PointSet->SurfacePoints();
  eval._Neighbors = _PointSet->SurfaceNeighbors(_Radius);
  eval._Distances = &_Distances;
  eval._Radius    = _Radius;
  parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);
}

// -----------------------------------------------------------------------------
double MetricDistortion::Evaluate()
{
  if (_NumberOfPoints == 0) return .0;
  MIRTK_START_TIMING();
  MetricDistortionUtils::Evaluate eval;
  eval._Points    = _PointSet->SurfacePoints();
  eval._Neighbors = _PointSet->SurfaceNeighbors(_Radius);
  eval._Distances =  &_Distances;
  eval._Radius    = _Radius;
  eval._Scale     = 1.0 / sqrt(_InitialArea / _PointSet->SurfaceArea());
  parallel_reduce(blocked_range<int>(0, _NumberOfPoints), eval);
  MIRTK_DEBUG_TIMING(3, "evaluation of metric distortion");
  return eval._Sum / _NumberOfPoints;
}

// -----------------------------------------------------------------------------
void MetricDistortion::EvaluateGradient(double *gradient, double step, double weight)
{
  if (_NumberOfPoints == 0) return;

  MIRTK_START_TIMING();
  memset(_Gradient, 0, _NumberOfPoints * sizeof(GradientType));

  MetricDistortionUtils::EvaluateGradient eval;
  eval._Points    = _PointSet->SurfacePoints();
  eval._Normals   = _PointSet->SurfaceNormals();
  eval._Status    = _PointSet->SurfaceStatus();
  eval._Neighbors = _PointSet->SurfaceNeighbors(_Radius);
  eval._Distances = &_Distances;
  eval._Gradient  = _Gradient;
  eval._Radius    = _Radius;
  eval._Scale     = 1.0 / sqrt(_InitialArea / _PointSet->SurfaceArea());
  parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);

  InternalForce::EvaluateGradient(gradient, step, 2.0 * weight / _NumberOfPoints);
  MIRTK_DEBUG_TIMING(3, "evaluation of metric distortion force");
}


} // namespace mirtk
