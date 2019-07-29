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

#include "mirtk/RepulsiveForce.h"

#include "mirtk/Math.h"
#include "mirtk/Point.h"
#include "mirtk/Vector3.h"
#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/PointSetUtils.h"

#include "vtkSmartPointer.h"
#include "vtkIdList.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkAbstractPointLocator.h"
#include "vtkOctreePointLocator.h"


namespace mirtk {


// Global debug flag defined in Options.cc
MIRTK_Common_EXPORT extern int debug;

// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(RepulsiveForce);


// =============================================================================
// Auxiliaries
// =============================================================================

namespace RepulsiveForceUtils {

// -----------------------------------------------------------------------------
inline double Repulsion(double d, double r)
{
  if (d > r) return 0.;
  d = r - d;
  return (d * d * d) / (r * r);
}

// -----------------------------------------------------------------------------
inline double RepulsionDerivative(double d, double r)
{
  if (d > r) return 0.;
  d = r - d;
  return 3. * (d * d) / (r * r);
}

// -----------------------------------------------------------------------------
/// Evaluate energy of repulsive force term
struct Evaluate
{
  vtkPoints               *_Points;
  vtkDataArray            *_Normals;
  vtkAbstractPointLocator *_Locator;
  double                   _FrontfaceDistance;
  double                   _BackfaceDistance;
  double                   _Penalty;

  Evaluate() : _Penalty(0.) {}

  Evaluate(const Evaluate &other, split)
  :
    _Points(other._Points),
    _Normals(other._Normals),
    _Locator(other._Locator),
    _FrontfaceDistance(other._FrontfaceDistance),
    _BackfaceDistance(other._BackfaceDistance),
    _Penalty(0.)
  {}

  void join(const Evaluate &other)
  {
    _Penalty += other._Penalty;
  }

  void operator ()(const blocked_range<int> &ptIds)
  {
    int     num;
    double  d, w, m, penalty;
    Point   p, q;
    Vector3 f, np, nq;

    const double r = max(_FrontfaceDistance, _BackfaceDistance);

    vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
    for (int ptId = ptIds.begin(), id; ptId != ptIds.end(); ++ptId) {
      _Points ->GetPoint(ptId, p);
      _Normals->GetTuple(ptId, np);
      _Locator->FindPointsWithinRadius(r, p, ids);
      if (ids->GetNumberOfIds() > 0) {
        penalty = 0., num = 0;
        for (vtkIdType i = 0; i < ids->GetNumberOfIds(); ++i) {
          id = static_cast<int>(ids->GetId(i));
          if (id != ptId) {
            _Points->GetPoint(id, q);
            _Normals->GetTuple(id, nq);
            f = q - p;
            d = f.Length();
            w = SShapedMembershipFunction(-np.Dot(nq), 0., 1.);
            if (np.Dot(f) < 0.) {
              m = w * Repulsion(d, _BackfaceDistance);
            } else {
              m = w * Repulsion(d, _FrontfaceDistance);
            }
            if (m != 0.) {
              penalty += m;
              ++num;
            }
          }
        }
        if (num > 0) _Penalty += penalty / num;
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate gradient of repulsive force term
struct EvaluateGradient
{
  typedef RepulsiveForce::GradientType GradientType;

  vtkPoints               *_Points;
  vtkDataArray            *_Status;
  vtkDataArray            *_Normals;
  vtkAbstractPointLocator *_Locator;
  double                   _FrontfaceDistance;
  double                   _BackfaceDistance;
  GradientType            *_Gradient;
  vtkDataArray            *_Magnitude;
  double                   _Weight;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    int          num;
    double       d, w, m, mag;
    Vector3      p, f, np, nq;
    GradientType gradient;

    const double r = max(_FrontfaceDistance, _BackfaceDistance);

    vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
    for (int ptId = ptIds.begin(), id; ptId != ptIds.end(); ++ptId) {
      if (_Status && _Status->GetComponent(ptId, 0) == .0) continue;
      _Points ->GetPoint(ptId, p);
      _Normals->GetTuple(ptId, np);
      _Locator->FindPointsWithinRadius(r, p, ids);
      if (ids->GetNumberOfIds() > 0) {
        gradient = 0., num = 0, mag = 0.;
        for (vtkIdType i = 0; i < ids->GetNumberOfIds(); ++i) {
          id = static_cast<int>(ids->GetId(i));
          if (id != ptId) {
            _Points->GetPoint(id, f), f -= p;
            _Normals->GetTuple(id, nq);
            d = f.Normalize();
            w = SShapedMembershipFunction(-np.Dot(nq), 0., 1.);
            if (np.Dot(f) < 0.) {
              m = w * RepulsionDerivative(d, _BackfaceDistance);
            } else {
              m = w * RepulsionDerivative(d, _FrontfaceDistance);
            }
            if (m > 0.) {
              mag += m;
              gradient._x += m * f._x;
              gradient._y += m * f._y;
              gradient._z += m * f._z;
              ++num;
            }
          }
        }
        if (num > 0) {
          _Gradient[ptId] = (gradient /= num);
          if (_Magnitude) _Magnitude->SetComponent(ptId, 0, _Weight * mag / num);
        } else {
          if (_Magnitude) _Magnitude->SetComponent(ptId, 0, 0.);
        }
      }
    }
  }
};


} // namespace RepulsiveForceUtils

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
RepulsiveForce::RepulsiveForce(const char *name, double weight)
:
  SurfaceConstraint(name, weight),
  _FrontfaceRadius(-1.),
  _BackfaceRadius(-1.)
{
  _ParameterPrefix.push_back("Repulsive force ");
  _ParameterPrefix.push_back("Node repulsion ");
}

// -----------------------------------------------------------------------------
void RepulsiveForce::CopyAttributes(const RepulsiveForce &other)
{
  _FrontfaceRadius = other._FrontfaceRadius;
  _BackfaceRadius  = other._BackfaceRadius;

  if (other._Locator) _Locator.TakeReference(other._Locator->NewInstance());
  else                _Locator = nullptr;
}

// -----------------------------------------------------------------------------
RepulsiveForce::RepulsiveForce(const RepulsiveForce &other)
:
  SurfaceConstraint(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
RepulsiveForce &RepulsiveForce::operator =(const RepulsiveForce &other)
{
  if (this != &other) {
    SurfaceConstraint::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
RepulsiveForce::~RepulsiveForce()
{
}

// -----------------------------------------------------------------------------
void RepulsiveForce::Initialize()
{
  _NumberOfPoints = 0;
  if (_FrontfaceRadius == .0 && _BackfaceRadius == 0.) return;

  SurfaceConstraint::Initialize();
  if (_FrontfaceRadius < 0. || _BackfaceRadius  < 0.) {
    const double edge_length = AverageEdgeLength(OriginalSurface());
    if (_FrontfaceRadius < 0.) {
      _FrontfaceRadius = abs(_FrontfaceRadius) * edge_length;
    }
    if (_BackfaceRadius < 0.) {
      _BackfaceRadius = abs(_BackfaceRadius) * edge_length;
    }
  }

  _Locator = vtkSmartPointer<vtkOctreePointLocator>::New();

  if (debug) AddPointData("Magnitude")->FillComponent(0, 0.);
}

// -----------------------------------------------------------------------------
void RepulsiveForce::Reinitialize()
{
  if (_FrontfaceRadius == .0 && _BackfaceRadius == 0.) {
    _NumberOfPoints = 0;
    _Locator        = nullptr;
  } else {
    SurfaceConstraint::Reinitialize();
  }
}

// =============================================================================
// Configuration
// =============================================================================

// -----------------------------------------------------------------------------
bool RepulsiveForce::SetWithoutPrefix(const char *name, const char *value)
{
  if (strcmp(name, "Radius") == 0) {
    double r;
    if (!FromString(value, r)) return false;
    _FrontfaceRadius = _BackfaceRadius = r;
    return true;
  }
  if (strcmp(name, "Frontface radius") == 0) {
    return FromString(value, _FrontfaceRadius);
  }
  if (strcmp(name, "Backface radius") == 0) {
    return FromString(value, _BackfaceRadius);
  }
  return PointSetForce::SetWithoutPrefix(name, value);
}

// -----------------------------------------------------------------------------
ParameterList RepulsiveForce::Parameter() const
{
  ParameterList params = PointSetForce::Parameter();
  if (fequal(_FrontfaceRadius, _BackfaceRadius)) {
    InsertWithPrefix(params, "Radius", _FrontfaceRadius);
  } else {
    InsertWithPrefix(params, "Frontface radius", _FrontfaceRadius);
    InsertWithPrefix(params, "Backface radius",  _BackfaceRadius);
  }
  return params;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void RepulsiveForce::Update(bool gradient)
{
  if (_NumberOfPoints == .0) return;

  // Update base class
  SurfaceConstraint::Update(gradient);

  // Make shallow copy without data arrays such that modified time of data set
  // does not depend on unused attributes that may be modified by other terms
  vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();
  surface->ShallowCopy(DeformedSurface());
  surface->GetPointData()->Initialize();
  surface->GetCellData ()->Initialize();
  _Locator->SetDataSet(surface);
  _Locator->BuildLocator();
}

// -----------------------------------------------------------------------------
double RepulsiveForce::Evaluate()
{
  if (_NumberOfPoints == 0) return 0.;

  RepulsiveForceUtils::Evaluate eval;
  eval._Points            = Points();
  eval._Normals           = Normals();
  eval._Locator           = _Locator;
  eval._FrontfaceDistance = _FrontfaceRadius;
  eval._BackfaceDistance  = _BackfaceRadius;
  parallel_reduce(blocked_range<int>(0, _NumberOfPoints), eval);

  return eval._Penalty / _NumberOfPoints;
}

// -----------------------------------------------------------------------------
void RepulsiveForce::EvaluateGradient(double *gradient, double step, double weight)
{
  if (_NumberOfPoints == 0) return;

  memset(_Gradient, 0, _NumberOfPoints * sizeof(GradientType));

  RepulsiveForceUtils::EvaluateGradient eval;
  eval._Points            = Points();
  eval._Status            = Status();
  eval._Normals           = Normals();
  eval._Locator           = _Locator;
  eval._FrontfaceDistance = _FrontfaceRadius;
  eval._BackfaceDistance  = _BackfaceRadius;
  eval._Gradient          = _Gradient;
  eval._Magnitude         = PointData("Magnitude", true);
  eval._Weight            = weight;
  parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);

  SurfaceConstraint::EvaluateGradient(gradient, step, weight / _NumberOfPoints);
}


} // namespace mirtk
