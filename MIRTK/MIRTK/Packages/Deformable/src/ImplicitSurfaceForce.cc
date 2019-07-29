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

#include "mirtk/ImplicitSurfaceForce.h"

#include "mirtk/Math.h"
#include "mirtk/Parallel.h"
#include "mirtk/List.h"
#include "mirtk/Stack.h"
#include "mirtk/UnorderedSet.h"
#include "mirtk/ImplicitSurfaceUtils.h"
#include "mirtk/MeshSmoothing.h"
#include "mirtk/DataStatistics.h"

#include "vtkType.h"
#include "vtkPoints.h"
#include "vtkDataArray.h"
#include "vtkPointData.h"


namespace mirtk {


// =============================================================================
// Auxiliary functions
// =============================================================================

namespace ImplicitSurfaceForceUtils {

// -----------------------------------------------------------------------------
/// Evaluate implicit surface distance function at mesh points
struct ComputeMinimumDistances
{
  ImplicitSurfaceForce *_Force;
  vtkPoints            *_Points;
  vtkDataArray         *_Status;
  vtkDataArray         *_Distances;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    double p[3];
    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      if (_Status && _Status->GetComponent(ptId, 0) == 0.) {
        _Distances->SetComponent(ptId, 0, 0.);
      } else {
        _Points->GetPoint(ptId, p);
        _Distances->SetComponent(ptId, 0, _Force->Distance(p));
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Compute distance to implicit surface along normal direction
struct ComputeNormalDistances
{
  ImplicitSurfaceForce *_Force;
  vtkPoints            *_Points;
  vtkDataArray         *_Status;
  vtkDataArray         *_Normals;
  vtkDataArray         *_Distances;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    double p[3], n[3];
    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      if (_Status && _Status->GetComponent(ptId, 0) == 0.) {
        _Distances->SetComponent(ptId, 0, 0.);
      } else {
        _Points ->GetPoint(ptId, p);
        _Normals->GetTuple(ptId, n);
        _Distances->SetComponent(ptId, 0, _Force->Distance(p, n));
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Find holes in implicit surface
int FindHoles(vtkPoints *points, const EdgeTable *edges, vtkDataArray *normals,
              vtkDataArray *mask, vtkDataArray *distances, vtkDataArray *labels,
              double d_min, double d_threshold, double max_radius, int max_size)
{
  const int    npoints     = points->GetNumberOfPoints();
  const double max_radius2 = max_radius * max_radius;

  Stack<int>        active;
  UnorderedSet<int> cluster, boundary;

  Point      p, c;
  Vector3    n1, n2, n;
  const int *adjIds;
  int        adjPts, ptId, nholes = 0;
  double     d0, distance, d_boundary;
  double     radius2, min_hole_radius2, max_hole_radius2, dp;
  bool       discard;

  labels->FillComponent(0, -1.);

  Array<double> dists(npoints);
  for (ptId = 0; ptId < npoints; ++ptId) {
    dists[ptId] = abs(distances->GetComponent(ptId, 0));
  }
  Array<int> order = DecreasingOrder(dists);
  for (auto seedId : order) {
    if (labels->GetComponent(seedId, 0) >= 0.) continue;
    d0 = abs(distances->GetComponent(seedId, 0));
    if (d0 < d_min) break;
    d_boundary = max(d_threshold, d0 / 3.);
    cluster.clear();
    boundary.clear();
    active.push(seedId);
    while (!active.empty()) {
      ptId = active.top(), active.pop();
      if (labels->GetComponent(ptId, 0) < 0.) {
        cluster.insert(ptId);
        normals->GetTuple(ptId, n1);
        labels->SetComponent(ptId, 0, 2.);
        edges->GetAdjacentPoints(ptId, adjPts, adjIds);
        for (int i = 0; i < adjPts; ++i) {
          auto &adjId = adjIds[i];
          if (labels->GetComponent(adjId, 0) < 0.) {
            distance = abs(distances->GetComponent(adjId, 0));
            if (distance < d_boundary) {
              boundary.insert(adjId);
            } else {
              normals->GetTuple(adjId, n2);
              if (n1.Dot(n2) < .2) {
                boundary.insert(adjId);
              } else {
                active.push(adjId);
              }
            }
          }
        }
      }
    }
    // Discard clusters containing masked points
    discard = false;
    if (mask) {
      for (auto ptId : cluster) {
        if (mask->GetComponent(ptId, 0) == 0.) {
          discard = true;
          break;
        }
      }
    }
    // Discard large clusters of distant points
    if (!discard) {
      discard = (cluster.size() > static_cast<size_t>(max_size) || boundary.empty());
    }
    if (!discard) {
      c = 0.;
      for (auto ptId : boundary) {
        points->GetPoint(ptId, p);
        c += p;
      }
      c /= boundary.size();
      min_hole_radius2 = inf;
      max_hole_radius2 = 0.;
      for (auto ptId : boundary) {
        points->GetPoint(ptId, p);
        radius2 = p.SquaredDistance(c);
        if (radius2 < min_hole_radius2) {
          min_hole_radius2 = radius2;
        }
        if (radius2 > max_hole_radius2) {
          max_hole_radius2 = radius2;
          if (radius2 > max_radius2) break;
        }
      }
      discard = (max_hole_radius2 > max_radius2 || 9. * min_hole_radius2 < max_hole_radius2);
    }
    // Discard clusters which are caused for example by bridges in a sulcus which
    // are yet to be deflated by the convex hull/sphere to white surface mesh deformation.
    if (!discard) {
      bool check_edges = false;
      normals->GetTuple(seedId, n1);
      for (auto ptId : cluster) {
        normals->GetTuple(ptId, n2);
        if (n1.Dot(n2) < -.5) {
          check_edges = true;
          break;
        }
      }
      if (check_edges) {
        size_t num_edge_points = 0;
        for (auto ptId : boundary) {
          n1 = 0., n2 = 0.;
          edges->GetAdjacentPoints(ptId, adjPts, adjIds);
          for (int i = 0; i < adjPts; ++i) {
            auto &adjId = adjIds[i];
            normals->GetTuple(adjId, n);
            if (labels->GetComponent(adjId, 0) == 2.) {
              n1 += n;
            } else if (boundary.find(adjId) == boundary.end()) {
              n2 += n;
            }
          }
          n1.Normalize();
          n2.Normalize();
          dp = n1.Dot(n2);
          if (abs(dp) < .5) ++num_edge_points;
        }
        discard = (num_edge_points > max(.1 * boundary.size(), 4.));
      }
    }
    // When clusters should be discarded, change label to zero and one otherwise
    if (discard) {
      for (auto ptId : cluster) {
        labels->SetComponent(ptId, 0, 0.);
      }
    } else {
      for (auto ptId : cluster) {
        labels->SetComponent(ptId, 0, 1.);
      }
      ++nholes;
    }
  }
  // Ensure all points have a non-negative label
  for (ptId = 0; ptId < npoints; ++ptId) {
    if (labels->GetComponent(ptId, 0) < 0.) {
      labels->SetComponent(ptId, 0, 0.);
    }
  }

  return nholes;
}

// -----------------------------------------------------------------------------
/// Dilate holes found in implicit surface
void DilateHoles(const EdgeTable *edges, vtkDataArray *distances, vtkDataArray *labels, int niter)
{
  if (niter < 1) return;
  const int npoints = labels->GetNumberOfTuples();

  int        adjPts;
  const int *adjIds;
  double     distance;

  UnorderedSet<int> boundary;
  for (int iter = 0; iter < niter; ++iter) {
    for (int ptId = 0; ptId < npoints; ++ptId) {
      if (labels->GetComponent(ptId, 0) != 0.) {
        distance = distances->GetComponent(ptId, 0);
        edges->GetAdjacentPoints(ptId, adjPts, adjIds);
        for (int i = 0; i < adjPts; ++i) {
          if (labels->GetComponent(adjIds[i], 0) == 0. &&
              distances->GetComponent(adjIds[i], 0) < distance) {
            boundary.insert(adjIds[i]);
            break;
          }
        }
      }
    }
    for (auto ptId : boundary) {
      labels->SetComponent(ptId, 0, 1.);
    }
    boundary.clear();
  }
}

// -----------------------------------------------------------------------------
/// Replace surface distance measurement of holes by average distances of hole boundary points
void FixHoles(vtkPoints *points, const EdgeTable *edges, vtkDataArray *normals,
              vtkDataArray *distances, vtkDataArray *labels, vtkDataArray *output, int niter)
{
  const int npoints = points->GetNumberOfPoints();

  int        adjPts;
  const int *adjIds;
  double     w, wsum, sum;
  Vector3    n1, n2;

  UnorderedSet<int> ptIds, active, boundary;
  for (int ptId = 0; ptId < npoints; ++ptId) {
    if (labels->GetComponent(ptId, 0) != 0.) {
      ptIds.insert(ptId);
      edges->GetAdjacentPoints(ptId, adjPts, adjIds);
      for (int i = 0; i < adjPts; ++i) {
        if (labels->GetComponent(adjIds[i], 0) == 0.) {
          boundary.insert(ptId);
          break;
        }
      }
    }
  }

  if (output != distances) {
    for (int ptId = 0; ptId < npoints; ++ptId) {
      output->SetComponent(ptId, 0, distances->GetComponent(ptId, 0));
    }
  }

  if (ptIds.empty()) return;

  while (!boundary.empty()) {
    active.swap(boundary);
    boundary.clear();
    for (auto ptId : active) {
      sum = wsum = 0.;
      normals->GetTuple(ptId, n1);
      edges->GetAdjacentPoints(ptId, adjPts, adjIds);
      for (int i = 0; i < adjPts; ++i) {
        if (ptIds.find(adjIds[i]) == ptIds.end()) {
          normals->GetTuple(adjIds[i], n2);
          w = clamp(n1.Dot(n2), 0., 1.);
          sum  += w * distances->GetComponent(adjIds[i], 0);
          wsum += w;
        } else {
          boundary.insert(adjIds[i]);
        }
      }
      if (wsum > 0.) sum /= wsum;
      output->SetComponent(ptId, 0, sum);
      ptIds.erase(ptId);
    }
  }

  vtkSmartPointer<vtkDataArray> input;
  input.TakeReference(output->NewInstance());
  for (int iter = 0; iter < niter; ++iter) {
    input->DeepCopy(output);
    for (int ptId = 0; ptId < npoints; ++ptId) {
      if (labels->GetComponent(ptId, 0) != 0.) {
        sum = wsum = 0.;
        normals->GetTuple(ptId, n1);
        edges->GetAdjacentPoints(ptId, adjPts, adjIds);
        for (int i = 0; i < adjPts; ++i) {
          normals->GetTuple(adjIds[i], n2);
          w = clamp(n1.Dot(n2), 0., 1.);
          sum  += w * input->GetComponent(adjIds[i], 0);
          wsum += w;
        }
        if (wsum > 0.) sum /= wsum;
        output->SetComponent(ptId, 0, sum);
      }
    }
  }
}


} // namespace ImplicitSurfaceForceUtils
using namespace ImplicitSurfaceForceUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
ImplicitSurfaceForce::ImplicitSurfaceForce(const char *name, double weight)
:
  SurfaceForce(name, weight),
  _DistanceMeasure(DM_Minimum),
  _Offset(0.),
  _MinStepLength(.1),
  _MaxDistance(0.),
  _Tolerance(1e-3),
  _DistanceSmoothing(1),
  _FillInHoles(false)
{
}

// -----------------------------------------------------------------------------
void ImplicitSurfaceForce::CopyAttributes(const ImplicitSurfaceForce &other)
{
  _DistanceMeasure   = other._DistanceMeasure;
  _Offset            = other._Offset;
  _MinStepLength     = other._MinStepLength;
  _MaxDistance       = other._MaxDistance;
  _Tolerance         = other._Tolerance;
  _DistanceSmoothing = other._DistanceSmoothing;
  _FillInHoles       = other._FillInHoles;
}

// -----------------------------------------------------------------------------
ImplicitSurfaceForce::ImplicitSurfaceForce(const ImplicitSurfaceForce &other)
:
  SurfaceForce(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
ImplicitSurfaceForce &ImplicitSurfaceForce::operator =(const ImplicitSurfaceForce &other)
{
  if (this != &other) {
    SurfaceForce::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
ImplicitSurfaceForce::~ImplicitSurfaceForce()
{
}

// =============================================================================
// Configuration
// =============================================================================

// -----------------------------------------------------------------------------
bool ImplicitSurfaceForce::SetWithPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Implicit surface distance measure") == 0) {
    return FromString(value, _DistanceMeasure);
  }
  if (strcmp(param, "Implicit surface distance offset") == 0) {
    return FromString(value, _Offset);
  }
  if (strcmp(param, "Implicit surface distance step length") == 0) {
    return FromString(value, _MinStepLength);
  }
  if (strcmp(param, "Implicit surface distance threshold") == 0) {
    return FromString(value, _MaxDistance);
  }
  if (strcmp(param, "Implicit surface distance tolerance") == 0) {
    return FromString(value, _Tolerance);
  }
  if (strcmp(param, "Implicit surface distance smoothing") == 0) {
    return FromString(value, _DistanceSmoothing);
  }
  if (strcmp(param, "Implicit surface distance hole filling") == 0) {
    return FromString(value, _FillInHoles);
  }
  return SurfaceForce::SetWithPrefix(param, value);
}

// -----------------------------------------------------------------------------
bool ImplicitSurfaceForce::SetWithoutPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Measure") == 0) {
    return FromString(value, _DistanceMeasure);
  }
  if (strcmp(param, "Offset") == 0) {
    return FromString(value, _Offset);
  }
  if (strcmp(param, "Step length") == 0) {
    return FromString(value, _MinStepLength);
  }
  if (strcmp(param, "Threshold") == 0) {
    return FromString(value, _MaxDistance);
  }
  if (strcmp(param, "Tolerance") == 0) {
    return FromString(value, _Tolerance);
  }
  if (strcmp(param, "Smoothing") == 0) {
    return FromString(value, _DistanceSmoothing);
  }
  if (strcmp(param, "Hole filling") == 0) {
    return FromString(value, _FillInHoles);
  }
  return SurfaceForce::SetWithoutPrefix(param, value);
}

// -----------------------------------------------------------------------------
ParameterList ImplicitSurfaceForce::Parameter() const
{
  ParameterList params = SurfaceForce::Parameter();
  InsertWithPrefix(params, "Measure",      _DistanceMeasure);
  InsertWithPrefix(params, "Offset",       _Offset);
  InsertWithPrefix(params, "Step length",  _MinStepLength);
  InsertWithPrefix(params, "Threshold",    _MaxDistance);
  InsertWithPrefix(params, "Tolerance",    _Tolerance);
  InsertWithPrefix(params, "Smoothing",    _DistanceSmoothing);
  InsertWithPrefix(params, "Hole filling", _FillInHoles);
  return params;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void ImplicitSurfaceForce::Initialize()
{
  // Initialize base class
  SurfaceForce::Initialize();

  // Initialize maximum distance
  if (_MaxDistance <= .0) {
    VoxelType mind, maxd;
    _Image->GetMinMax(mind, maxd); // ignoring background
    _MaxDistance = static_cast<double>(max(abs(mind), abs(maxd)));
  }

  // Initialize input image interpolators
  _Distance.Input(_Image);
  _Distance.Initialize();
  _DistanceGradient.Input(_Image);
  _DistanceGradient.Initialize();
}

// =============================================================================
// Surface distance
// =============================================================================

// -----------------------------------------------------------------------------
double ImplicitSurfaceForce::SelfDistance(const double p[3], const double n[3]) const
{
  return SurfaceForce::SelfDistance(p, n, _MaxDistance);
}

// -----------------------------------------------------------------------------
double ImplicitSurfaceForce::Distance(const double p[3]) const
{
  return ImplicitSurfaceUtils::Evaluate(_Distance, p, _Offset);
}

// -----------------------------------------------------------------------------
double ImplicitSurfaceForce::Distance(const double p[3], const double n[3]) const
{
  const double mind = ImplicitSurfaceUtils::Evaluate(_Distance, p, _Offset);
  return ImplicitSurfaceUtils::SignedDistance(p, n, mind, _MinStepLength, _MaxDistance, _Distance, _Offset, _Tolerance);
}

// -----------------------------------------------------------------------------
void ImplicitSurfaceForce::DistanceGradient(const double p[3], double g[3], bool normalize) const
{
  ImplicitSurfaceUtils::Evaluate(_DistanceGradient, p, g, normalize);
}

// -----------------------------------------------------------------------------
vtkDataArray *ImplicitSurfaceForce::MinimumDistances() const
{
  return PointData("MinimumImplicitSurfaceDistance");
}

// -----------------------------------------------------------------------------
void ImplicitSurfaceForce::InitializeMinimumDistances()
{
  vtkDataArray *d = AddPointData("MinimumImplicitSurfaceDistance", 1, VTK_FLOAT, true);
  d->FillComponent(0, numeric_limits<double>::infinity());
}

// -----------------------------------------------------------------------------
void ImplicitSurfaceForce::UpdateMinimumDistances()
{
  vtkDataArray * const distances = MinimumDistances();
  if (distances->GetMTime() < _PointSet->Surface()->GetMTime()) {
    ComputeMinimumDistances eval;
    eval._Force     = this;
    eval._Points    = Points();
    eval._Status    = InitialStatus();
    eval._Distances = distances;
    parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);
  }
}

// -----------------------------------------------------------------------------
vtkDataArray *ImplicitSurfaceForce::NormalDistances() const
{
  return PointData("NormalImplicitSurfaceDistance");
}

// -----------------------------------------------------------------------------
void ImplicitSurfaceForce::InitializeNormalDistances()
{
  vtkDataArray *d = AddPointData("NormalImplicitSurfaceDistance", 1, VTK_FLOAT, true);
  if (_FillInHoles) {
    AddPointData("ImplicitSurfaceHoleMask", 1, VTK_CHAR, true);
  }
  if (debug && (_FillInHoles || _DistanceSmoothing > 0)) {
    AddPointData("OriginalNormalImplicitSurfaceDistance", 1, VTK_FLOAT, true);
  }
  d->FillComponent(0, _MaxDistance);
}

// -----------------------------------------------------------------------------
void ImplicitSurfaceForce::UpdateNormalDistances()
{
  vtkPolyData  * const surface   = DeformedSurface();
  vtkDataArray * const distances = NormalDistances();
  if (distances->GetMTime() < surface->GetMTime()) {

    vtkDataArray * const status = InitialStatus();

    ComputeNormalDistances eval;
    eval._Force     = this;
    eval._Points    = Points();
    eval._Status    = status;
    eval._Normals   = Normals();
    eval._Distances = distances;
    parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);

    if (debug) {
      vtkDataArray * const orig = PointData("OriginalNormalImplicitSurfaceDistance", true);
      if (orig) orig->CopyComponent(0, distances, 0);
    }
    if (_FillInHoles) {
      int    num = 0;
      double mean = 0., sigma = 0., d;
      for (int ptId = 0; ptId < _NumberOfPoints; ++ptId) {
        if (!status || status->GetComponent(ptId, 0) != 0.) {
          d = abs(distances->GetComponent(ptId, 0));
          if (d < _MaxDistance) {
            mean  += d;
            sigma += d * d;
            ++num;
          }
        }
      }
      if (num > 0) {
        mean /= num;
        sigma = sqrt(sigma / num - mean * mean);
        const double d_min = max(mean + 3. * sigma, .25 * _MaxDistance);
        const double d_threshold = mean + sigma;
        if (d_threshold + sigma < d_min) {
          const double edge_length = _PointSet->AverageInputSurfaceEdgeLength();
          const double max_radius  = 5. * edge_length;
          const size_t max_size    = 100;
          const bool   optional    = true;
          vtkDataArray * const mask  = PointData("ImplicitSurfaceFillMask", optional);
          vtkDataArray * const holes = PointData("ImplicitSurfaceHoleMask");
          if (FindHoles(Points(), Edges(), Normals(), mask, distances, holes, d_min, d_threshold, max_radius, max_size) > 0) {
            DilateHoles(Edges(), distances, holes, 2);
            FixHoles(Points(), Edges(), Normals(), distances, holes, distances, 3);
          }
        }
      }
    }
    if (_DistanceSmoothing > 0) {
      MeshSmoothing smoother;
      smoother.Input(surface);
      smoother.Mask(status);
      smoother.SmoothPointsOff();
      smoother.SmoothArray(distances->GetName());
      smoother.Weighting(MeshSmoothing::NormalDeviation);
      smoother.NumberOfIterations(_DistanceSmoothing);
      smoother.Run();
      vtkPointData *outputPD = smoother.Output()->GetPointData();
      distances->CopyComponent(0, outputPD->GetArray(distances->GetName()), 0);
    }

    distances->Modified();
  }
}

// -----------------------------------------------------------------------------
vtkDataArray *ImplicitSurfaceForce::Distances() const
{
  switch (_DistanceMeasure) {
    case DM_Minimum: return MinimumDistances();
    case DM_Normal:  return NormalDistances();
    default:
      cerr << "ImplicitSurfaceForce::Distances: Invalid distance measure type: " << _DistanceMeasure << endl;
      exit(1);
  }
}

// -----------------------------------------------------------------------------
void ImplicitSurfaceForce::InitializeDistances()
{
  switch (_DistanceMeasure) {
    case DM_Minimum: InitializeMinimumDistances(); break;
    case DM_Normal:  InitializeNormalDistances();  break;
    default:
      cerr << "ImplicitSurfaceForce::InitializeDistances: Invalid distance measure type: " << _DistanceMeasure << endl;
      exit(1);
  }
}

// -----------------------------------------------------------------------------
void ImplicitSurfaceForce::UpdateDistances()
{
  switch (_DistanceMeasure) {
    case DM_Minimum: UpdateMinimumDistances(); break;
    case DM_Normal:  UpdateNormalDistances();  break;
    default:
      cerr << "ImplicitSurfaceForce::UpdateDistances: Invalid distance measure type: " << _DistanceMeasure << endl;
      exit(1);
  }
}


} // namespace mirtk
