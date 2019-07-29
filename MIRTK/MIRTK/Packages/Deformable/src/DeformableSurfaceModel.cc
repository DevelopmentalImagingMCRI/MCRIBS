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

#include "mirtk/DeformableSurfaceModel.h"

#include "mirtk/Config.h" // WINDOWS
#include "mirtk/Array.h"
#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/Triangle.h"

#include "mirtk/Vtk.h"
#include "mirtk/VtkMath.h"

#include "mirtk/MeshSmoothing.h"
#include "mirtk/SurfaceRemeshing.h"
#include "mirtk/SurfaceCurvature.h"
#include "mirtk/SurfaceCollisions.h"
#include "mirtk/PointSamples.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/PointSetIO.h"

#include "mirtk/ImplicitSurfaceForce.h"
#include "mirtk/ImplicitSurfaceUtils.h"

#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkUnsignedCharArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkPolyDataNormals.h"
#include "vtkAbstractCellLocator.h"
#include "vtkGenericCell.h"
#include "vtkCellLocator.h"
#include "vtkWindowedSincPolyDataFilter.h"


namespace mirtk {


using namespace ImplicitSurfaceUtils;


// =============================================================================
// Auxiliary functor classes for parallel execution
// =============================================================================

namespace DeformableSurfaceModelUtils {


// -----------------------------------------------------------------------------
/// Compute centroid of point set
void GetCentroid(vtkPolyData *mesh, double centroid[3])
{
  double p[3];
  centroid[0] = centroid[1] = centroid[2] = .0;
  const vtkIdType npoints = mesh->GetNumberOfPoints();
  for (vtkIdType ptId = 0; ptId < npoints; ++ptId) {
    mesh->GetPoint(ptId, p);
    centroid[0] += p[0];
    centroid[1] += p[1];
    centroid[2] += p[2];
  }
  centroid[0] /= npoints;
  centroid[1] /= npoints;
  centroid[2] /= npoints;
}

// -----------------------------------------------------------------------------
/// Get approximate scale of point set
void GetScale(vtkPolyData *mesh, const double centroid[3], double scale[3])
{
  double p[3];
  scale[0] = scale[1] = scale[2] = .0;
  const vtkIdType npoints = mesh->GetNumberOfPoints();
  for (vtkIdType ptId = 0; ptId < npoints; ++ptId) {
    mesh->GetPoint(ptId, p);
    scale[0] += abs(p[0] - centroid[0]);
    scale[1] += abs(p[1] - centroid[1]);
    scale[2] += abs(p[2] - centroid[2]);
  }
  scale[0] /= npoints;
  scale[1] /= npoints;
  scale[2] /= npoints;
}

// -----------------------------------------------------------------------------
/// Determine maximum vertex displacement
class MaxVertexDisplacement
{
private:

  const double *_Gradient;
  double        _MaxNorm;
  bool          _FoundNaN;

public:

  /// Constructor
  MaxVertexDisplacement(const double *gradient)
  :
    _Gradient(gradient), _MaxNorm(.0), _FoundNaN(false)
  {}

  /// Copy constructor
  MaxVertexDisplacement(const MaxVertexDisplacement &other)
  :
    _Gradient(other._Gradient),
    _MaxNorm (other._MaxNorm),
    _FoundNaN(other._FoundNaN)
  {}

  /// Split constructor
  MaxVertexDisplacement(const MaxVertexDisplacement &other, split)
  :
    _Gradient(other._Gradient),
    _MaxNorm (other._MaxNorm),
    _FoundNaN(other._FoundNaN)
  {}

  /// Join results
  void join(const MaxVertexDisplacement &other)
  {
    if (other._MaxNorm > _MaxNorm) _MaxNorm = other._MaxNorm;
    if (other._FoundNaN) _FoundNaN = true;
  }

  /// Maximum norm
  double Norm() const { return sqrt(_MaxNorm); }

  /// Whether gradient vector contains at least one NaN value
  bool FoundNaN() const { return _FoundNaN; }

  /// Determine maximum norm of specified vertex displacements
  void operator()(const blocked_range<int> &re)
  {
    double norm;
    const double *g = _Gradient + 3 * re.begin();
    for (int i = re.begin(); i != re.end(); ++i, g += 3) {
      norm = pow(g[0], 2) + pow(g[1], 2) + pow(g[2], 2);
      if (norm > _MaxNorm) _MaxNorm = norm;
      if (IsNaN(norm)) _FoundNaN = true;
    }
  }
};

// -----------------------------------------------------------------------------
/// Move points of deformable surface model
struct MovePoints
{
  vtkPoints    *_InputPoints;
  vtkPoints    *_OutputPoints;
  const double *_Displacement;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    double x[3];
    const double *dx = _Displacement + 3 * re.begin();
    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId, dx += 3) {
      _InputPoints->GetPoint(ptId, x);
      x[0] += dx[0], x[1] += dx[1], x[2] += dx[2];
      _OutputPoints->SetPoint(ptId, x);
    }
  }

  static void Run(vtkPoints *input, const double *dx, vtkPoints *output = NULL)
  {
    MovePoints move;
    move._InputPoints  = input;
    move._OutputPoints = output ? output : input;
    move._Displacement = dx;
    parallel_for(blocked_range<vtkIdType>(0, input->GetNumberOfPoints()), move);
  }
};

// -----------------------------------------------------------------------------
/// Perform one iteration of gradient averaging
struct AverageGradient
{
  const EdgeTable *_EdgeTable;
  double          *_Input;
  double          *_Output;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    int n;
    const int *adjPts;
    double *in, *out = _Output + 3 * re.begin();

    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId, out += 3) {
      in = _Input + 3 * ptId;
      out[0] = in[0], out[1] = in[1], out[2] = in[2];
      _EdgeTable->GetAdjacentPoints(ptId, n, adjPts);
      for (int i = 0; i < n; ++i) {
        in = _Input + 3 * adjPts[i];
        out[0] += in[0], out[1] += in[1], out[2] += in[2];
      }
      n += 1;
      out[0] /= n, out[1] /= n, out[2] /= n;
    }
  }
};

// -----------------------------------------------------------------------------
/// Perform one iteration of gradient magnitude averaging
struct AverageGradientMagnitude
{
  const EdgeTable *_EdgeTable;
  double          *_Input;
  double          *_Output;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    int        numAdjPts;
    const int *adjPts;
    double     norm, avg_norm, *adj;

    double *in  = _Input  + 3 * re.begin();
    double *out = _Output + 3 * re.begin();

    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId, in += 3, out += 3) {
      norm = sqrt(vtkMath::Dot(in, in));
      if (norm) {
        avg_norm = norm;
        _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPts);
        for (int i = 0; i < numAdjPts; ++i) {
          adj = _Input + 3 * adjPts[i];
          avg_norm += sqrt(vtkMath::Dot(adj, adj));
        }
        avg_norm /= (numAdjPts + 1);
        avg_norm /= norm;
        out[0] = in[0] * avg_norm;
        out[1] = in[1] * avg_norm;
        out[2] = in[2] * avg_norm;
      } else {
        out[0] = out[1] = out[2] = .0;
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Perform one iteration of signed gradient averaging
struct AverageSignedGradient
{
  const EdgeTable *_EdgeTable;
  double          *_Input;
  double          *_Output;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    int        numAdjPts, n;
    const int *adjPts;
    double    *adj;

    double *in  = _Input  + 3 * re.begin();
    double *out = _Output + 3 * re.begin();

    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId, in += 3, out += 3) {
      out[0] = in[0], out[1] = in[1], out[2] = in[2], n = 1;
      _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPts);
      for (int i = 0; i < numAdjPts; ++i) {
        adj = _Input + 3 * adjPts[i];
        if (adj[0]*in[0] + adj[1]*in[1] + adj[2]*in[2] > .0) {
          out[0] += adj[0], out[1] += adj[1], out[2] += adj[2], ++n;
        }
      }
      out[0] /= n, out[1] /= n, out[2] /= n;
    }
  }
};

// -----------------------------------------------------------------------------
/// Perform one iteration of signed gradient magnitude averaging
struct AverageSignedGradientMagnitude
{
  const EdgeTable *_EdgeTable;
  double          *_Input;
  double          *_Output;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    int        numAdjPts, n;
    const int *adjPts;
    double    *adj, norm, avg_norm;

    double *in  = _Input  + 3 * re.begin();
    double *out = _Output + 3 * re.begin();

    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId, in += 3, out += 3) {
      norm = sqrt(vtkMath::Dot(in, in));
      if (norm) {
        avg_norm = norm, n = 1;
        _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPts);
        for (int i = 0; i < numAdjPts; ++i) {
          adj = _Input + 3 * adjPts[i];
          if (adj[0]*in[0] + adj[1]*in[1] + adj[2]*in[2] > .0) {
            avg_norm += sqrt(vtkMath::Dot(adj, adj)), ++n;
          }
        }
        avg_norm /= n;
        avg_norm /= norm;
        out[0] = in[0] * avg_norm;
        out[1] = in[1] * avg_norm;
        out[2] = in[2] * avg_norm;
      } else {
        out[0] = out[1] = out[2] = .0;
      }
    }
  }
};

// -----------------------------------------------------------------------------
template <class AvgFunc>
void AverageGradientVectors(double *dx, int n, const EdgeTable *edgeTable, int niter)
{
  AvgFunc avg;
  avg._Input     = dx;
  avg._Output    = Allocate<double>(n);
  avg._EdgeTable = edgeTable;
  blocked_range<vtkIdType> ptIds(0, edgeTable->NumberOfPoints());
  for (int iter = 0; iter < niter; ++iter) {
    parallel_for(ptIds, avg);
    swap(avg._Input, avg._Output);
  }
  if (avg._Output == dx) {
    delete[] avg._Input;
  } else {
    memcpy(dx, avg._Output, n * sizeof(double));
    delete[] avg._Output;
  }
}

// -----------------------------------------------------------------------------
/// Smooth gyral points along maximum curvature direction
vtkSmartPointer<vtkPolyData>
SmoothGyri(vtkSmartPointer<vtkPolyData> surface,
           int niter = 100, double lambda = .75, double mu = -.751)
{
  if (IsNaN(mu)) mu = lambda;

  SharedPtr<EdgeTable> edgeTable = NewShared<EdgeTable>(surface);

  int        numAdjPts;
  const int *adjPtIds;
  double     k2, k2_range[2], e2[3], p1[3], p2[3], p[3], e[3], d, w, wsum;
  double     current_lambda, alpha, beta;

  vtkSmartPointer<vtkDataArray> k2_array, e2_array;
  {
    SurfaceCurvature curvature;
    int curvature_type = SurfaceCurvature::Maximum;
    curvature_type    |= SurfaceCurvature::MaximumDirection;
    curvature.Input(surface);
    curvature.EdgeTable(edgeTable);
    curvature.CurvatureType(curvature_type);
    curvature.VtkCurvatures(false);
    curvature.TensorAveraging(3);
    curvature.Normalize(false);
    curvature.Run();
    vtkPointData *curvaturePD = curvature.Output()->GetPointData();
    k2_array = curvaturePD->GetArray(SurfaceCurvature::MAXIMUM);
    e2_array = curvaturePD->GetArray(SurfaceCurvature::MAXIMUM_DIRECTION);
  }

  vtkSmartPointer<vtkPoints> points;
  points = vtkSmartPointer<vtkPoints>::New();
  points->DeepCopy(surface->GetPoints());

  vtkSmartPointer<vtkPolyData> output;
  output.TakeReference(surface->NewInstance());
  output->ShallowCopy(surface);
  output->SetPoints(points);

  points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(surface->GetNumberOfPoints());

  k2_array->GetRange(k2_range);

  for (int iter = 1; iter <= niter; ++iter) {
    if (iter % 2 == 1) current_lambda = lambda;
    else               current_lambda = mu;
    for (vtkIdType ptId = 0; ptId < output->GetNumberOfPoints(); ++ptId) {
      output->GetPoint(ptId, p1);
      k2 = k2_array->GetComponent(ptId, 0);
      //if (k2 > .0) {
      {
        e2_array->GetTuple(ptId, e2);
        vtkMath::Normalize(e2);
        p[0] = p[1] = p[2] = wsum = .0;
        edgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPtIds);
        for (int i = 0; i < numAdjPts; ++i) {
          output->GetPoint(adjPtIds[i], p2);
          vtkMath::Subtract(p2, p1, e);
          d = vtkMath::Norm(e);
          vtkMath::MultiplyScalar(e, 1.0 / d);
          w = abs(vtkMath::Dot(e, e2));
          p[0] += w * p2[0];
          p[1] += w * p2[1];
          p[2] += w * p2[2];
          wsum += w;
        }
        if (wsum > .0) {
          alpha = current_lambda * (k2 / (k2 < .0 ? k2_range[0] : k2_range[1]));
          beta  = alpha / wsum;
          alpha = 1.0 - alpha;
          p1[0] = alpha * p1[0] + beta * p[0];
          p1[1] = alpha * p1[1] + beta * p[1];
          p1[2] = alpha * p1[2] + beta * p[2];
        }
      }
      points->SetPoint(ptId, p1);
    }
    output->GetPoints()->DeepCopy(points);
  }

  return output;
}

// -----------------------------------------------------------------------------
/// Compute an adaptive edge length for the local remeshing given an implicit target surface
vtkSmartPointer<vtkDataArray>
ComputeEdgeLengthRange(vtkPolyData *surface, double minl, double maxl, const DistanceImage *dmap, double offset = .0)
{
  static int iter = 0; ++iter;

#define CELL_EDGE_LENGTH_ARRAY 0

#if CELL_EDGE_LENGTH_ARRAY
  const int ntuples = static_cast<int>(surface->GetNumberOfCells());
#else
  const int ntuples = static_cast<int>(surface->GetNumberOfPoints());
#endif

  vtkSmartPointer<vtkDataArray> width = vtkSmartPointer<vtkFloatArray>::New();
  width->SetName("Width");
  width->SetNumberOfComponents(1);
  width->SetNumberOfTuples(ntuples);

  vtkSmartPointer<vtkDataArray> distance = vtkSmartPointer<vtkFloatArray>::New();
  distance->SetName("Distance");
  distance->SetNumberOfComponents(1);
  distance->SetNumberOfTuples(ntuples);

  vtkSmartPointer<vtkDataArray> min_edge_length = vtkSmartPointer<vtkFloatArray>::New();
  min_edge_length->SetName("MinEdgeLength");
  min_edge_length->SetNumberOfComponents(1);
  min_edge_length->SetNumberOfTuples(ntuples);

  vtkSmartPointer<vtkDataArray> max_edge_length = vtkSmartPointer<vtkFloatArray>::New();
  max_edge_length->SetName("MaxEdgeLength");
  max_edge_length->SetNumberOfComponents(1);
  max_edge_length->SetNumberOfTuples(ntuples);

  vtkSmartPointer<vtkDataArray> normals, cell_normals;
  normals      = surface->GetPointData()->GetNormals();
  cell_normals = surface->GetCellData ()->GetNormals();
  if (!normals || !cell_normals) {
    vtkNew<vtkPolyDataNormals> calc_normals;
    SetVTKInput(calc_normals, surface);
    calc_normals->ComputeCellNormalsOn();
    calc_normals->ComputePointNormalsOn();
    calc_normals->AutoOrientNormalsOn();
    calc_normals->SplittingOff();
    calc_normals->Update();
    normals      = calc_normals->GetOutput()->GetPointData()->GetNormals();
    cell_normals = calc_normals->GetOutput()->GetCellData ()->GetNormals();
  }

  const double ds   = min(min(dmap->GetXSize(), dmap->GetYSize()), dmap->GetZSize());
  const double maxw = 10.0;// * maxl;//10.0 * ds;
  const double maxd = maxl;
  const double dtol =  0.25  * ds;
  const double tol  =  0.01  * ds;
  const double minh =  0.001 * ds;

  // Parameters of logistic functions
  const double k  = 8.0;
  const double h1 = 0.4;
  const double h2 = 0.6;

  static DistanceImage width_image;
  if (iter == 1) {
//    MIRTK_START_TIMING();
//    width_image = Width<MinWidth>(maxw, *dmap, offset);
//    MIRTK_END_TIMING("computation of width of gaps");
//    width_image.Write("dmap_width.nii.gz");
    width_image.Read("dmap_width.nii.gz");
  }

  DistanceFunction width_map;
  width_map.Input(&width_image);
  width_map.Initialize();

  MIRTK_START_TIMING();
  double p[3], n[3];
  double mind, d, w, l1, l2, a, b, h;
  DistanceFunction distance_map;
  distance_map.Input(dmap);
  distance_map.Initialize();

#if CELL_EDGE_LENGTH_ARRAY
  int    subId;
  double pcoords[3], *weights = new double[surface->GetMaxCellSize()];
  double v0[3], v1[3], v2[3], e1[3], e2[3], e3[3], c[3], w1, w2, w3;
  vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
  for (vtkIdType cellId = 0; cellId < ntuples; ++cellId) {

    surface->GetCell(cellId, cell);
    cell_normals->GetTuple(cellId, n);
    subId = cell->GetParametricCenter(pcoords);
    cell->EvaluateLocation(subId, pcoords, c, weights);
    surface->GetPoint(cell->GetPointId(0), v0);
    surface->GetPoint(cell->GetPointId(1), v1);
    surface->GetPoint(cell->GetPointId(2), v2);
    vtkMath::Subtract(v1, v0, e1), vtkMath::Normalize(e1);
    vtkMath::Subtract(v2, v0, e2), vtkMath::Normalize(e2);
    vtkMath::Subtract(v2, v1, e3), vtkMath::Normalize(e3);

    mind = Evaluate(distance_map, c, offset);

    d = Distance(c, n, mind, minh, maxd, distance_map, offset, tol);
    p[0] = c[0] + .5 * (v0[0] - c[0]);
    p[1] = c[1] + .5 * (v0[1] - c[1]);
    p[2] = c[2] + .5 * (v0[2] - c[2]);
    d = max(d, Distance(p, n, minh, maxd, distance_map, offset, tol));
    p[0] = c[0] + .5 * (v1[0] - c[0]);
    p[1] = c[1] + .5 * (v1[1] - c[1]);
    p[2] = c[2] + .5 * (v1[2] - c[2]);
    d = max(d, Distance(p, n, minh, maxd, distance_map, offset, tol));
    p[0] = c[0] + .5 * (v2[0] - c[0]);
    p[1] = c[1] + .5 * (v2[1] - c[1]);
    p[2] = c[2] + .5 * (v2[2] - c[2]);
    d = max(d, Distance(p, n, minh, maxd, distance_map, offset, tol));

    if (mind <= .0 || d < dtol) {
      w = .0;
    } else {
      w1 = Width(c, e1, mind, minh, maxw, distance_map, offset, tol);
      w2 = Width(c, e2, mind, minh, maxw, distance_map, offset, tol);
      w3 = Width(c, e3, mind, minh, maxw, distance_map, offset, tol);
      w  = min(min(w1, w2), w3);
    }

    distance->SetComponent(cellId, 0, d);
    width   ->SetComponent(cellId, 0, w);

    // Desired edge length range
    h = d / (.95 * maxd);
    l1 = clamp(w / 6.0, minl, maxl);
    l2 = clamp(w / 4.0, minl, maxl);
    a = 1.0 / (1.0 + exp(- k * (h - h1))), b = 1.0 - a;
    min_edge_length->SetComponent(cellId, 0, clamp(a * l1 + b * minl, minl, maxl));
    a = 1.0 / (1.0 + exp(- k * (h - h2))), b = 1.0 - a;
    max_edge_length->SetComponent(cellId, 0, clamp(a * l2 + b * maxl, minl, maxl));
  }
  delete[] weights;
#else
  EdgeTable edgeTable(surface);

//  const int ndirs = 8;
//
//  const int    num_alpha = ndirs;
//  const double min_alpha = .0;
//  const double max_alpha = M_PI;
//  const double del_alpha = (max_alpha - min_alpha) / num_alpha;
//
//  const int    num_beta  = 5;
//  const double min_beta  = (num_beta == 1 ?  .0 : - M_PI * 20.0 / 180.0);
//  const double max_beta  = (num_beta == 1 ?  .0 :   M_PI * 20.0 / 180.0);
//  const double del_beta  = (num_beta == 1 ? 1.0 : (max_beta - min_beta) / (num_beta - 1));
//
//  double cos_beta, sin_beta, cos_alpha, sin_alpha, e[3], e1[3], e2[3];
//
//  PointSet dirs(num_alpha * num_beta);

  int        numAdjPts;
  const int *adjPtIds;
  double     p2[3], e1[3], e2[3];

  for (vtkIdType ptId = 0; ptId < ntuples; ++ptId) {
    surface->GetPoint(ptId, p);
    normals->GetTuple(ptId, n);
    ComputeTangents(n, e1, e2);
    edgeTable.GetAdjacentPoints(ptId, numAdjPts, adjPtIds);

    mind = Evaluate(distance_map, p, offset);

    d = Distance(p, n, mind, minh, maxw, distance_map, offset, tol);
    for (int i = 0; i < numAdjPts; ++i ) {
      surface->GetPoint(adjPtIds[i], p2);
      d = max(d, Distance(p2, n, minh, maxw, distance_map, offset, tol));
    }

    if (mind <= .0 || d < dtol) {
      w = .0;
    } else {
      double x = p[0], y = p[1], z = p[2];
      width_map.WorldToImage(x, y, z);
      w = width_map.Evaluate(x, y, z);
//      int i = 0;
//      ComputeTangents(n, e1, e2);
//      for (double alpha = min_alpha; alpha < max_alpha; alpha += del_alpha) {
//        cos_alpha = cos(alpha);
//        sin_alpha = sin(alpha);
//        e[0] = e1[0] * cos_alpha + e2[0] * sin_alpha;
//        e[1] = e1[1] * cos_alpha + e2[1] * sin_alpha;
//        e[2] = e1[2] * cos_alpha + e2[2] * sin_alpha;
//        for (double beta = min_beta; beta <= max_beta; beta += del_beta) {
//          cos_beta = cos(beta);
//          sin_beta = sin(beta);
//          dirs(i)._x = e[0] * cos_beta + n[0] * sin_beta;
//          dirs(i)._y = e[1] * cos_beta + n[1] * sin_beta;
//          dirs(i)._z = e[2] * cos_beta + n[2] * sin_beta;
//          ++i;
//        }
//      }
//      w = Evaluate<MinWidth>(p, dirs, mind, minh, maxw, distance_map, offset, tol);
      //w = EvaluateInTangentPlane<MinWidth>(p, n, mind, minh, maxw, distance_map, offset, tol, ndirs, 20.0);
    }

    distance->SetComponent(ptId, 0, d);
    width   ->SetComponent(ptId, 0, w);

    // Desired edge length range
    h = d / (.95 * maxd);
    l1 = clamp(w / 6.0, minl, maxl);
    l2 = clamp(w / 4.0, minl, maxl);
    a = 1.0 / (1.0 + exp(- k * (h - h1))), b = 1.0 - a;
    min_edge_length->SetComponent(ptId, 0, clamp(a * l1 + b * minl, minl, maxl));
    a = 1.0 / (1.0 + exp(- k * (h - h2))), b = 1.0 - a;
    max_edge_length->SetComponent(ptId, 0, clamp(a * l2 + b * maxl, minl, maxl));
  }
#endif
  MIRTK_DEBUG_TIMING(5, "computation of adaptive edge length");

  vtkDataSetAttributes *attrs;
#if CELL_EDGE_LENGTH_ARRAY
  attrs = surface->GetCellData();
#else
  attrs = surface->GetPointData();
#endif
  attrs->RemoveArray(distance->GetName());
  attrs->RemoveArray(width->GetName());
  attrs->RemoveArray(min_edge_length->GetName());
  attrs->RemoveArray(max_edge_length->GetName());

  attrs->AddArray(distance);
  attrs->AddArray(width);
  attrs->AddArray(min_edge_length);
  attrs->AddArray(max_edge_length);

  char fname[64];
  #ifdef WINDOWS
    sprintf_s(fname, 64, "surface_distance_%03d.vtp", iter);
  #else
    sprintf(fname, "surface_distance_%03d.vtp", iter);
  #endif
  WritePolyData(fname, surface);

  attrs->RemoveArray(min_edge_length->GetName());
  attrs->RemoveArray(max_edge_length->GetName());

  return NULL;
}


} // namespace DeformableSurfaceModelUtils
using namespace DeformableSurfaceModelUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
DeformableSurfaceModel::DeformableSurfaceModel()
:
  _Image(NULL),
  _ImplicitSurface(NULL),
  _Transformation(NULL),
  _NumberOfTerms(0),
  _NeighborhoodRadius(2),
  _GradientAveraging(0),
  _GradientWeighting(MeshSmoothing::Default),
  _AverageSignedGradients(false),
  _AverageGradientMagnitude(false),
  _MinEdgeLength(-1.0),
  _MaxEdgeLength(-1.0),
  _MinFeatureAngle(180.0),
  _MaxFeatureAngle(180.0),
  _AllowTriangleInversion(true),
  _RemeshInterval(0),
  _RemeshCounter(0),
  _RemeshAdaptively(false),
  _LowPassInterval(0),
  _LowPassIterations(100),
  _LowPassBand(.75),
  _MaxInputDistance(inf),
  _HardNonSelfIntersection(false),
  _MinFrontfaceDistance(.0),
  _MinBackfaceDistance(.0),
  _MaxCollisionAngle(45.0),
  _FastCollisionTest(false),
  _FixPassivePoints(true),
  _AllowExpansion(true),
  _AllowContraction(true),
  _IsSurfaceMesh(false),
  _MinimizeExtrinsicEnergy(false),
  _LowPassCounter(0)
{
}

// -----------------------------------------------------------------------------
DeformableSurfaceModel::~DeformableSurfaceModel()
{
  Clear();
}

// =============================================================================
// Energy terms
// =============================================================================

// -----------------------------------------------------------------------------
void DeformableSurfaceModel::Initialize()
{
  // Check input
  if (_Input == NULL || _Input->GetNumberOfPoints() == 0) {
    cerr << "DeformableSurfaceModel::Initialize: Missing initial surface mesh" << endl;
    exit(1);
  }
  if (_NumberOfTerms == 0) {
    cerr << "DeformableSurfaceModel::Initialize: No internal and/or external forces added" << endl;
    exit(1);
  }
  if (!_Transformation) {
    for (size_t i = 0; i < _Constraint.size(); ++i) _Constraint[i]->Weight(.0);
  }

  // Whether deformable model is a surface mesh
  _IsSurfaceMesh = mirtk::IsSurfaceMesh(_Input);

  // Initialize input locator if needed
  if (_MaxInputDistance <= 0.) _MaxInputDistance = inf;
  if (_IsSurfaceMesh && !IsInf(_MaxInputDistance)) {
    _InputCellLocator = vtkSmartPointer<vtkCellLocator>::New();
    _InputCellLocator->SetDataSet(_Input);
    _InputCellLocator->BuildLocator();
  }

  // Gradient smoothing
  if (_Transformation && (_GradientAveraging > 0 || _LowPassInterval > 0)) {
    cerr << "DeformableSurfaceMode::Initialize: "
                 " Gradient averaging and/or low-pass filtering only allowed for"
                 " non-parametric deformable surface models" << endl;
    exit(1);
  }
  if (_LowPassInterval > 0 && _LowPassIterations > 0 && !_IsSurfaceMesh) {
    cerr << "DeformableSurfaceMode::Initialize: "
                 " Low-pass filtering only allowed for non-parametric deformable"
                 " surface models, but input is not a surface mesh" << endl;
    exit(1);
  }
  _LowPassCounter = 0;

  // Local adaptive remeshing settings
  if (_IsSurfaceMesh && IsTriangularMesh(_Input) && _RemeshInterval > 0) {
    if (_MinEdgeLength < .0 || _MaxEdgeLength <= .0) {
      double mean = AverageEdgeLength(_Input);
      if (_MinEdgeLength <  .0) _MinEdgeLength =  .5 * mean;
      if (_MaxEdgeLength <= .0) _MaxEdgeLength = 2.0 * mean;
    }
  } else {
    _RemeshInterval = 0;
  }
  _RemeshCounter = 0;

  // Initialize output surface mesh
  _PointSet.InputPointSet(_Input);
  _PointSet.Transformation(_Transformation);
  _PointSet.SelfUpdate(false);
  _PointSet.NeighborhoodRadius(_NeighborhoodRadius);
  _PointSet.Initialize(_Transformation == NULL);
  this->Changed(true);

  // Initialize energy terms
  for (size_t i = 0; i < _ExternalForce.size(); ++i) {
    _ExternalForce[i]->PointSet(&_PointSet);
    if (IsImplicitSurfaceForce(_ExternalForce[i])) {
      _ExternalForce[i]->Image(_ImplicitSurface);
    } else {
      _ExternalForce[i]->Image(_Image);
    }
  }
  for (size_t i = 0; i < _InternalForce.size(); ++i) {
    _InternalForce[i]->PointSet(&_PointSet);
  }
  for (int i = 0; i < _NumberOfTerms; ++i) {
    EnergyTerm *term = Term(i);
    if (term->Weight() != .0) {
      term->Transformation(_Transformation);
      term->Initialize();
    }
  }
}

// -----------------------------------------------------------------------------
void DeformableSurfaceModel::Clear()
{
  for (size_t i = 0; i < _ExternalForce.size(); ++i) {
    if (_ExternalForceOwner[i]) Delete(_ExternalForce[i]);
  }
  for (size_t i = 0; i < _InternalForce.size(); ++i) {
    if (_InternalForceOwner[i]) Delete(_InternalForce[i]);
  }
  for (size_t i = 0; i < _Constraint.size(); ++i) {
    if (_ConstraintOwner[i]) Delete(_Constraint[i]);
  }
  _ExternalForce.clear();
  _ExternalForceOwner.clear();
  _InternalForce.clear();
  _InternalForceOwner.clear();
  _Constraint.clear();
  _ConstraintOwner.clear();
  _NumberOfTerms = 0;
}

// -----------------------------------------------------------------------------
bool DeformableSurfaceModel::Empty() const
{
  return _NumberOfTerms == 0;
}

// -----------------------------------------------------------------------------
int DeformableSurfaceModel::NumberOfForces() const
{
  return NumberOfInternalForces() + NumberOfExternalForces();
}

// -----------------------------------------------------------------------------
int DeformableSurfaceModel::NumberOfInternalForces() const
{
  return static_cast<int>(_InternalForce.size());
}

// -----------------------------------------------------------------------------
int DeformableSurfaceModel::NumberOfExternalForces() const
{
  return static_cast<int>(_ExternalForce.size());
}

// -----------------------------------------------------------------------------
void DeformableSurfaceModel::Add(class ExternalForce *term, bool ownit)
{
  _ExternalForce.push_back(term);
  _ExternalForceOwner.push_back(ownit);
  ++_NumberOfTerms;
}

// -----------------------------------------------------------------------------
void DeformableSurfaceModel::Sub(class ExternalForce *term)
{
  Array<class ExternalForce *>::iterator it = _ExternalForce.begin();
  while (it != _ExternalForce.end()) {
    if (*it == term) {
      _ExternalForce.erase(it);
      _ExternalForceOwner.erase(_ExternalForceOwner.begin() + distance(_ExternalForce.begin(), it));
      --_NumberOfTerms;
      break;
    }
    ++it;
  }
}

// -----------------------------------------------------------------------------
void DeformableSurfaceModel::Add(class InternalForce *term, bool ownit)
{
  _InternalForce.push_back(term);
  _InternalForceOwner.push_back(ownit);
  ++_NumberOfTerms;
}

// -----------------------------------------------------------------------------
void DeformableSurfaceModel::Sub(class InternalForce *term)
{
  for (auto it = _InternalForce.begin(); it != _InternalForce.end(); ++it) {
    if (*it == term) {
      _InternalForce.erase(it);
      _InternalForceOwner.erase(_InternalForceOwner.begin() + distance(_InternalForce.begin(), it));
      --_NumberOfTerms;
      break;
    }
  }
}

// -----------------------------------------------------------------------------
void DeformableSurfaceModel::Add(TransformationConstraint *term, bool ownit)
{
  _Constraint.push_back(term);
  _ConstraintOwner.push_back(ownit);
  ++_NumberOfTerms;
}

// -----------------------------------------------------------------------------
void DeformableSurfaceModel::Sub(TransformationConstraint *term)
{
  for (auto it = _Constraint.begin(); it != _Constraint.end(); ++it) {
    if (*it == term) {
      _Constraint.erase(it);
      _ConstraintOwner.erase(_ConstraintOwner.begin() + distance(_Constraint.begin(), it));
      --_NumberOfTerms;
      break;
    }
  }
}

// -----------------------------------------------------------------------------
EnergyTerm *DeformableSurfaceModel::Term(int i)
{
  if (i < NumberOfExternalForces()) return _ExternalForce[i];
  i -= NumberOfExternalForces();
  if (i < NumberOfInternalForces()) return _InternalForce[i];
  i -= NumberOfInternalForces();
  return _Constraint[i];
}

// -----------------------------------------------------------------------------
const EnergyTerm *DeformableSurfaceModel::Term(int i) const
{
  if (i < NumberOfExternalForces()) return _ExternalForce[i];
  i -= NumberOfExternalForces();
  if (i < NumberOfInternalForces()) return _InternalForce[i];
  i -= NumberOfInternalForces();
  return _Constraint[i];
}

// -----------------------------------------------------------------------------
ExternalForce *DeformableSurfaceModel::ExternalForce(int i)
{
  return _ExternalForce[i];
}

// -----------------------------------------------------------------------------
const ExternalForce *DeformableSurfaceModel::ExternalForce(int i) const
{
  return _ExternalForce[i];
}

// -----------------------------------------------------------------------------
InternalForce *DeformableSurfaceModel::InternalForce(int i)
{
  return _InternalForce[i];
}

// -----------------------------------------------------------------------------
const InternalForce *DeformableSurfaceModel::InternalForce(int i) const
{
  return _InternalForce[i];
}

// -----------------------------------------------------------------------------
bool DeformableSurfaceModel::IsExternalForce(const EnergyTerm *term)
{
  return dynamic_cast<const class ExternalForce *>(term) != NULL;
}

// -----------------------------------------------------------------------------
bool DeformableSurfaceModel::IsImplicitSurfaceForce(const EnergyTerm *term)
{
  return dynamic_cast<const ImplicitSurfaceForce *>(term) != NULL;
}

// -----------------------------------------------------------------------------
bool DeformableSurfaceModel::IsInternalForce(const EnergyTerm *term)
{
  return dynamic_cast<const class InternalForce *>(term) != NULL;
}

// -----------------------------------------------------------------------------
bool DeformableSurfaceModel::IsExternalForce(int i) const
{
  return IsExternalForce(Term(i));
}

// -----------------------------------------------------------------------------
bool DeformableSurfaceModel::IsImplicitSurfaceForce(int i) const
{
  return IsImplicitSurfaceForce(Term(i));
}

// -----------------------------------------------------------------------------
bool DeformableSurfaceModel::IsInternalForce(int i) const
{
  return IsInternalForce(Term(i));
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool DeformableSurfaceModel::Set(const char *name, const char *value)
{
  if (strcmp(name, "No. of gradient averaging iterations") == 0) {
    return FromString(value, _GradientAveraging);
  }
  if (strcmp(name, "Average gradient vectors with same sign") == 0) {
    return FromString(value, _AverageSignedGradients);
  }
  if (strcmp(name, "Average magnitude of gradient vectors") == 0) {
    return FromString(value, _AverageGradientMagnitude);
  }
  if (strcmp(name, "Minimum edge length") == 0) {
    return FromString(value, _MinEdgeLength);
  }
  if (strcmp(name, "Maximum edge length") == 0) {
    return FromString(value, _MaxEdgeLength);
  }
  if (strcmp(name, "Minimum feature angle") == 0 || strcmp(name, "Minimum edge angle") == 0) {
    return FromString(value, _MinFeatureAngle);
  }
  if (strcmp(name, "Maximum feature angle") == 0 || strcmp(name, "Maximum edge angle") == 0) {
    return FromString(value, _MaxFeatureAngle);
  }
  if (strcmp(name, "Remesh interval") == 0) {
    return FromString(value, _RemeshInterval);
  }
  if (strcmp(name, "Adatpive remeshing") == 0 || strcmp(name, "Remesh adaptively") == 0) {
    return FromString(value, _RemeshAdaptively);
  }
  if (strcmp(name, "Maximum distance from input surface") == 0) {
    return FromString(value, _MaxInputDistance);
  }
  if (strcmp(name, "Hard non-self-intersection constraint") == 0) {
    return FromString(value, _HardNonSelfIntersection);
  }
  if (strcmp(name, "Minimum frontface distance") == 0) {
    return FromString(value, _MinFrontfaceDistance);
  }
  if (strcmp(name, "Minimum backface distance") == 0) {
    return FromString(value, _MinBackfaceDistance);
  }
  if (strcmp(name, "Maximum collision angle") == 0) {
    return FromString(value, _MaxCollisionAngle);
  }
  if (strcmp(name, "Fast collision test") == 0) {
    return FromString(value, _FastCollisionTest);
  }
  if (strcmp(name, "Allow triangle inversion") == 0) {
    return FromString(value, _AllowTriangleInversion);
  }
  if (strcmp(name, "Allow surface expansion") == 0) {
    return FromString(value, _AllowExpansion);
  }
  if (strcmp(name, "Allow surface contraction") == 0) {
    return FromString(value, _AllowContraction);
  }

  bool known = false;
  for (int i = 0; i < _NumberOfTerms; ++i) {
    known = Term(i)->Set(name, value) || known;
  }
  return known;
}

// -----------------------------------------------------------------------------
ParameterList DeformableSurfaceModel::Parameter() const
{
  ParameterList params;
  for (int i = 0; i < _NumberOfTerms; ++i) {
    Insert(params, Term(i)->Parameter());
  }
  Insert(params, "No. of gradient averaging iterations", _GradientAveraging);
  Insert(params, "Average gradient vectors with same sign", _AverageSignedGradients);
  Insert(params, "Average magnitude of gradient vectors", _AverageGradientMagnitude);
  Insert(params, "Minimum edge length", _MinEdgeLength);
  Insert(params, "Maximum edge length", _MaxEdgeLength);
  Insert(params, "Minimum feature angle", _MinFeatureAngle);
  Insert(params, "Maximum feature angle", _MaxFeatureAngle);
  Insert(params, "Remesh interval", _RemeshInterval);
  Insert(params, "Adaptive remeshing", _RemeshAdaptively);
  Insert(params, "Maximum distance from input surface", _MaxInputDistance);
  Insert(params, "Hard non-self-intersection constraint", _HardNonSelfIntersection);
  Insert(params, "Minimum frontface distance", _MinFrontfaceDistance);
  Insert(params, "Minimum backface distance", _MinBackfaceDistance);
  Insert(params, "Maximum collision angle", _MaxCollisionAngle);
  Insert(params, "Fast collision test", _FastCollisionTest);
  Insert(params, "Allow triangle inversion", _AllowTriangleInversion);
  Insert(params, "Allow surface expansion", _AllowExpansion);
  Insert(params, "Allow surface contraction", _AllowContraction);
  return params;
}

// =============================================================================
// Degrees of freedom
// =============================================================================

// -----------------------------------------------------------------------------
int DeformableSurfaceModel::NumberOfDOFs() const
{
  return _Transformation ? _Transformation->NumberOfDOFs() : 3 * _PointSet.NumberOfPoints();
}

// -----------------------------------------------------------------------------
void DeformableSurfaceModel::Put(const double *x)
{
  if (_Transformation) {
    _Transformation->Put(x);
  } else {
    // Set vertex positions of deformable surface mesh
    vtkPoints *points = _PointSet.Points();
    for (int i = 0; i < _PointSet.NumberOfPoints(); ++i, x += 3) {
      points->SetPoint(i, x);
    }
    _PointSet.PointsChanged();
  }
  // Mark deformable surface model as changed
  this->Changed(true);
  for (int i = 0; i < _NumberOfTerms; ++i) Term(i)->ResetValue();
}

// -----------------------------------------------------------------------------
void DeformableSurfaceModel::Get(double *x) const
{
  if (_Transformation) {
    _Transformation->Get(x);
  } else {
    for (int i = 0; i < _PointSet.NumberOfPoints(); ++i, x += 3) {
      _PointSet.GetPoint(i, x);
    }
  }
}

// -----------------------------------------------------------------------------
double DeformableSurfaceModel::Get(int dof) const
{
  if (_Transformation) {
    return _Transformation->Get(dof);
  } else {
    double p[3];
    _PointSet.GetPoint(dof / 3, p);
    return p[dof % 3];
  }
}

// -----------------------------------------------------------------------------
double DeformableSurfaceModel::Step(double *dx)
{
  double delta;
  if (_Transformation) {

    delta = _Transformation->Update(dx);

  } else {

    // Perform low-pass filtering
    //
    // The translation and scale "fix" is due to a bug in vtkWindowedSincPolyDataFilter:
    // http://vtk.1045678.n5.nabble.com/Bug-in-vtkWindowedSincPolyDataFilter-td1234055.html
    if (_IsSurfaceMesh && _LowPassInterval > 0 && _LowPassIterations > 0) {
      ++_LowPassCounter;
      if (_LowPassCounter >= _LowPassInterval) {
        _LowPassCounter = 0;

        MIRTK_START_TIMING();
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        points->SetNumberOfPoints(_PointSet.NumberOfPoints());
        vtkSmartPointer<vtkPolyData> surface;
        surface = vtkSmartPointer<vtkPolyData>::New();
        surface->ShallowCopy(_PointSet.Surface());
        surface->SetPoints(points);
        MovePoints::Run(_PointSet.Points(), dx, surface->GetPoints());

        double c1[3], s1[3];
        GetCentroid(surface, c1);
        GetScale(surface, c1, s1);

        vtkNew<vtkWindowedSincPolyDataFilter> filter;
        filter->SetPassBand(_LowPassBand);
        filter->SetNumberOfIterations(_LowPassIterations);
        filter->NormalizeCoordinatesOn();
        filter->FeatureEdgeSmoothingOff();
        SetVTKInput(filter, surface);
        filter->Update();

        double c2[3], s2[3];
        GetCentroid(filter->GetOutput(), c2);
        GetScale(filter->GetOutput(), c2, s2);

        double p1[3], p2[3], *dp = dx;
        for (vtkIdType ptId = 0; ptId < surface->GetNumberOfPoints(); ++ptId, dp += 3) {
          surface->GetPoint(ptId, p1);
          filter->GetOutput()->GetPoint(ptId, p2);
          dp[0] += (c1[0] + s1[0] * (p2[0] - c2[0]) / s2[0]) - p1[0];
          dp[1] += (c1[1] + s1[1] * (p2[1] - c2[1]) / s2[1]) - p1[1];
          dp[2] += (c1[2] + s1[2] * (p2[2] - c2[2]) / s2[2]) - p1[2];
        }
        MIRTK_DEBUG_TIMING(3, "low-pass filtering");
      }
    }
    // Enforce hard constraints
    this->EnforceHardConstraints(dx);
    // Determine maximum vertex displacement
    MaxVertexDisplacement max(dx);
    parallel_reduce(blocked_range<int>(0, _PointSet.NumberOfPoints()), max);
    delta = max.Norm();
    // Update points of output surface
    if (delta != .0) {
      MovePoints::Run(_PointSet.Points(), dx);
      _PointSet.PointsChanged();
    }
  }
  // Mark deformable surface model as changed
  if (delta != .0) {
    this->Changed(true);
    for (int i = 0; i < _NumberOfTerms; ++i) Term(i)->ResetValue();
  }
  // Return maximum vertex displacement
  return delta;
}

// -----------------------------------------------------------------------------
void DeformableSurfaceModel::Update(bool gradient)
{
  if (this->Changed() || gradient) {
    MIRTK_START_TIMING();
    // Update deformed point set
    if (_Transformation) {
      _PointSet.Update(true);
    }
    // 1. Update external forces
    for (auto force : _ExternalForce) {
      if (force->Weight() != 0.) {
        force->Update(gradient);
        force->ResetValue(); // in case energy term does not do this
      }
    }
    // 2. Update internal forces
    for (auto force : _InternalForce) {
      if (force->Weight() != 0.) {
        force->Update(gradient);
        force->ResetValue(); // in case energy term does not do this
      }
    }
    // 3. Update transformation constraints
    for (auto constraint : _Constraint) {
      if (constraint->Weight() != 0.) {
        constraint->Update(gradient);
        constraint->ResetValue(); // in case energy term does not do this
      }
    }
    // Mark deformable surface model as up-to-date
    this->Changed(false);
    MIRTK_DEBUG_TIMING(3, "update of energy function");
  }
}

// -----------------------------------------------------------------------------
bool DeformableSurfaceModel::Upgrade()
{
  return false;
}

// -----------------------------------------------------------------------------
bool DeformableSurfaceModel::Remesh()
{
  // Currently only remeshing of a triangulated surface mesh is supported
  if (_RemeshInterval == 0) return false;

  // Increment counter counting the number of iterations since the last remeshing
  ++_RemeshCounter;
  if (_RemeshCounter < _RemeshInterval) return false;
  _RemeshCounter = 0;

  MIRTK_START_TIMING();

  // Shallow copy of surface mesh
  vtkSmartPointer<vtkPolyData> input;

  if (_Transformation) {
    input = _PointSet.InputSurface();
  } else {
    input = _PointSet.Surface();
    if (!input->GetPointData()->HasArray("InitialPoints")) {
      vtkSmartPointer<vtkDataArray> initial_points;
      initial_points = vtkSmartPointer<vtkFloatArray>::New();
      initial_points->SetName("InitialPoints");
      initial_points->SetNumberOfComponents(3);
      initial_points->SetNumberOfTuples(input->GetNumberOfPoints());
      vtkPoints *points = input->GetPoints();
      for (vtkIdType ptId = 0; ptId < input->GetNumberOfPoints(); ++ptId) {
        initial_points->SetTuple(ptId, points->GetPoint(ptId));
      }
      input->GetPointData()->AddArray(initial_points);
    }
  }

  // Compute local edge length intervals
  SurfaceRemeshing remesher;
  if (_RemeshAdaptively) {
    vtkSmartPointer<vtkDataArray> adaptive_edge_length;
    if (false && _ImplicitSurface) {
      // FIXME: Experimental code that is not ready for use
      adaptive_edge_length = ComputeEdgeLengthRange(_PointSet.Surface(), _MinEdgeLength, _MaxEdgeLength, _ImplicitSurface);
    } else {
      SurfaceCurvature curv;
      curv.Input(_PointSet.Surface());
      curv.CurvatureType(SurfaceCurvature::Curvedness);
      curv.Run();
      adaptive_edge_length = curv.GetCurvedness();
      _PointSet.Surface()->GetPointData()->AddArray(adaptive_edge_length);
    }
    if (adaptive_edge_length) {
      remesher.AdaptiveEdgeLengthArray(adaptive_edge_length);
    } else {
      remesher.MinCellEdgeLengthArray(input->GetCellData()->GetArray("MinEdgeLength"));
      remesher.MaxCellEdgeLengthArray(input->GetCellData()->GetArray("MaxEdgeLength"));
    }
  }

  // Remesh surface
  remesher.Input(input);
  remesher.PointMask(_PointSet.InitialSurfaceStatus());
  remesher.MeltingOrder(SurfaceRemeshing::AREA);
  remesher.MeltNodesOn();
  remesher.MeltTrianglesOff();
  remesher.InvertTrianglesSharingOneLongEdge(_AllowTriangleInversion);
  remesher.InvertTrianglesToIncreaseMinHeight(_AllowTriangleInversion);
  remesher.MinEdgeLength(_MinEdgeLength);
  remesher.MaxEdgeLength(_MaxEdgeLength);
  remesher.MinFeatureAngle(_MinFeatureAngle);
  remesher.MaxFeatureAngle(_MaxFeatureAngle);
  remesher.Transformation(_Transformation);
  remesher.Run();

  vtkSmartPointer<vtkPolyData> output = remesher.Output();

  if (output != input) {
    // Update deformable surface mesh
    _PointSet.InputPointSet(output);
    if (_Transformation) {
      _PointSet.Initialize(false);
    } else {
      vtkSmartPointer<vtkPoints> points;
      if (output->GetPointData()->HasArray("InitialPoints")) {
        points = output->GetPoints();
        // Reset points of input surface mesh
        vtkSmartPointer<vtkPoints> initial_points;
        initial_points = vtkSmartPointer<vtkPoints>::New();
        initial_points->SetNumberOfPoints(points->GetNumberOfPoints());
        vtkDataArray *initial_pos = output->GetPointData()->GetArray("InitialPoints");
        for (vtkIdType ptId = 0; ptId < points->GetNumberOfPoints(); ++ptId) {
          initial_points->SetPoint(ptId, initial_pos->GetTuple(ptId));
        }
        output->SetPoints(initial_points);
      }
      // Initialize surface mesh and set new output points
      _PointSet.Initialize(false);
      if (points) _PointSet.Points()->DeepCopy(points);
    }

    // Reinitialize internal and external force terms
    for (size_t i = 0; i < _ExternalForce.size(); ++i) {
      if (_ExternalForce[i]->Weight() != .0) {
        _ExternalForce[i]->Reinitialize();
      }
    }
    for (size_t i = 0; i < _InternalForce.size(); ++i) {
      if (_InternalForce[i]->Weight() != .0) {
        _InternalForce[i]->Reinitialize();
      }
    }

    // Mark deformable surface model as modified
    this->Changed(true);
  }

  MIRTK_DEBUG_TIMING(3, "local adaptive remeshing");
  return (input != output);
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double DeformableSurfaceModel::RawValue(int i)
{
  EnergyTerm *term = Term(i);
  return (term->Weight() != .0 ? term->RawValue() : .0);
}

// -----------------------------------------------------------------------------
double DeformableSurfaceModel::InitialValue()
{
  MIRTK_START_TIMING();

  double value, sum = .0;
  for (int i = 0; i < _NumberOfTerms; ++i) {
    EnergyTerm *term = Term(i);
    value = (term->Weight() != .0 ? term->Value() : .0);
    if (IsNaN(value)) {
      string name = term->Name();
      if (name.empty()) name = ToString(i + 1);
      cerr << "DeformableSurfaceModel::InitialValue: Value of term " << name << " is NaN!" << endl;
      exit(1);
    }
    if (_MinimizeExtrinsicEnergy && i >= NumberOfExternalForces()) continue;
    sum += value;
  }

  MIRTK_DEBUG_TIMING(3, "initial evaluation of energy function");
  return sum;
}

// -----------------------------------------------------------------------------
double DeformableSurfaceModel::InitialValue(int i)
{
  return Term(i)->InitialValue();
}

// -----------------------------------------------------------------------------
double DeformableSurfaceModel::Value()
{
  MIRTK_START_TIMING();

  double value, sum = .0;
  for (int i = 0; i < _NumberOfTerms; ++i) {
    EnergyTerm *term = Term(i);
    value = (term->Weight() != .0 ? term->Value() : .0);
    if (IsNaN(value)) {
      string name = term->Name();
      if (name.empty()) name = ToString(i + 1);
      cerr << "DeformableSurfaceModel::Value: Value of term " << name << " is NaN!" << endl;
      exit(1);
    }
    if (_MinimizeExtrinsicEnergy && i >= NumberOfExternalForces()) continue;
    sum += value;
  }

  MIRTK_DEBUG_TIMING(3, "evaluation of energy function");
  return sum;
}

// -----------------------------------------------------------------------------
double DeformableSurfaceModel::Value(int i)
{
  return Term(i)->Value();
}

// -----------------------------------------------------------------------------
void DeformableSurfaceModel::Gradient(double *gradient, double step, bool *sgn_chg)
{
  MIRTK_START_TIMING();

  // Use default step length if none specified
  if (step <= .0) step = _StepLength;

  // Initialize output variables
  const int ndofs = this->NumberOfDOFs();
  memset(gradient, 0, ndofs * sizeof(double));
  if (sgn_chg) {
    for (int dof = 0; dof < ndofs; ++dof) {
      sgn_chg[dof] = true;
    }
  }

  // Sum (weighted) internal and external forces
  for (int i = 0; i < _NumberOfTerms; ++i) {
    EnergyTerm *term = Term(i);
    if (term->Weight() != .0) {
      term->Gradient(gradient, step);
    }
  }

  // Smooth gradient
  this->SmoothGradient(gradient);

  MIRTK_DEBUG_TIMING(3, "evaluation of energy gradient");
}

// -----------------------------------------------------------------------------
double DeformableSurfaceModel::GradientNorm(const double *dx) const
{
  if (_Transformation) return _Transformation->DOFGradientNorm(dx);
  MaxVertexDisplacement max(dx);
  parallel_reduce(blocked_range<int>(0, _PointSet.NumberOfPoints()), max);
  if (max.FoundNaN()) {
    cerr << "DeformableSurfaceModel::GradientNorm: Gradient vector contains NaNs" << endl;
    exit(1);
  }
  return max.Norm();
}

// -----------------------------------------------------------------------------
void DeformableSurfaceModel::GradientStep(const double *dx, double &min, double &max) const
{
  for (int i = 0; i < _NumberOfTerms; ++i) {
    const EnergyTerm *term = Term(i);
    if (term->Weight() != .0) {
      term->GradientStep(dx, min, max);
    }
  }
}

// -----------------------------------------------------------------------------
double DeformableSurfaceModel::Evaluate(double *dx, double step, bool *sgn_chg)
{
  // Update energy function
  if (this->Changed()) this->Update(dx != NULL);

  // Evaluate gradient
  if (dx) this->Gradient(dx, step, sgn_chg);

  // Evaluate energy
  return this->Value();
}

// -----------------------------------------------------------------------------
void DeformableSurfaceModel::SmoothGradient(double *dx) const
{
  // Only done when DoFs are the node positions themselves
  if (_Transformation) return;

  // Smooth vertex displacements such that adjacent nodes move coherently.
  // Can also be viewed as an averaging of the gradient vectors in a local
  // neighborhood. With decreasing smoothing iterations, a multi-resolution
  // optimization of the deformable surface model can be mimicked.
  if (_GradientAveraging > 0) {
    MIRTK_START_TIMING();
    const int ndofs = this->NumberOfDOFs();

    if (_IsSurfaceMesh && _GradientWeighting != MeshSmoothing::Combinatorial) {

      vtkSmartPointer<vtkPolyData> surface;
      surface = vtkSmartPointer<vtkPolyData>::New();
      surface->ShallowCopy(_PointSet.Surface());
      surface->GetCellData()->Initialize();
      surface->GetPointData()->Initialize();

      vtkSmartPointer<vtkDataArray> gradient;
      gradient = vtkSmartPointer<vtkDoubleArray>::New();
      gradient->SetName("Gradient");
      gradient->SetNumberOfComponents(3);
      gradient->SetNumberOfTuples(ndofs / 3);
      memcpy(gradient->GetVoidPointer(0), dx, ndofs * sizeof(double));
      surface->GetPointData()->AddArray(gradient);

      MeshSmoothing smoother;
      smoother.Input(surface);
      smoother.Lambda(1.0);
      smoother.SmoothPointsOff();
      smoother.AdjacentValuesOnlyOff();
      smoother.SignedSmoothing(_AverageSignedGradients);
      smoother.SmoothMagnitude(_AverageGradientMagnitude);
      smoother.SmoothArray("Gradient");
      if (_GradientWeighting == MeshSmoothing::Default) {
        smoother.Weighting(MeshSmoothing::InverseDistance);
      } else {
        smoother.Weighting(_GradientWeighting);
      }
      smoother.NumberOfIterations(_GradientAveraging);
      smoother.Run();

      gradient = smoother.Output()->GetPointData()->GetArray("Gradient");
      memcpy(dx, gradient->GetVoidPointer(0), ndofs * sizeof(double));

    } else {

      if (_GradientWeighting != MeshSmoothing::Combinatorial &&
          _GradientWeighting != MeshSmoothing::Default) {
        cerr << this->NameOfType() << "::SmoothGradient: Only combinatorial weighting available"
                                      " for non-polygonal point sets" << endl;
        exit(1);
      }

      const EdgeTable *edgeTable = _PointSet.Edges();
      if (_AverageSignedGradients) {
        if (_AverageGradientMagnitude) {
          typedef struct AverageSignedGradientMagnitude AvgOp;
          AverageGradientVectors<AvgOp>(dx, ndofs, edgeTable, _GradientAveraging);
        } else {
          typedef struct AverageSignedGradient AvgOp;
          AverageGradientVectors<AvgOp>(dx, ndofs, edgeTable, _GradientAveraging);
        }
      } else {
        if (_AverageGradientMagnitude) {
          typedef struct AverageGradientMagnitude AvgOp;
          AverageGradientVectors<AvgOp>(dx, ndofs, edgeTable, _GradientAveraging);
        } else {
          typedef struct AverageGradient AvgOp;
          AverageGradientVectors<AvgOp>(dx, ndofs, edgeTable, _GradientAveraging);
        }
      }

    }
    MIRTK_DEBUG_TIMING(3, "averaging of"
        << (_AverageSignedGradients ? " signed " : " ")
        << "energy gradient" << (_AverageGradientMagnitude ? " magnitude " : " ")
        << "(#iter=" << _GradientAveraging << ")");
  }
}

// -----------------------------------------------------------------------------
// Enforce non-self-intersection / minimum cell distance constraints
void DeformableSurfaceModel
::ResolveSurfaceCollisions(double *dx, bool nsi, double mind, double minw) const
{
  MIRTK_START_TIMING();

  const int    fix_attempt      =  2; // No. of attempts displacement is halfed
  const int    max_attempt      = 10; // Maximum no. of attempts to fix collisions
  const double col_check_radius = 1.1; // Collided triangle bounding sphere radius factor
  const double min_check_radius = 1.1 * max(mind, minw); // Minimum collision re-check radius

  vtkPolyData  * const current = _PointSet.Surface();
  vtkDataArray * const status  = _PointSet.SurfaceStatus();

  vtkNew<vtkIdList> cellIds;
  vtkIdType         npts, *pts;
  double            a[3], b[3], c[3], p[3], bounds[6], r, *d, alpha;

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(current->GetNumberOfPoints());

  vtkSmartPointer<vtkDataArray> mask = vtkSmartPointer<vtkUnsignedCharArray>::New();
  mask->SetName("CollisionsMask");
  mask->SetNumberOfComponents(1);
  mask->SetNumberOfTuples(current->GetNumberOfCells());
  mask->FillComponent(0, 1.);

  vtkSmartPointer<vtkDataArray> scale = vtkSmartPointer<vtkFloatArray>::New();
  scale->SetName("DisplacementScale");
  scale->SetNumberOfComponents(1);
  scale->SetNumberOfTuples(current->GetNumberOfPoints());

  vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();
  surface->ShallowCopy(current);
  surface->SetPoints(points);
  surface->GetPointData()->AddArray(scale);
  surface->GetCellData ()->AddArray(mask);

  vtkNew<vtkCellLocator> locator;
  locator->SetDataSet(current);
  locator->BuildLocator();

  SurfaceCollisions check;
  check.AdjacentIntersectionTest(nsi);
  check.NonAdjacentIntersectionTest(nsi);
  check.FrontfaceCollisionTest(mind > .0);
  check.BackfaceCollisionTest(minw > .0);
  check.FastCollisionTest(_FastCollisionTest);
  check.MinFrontfaceDistance(mind);
  check.MinBackfaceDistance(minw);
  check.MaxAngle(_MaxCollisionAngle);
  check.MaxSearchRadius(2. * _MaxEdgeLength);
  check.StoreIntersectionDetailsOff();
  check.StoreCollisionDetailsOff();
  check.ResetCollisionTypeOn();
  check.Input(surface);
  check.Mask(mask);

  bool modified = true;
  for (int attempt = 1; modified && attempt <= max_attempt; ++attempt) {

    // Re-evaluate collisions only for cells near previously detected collisions
    if (attempt > 1) {
      mask->FillComponent(0, 0.);
      vtkDataArray * const collisions = check.GetCollisionTypeArray();
      for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); ++cellId) {
        if (collisions->GetComponent(cellId, 0) != 0.) {
          mask->SetComponent(cellId, 0, 1.);
          current->GetCellPoints(cellId, npts, pts);
          current->GetPoint(pts[0], a);
          current->GetPoint(pts[1], b);
          current->GetPoint(pts[2], c);
          r = max(min_check_radius, col_check_radius * Triangle::BoundingSphereRadius(a, b, c, p));
          bounds[0] = p[0] - r, bounds[1] = p[0] + r;
          bounds[2] = p[1] - r, bounds[3] = p[1] + r;
          bounds[4] = p[2] - r, bounds[5] = p[2] + r;
          locator->FindCellsWithinBounds(bounds, cellIds.GetPointer());
          for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
            mask->SetComponent(cellIds->GetId(i), 0, 1.);
          }
        }
      }
    }

    // Evaluate collisions only for cells with a non-zero vertex displacement
    // as we cannot fix collisions for these cells by scaling the displacements
    for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); ++cellId) {
      if (mask->GetComponent(cellId, 0) != 0.) {
        mask->SetComponent(cellId, 0, 0.);
        current->GetCellPoints(cellId, npts, pts);
        for (vtkIdType i = 0; i < npts; ++i) {
          d = dx + 3 * pts[i];
          if (d[0] != 0. || d[1] != 0. || d[2] != 0.) {
            mask->SetComponent(cellId, 0, 1.);
            break;
          }
        }
      }
    }

    // Update surface points and (re-)check for collisions
    MovePoints::Run(_PointSet.Points(), dx, points);
    check.Run();

    // Add copy of collisions array to deformable surface for debug output
    // (before actually fixing the collisions as none should be left afterwards)
    if (attempt == 1 && debug > 0) {
      vtkSmartPointer<vtkDataArray> collisions;
      collisions.TakeReference(check.GetCollisionTypeArray()->NewInstance());
      collisions->DeepCopy(check.GetCollisionTypeArray());
      collisions->SetName(check.GetCollisionTypeArray()->GetName());
      current->GetCellData()->RemoveArray(collisions->GetName());
      current->GetCellData()->AddArray(collisions);
    }

    // Set displacement scaling factor for points of collided cells
    alpha = (attempt <= fix_attempt ? .5 : 0.);
    scale->FillComponent(0, 1.);
    for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); ++cellId) {
      if (check.GetCollisionType(cellId) != SurfaceCollisions::NoCollision) {
        surface->GetCellPoints(cellId, npts, pts);
        for (vtkIdType i = 0; i < npts; ++i) {
          if (status == nullptr || status->GetComponent(pts[i], 0) != 0.) {
            scale->SetComponent(pts[i], 0, alpha);
          }
        }
      }
    }

    // Adjust magnitude of displacement vectors
    double *d = dx;
    modified  = false;
    for (vtkIdType ptId = 0; ptId < surface->GetNumberOfPoints(); ++ptId, d += 3) {
      alpha = max(scale->GetComponent(ptId, 0), 0.);
      if (alpha < 1. && (d[0] != 0. || d[1] != 0. || d[2] != 0.)) {
        d[0] *= alpha, d[1] *= alpha, d[2] *= alpha;
        modified = true;
      }
    }
  }

  MIRTK_DEBUG_TIMING(3, "resolving collisions");
}

// -----------------------------------------------------------------------------
void DeformableSurfaceModel::EnforceHardConstraints(double *dx) const
{
  // Hard constraints only apply to non-parametric deformable surface models
  if (_Transformation) return;

  // Hard surface mesh constraints
  if (_IsSurfaceMesh) {

    // Disallow movement of passive nodes
    if (_PointSet.SurfaceStatus() && _FixPassivePoints) {
      double *d = dx;
      vtkDataArray *status = _PointSet.SurfaceStatus();
      for (int ptId = 0; ptId < _PointSet.NumberOfSurfacePoints(); ++ptId, d += 3) {
        if (status->GetComponent(ptId, 0) == .0) {
          d[0] = d[1] = d[2] = 0.;
        }
      }
    }

    // Disallow expansion / contraction
    if (!_AllowExpansion || !_AllowContraction) {
      double n[3], dp, *d = dx;
      vtkDataArray *normals = _PointSet.SurfaceNormals();
      for (int ptId = 0; ptId < _PointSet.NumberOfSurfacePoints(); ++ptId, d += 3) {
        normals->GetTuple(ptId, n);
        dp = d[0]*n[0] + d[1]*n[1] + d[2]*n[2];
        if ((!_AllowExpansion && dp > 0.) || (!_AllowContraction && dp < 0.)) {
          d[0] = d[1] = d[2] = 0.;
        }
      }
    }

    // Disallow points to travel further away from input surface than the given threshold
    if (!IsInf(_MaxInputDistance)) {
      int subId;
      vtkIdType cellId;
      vtkNew<vtkGenericCell> cell;
      double p[3], x[3], dist2, *d = dx;
      const double max_dist2 = _MaxInputDistance * _MaxInputDistance;
      for (int ptId = 0; ptId < _PointSet.NumberOfSurfacePoints(); ++ptId, d += 3) {
        _PointSet.Surface()->GetPoint(ptId, p);
        p[0] += d[0], p[1] += d[1], p[2] += d[2];
        _InputCellLocator->FindClosestPoint(p, x, cell.GetPointer(), cellId, subId, dist2);
        if (dist2 > max_dist2) d[0] = d[1] = d[2] = 0.;
      }
    }

    // Enforce non-self-intersection / minimum cell distance constraints
    if (_HardNonSelfIntersection || _MinFrontfaceDistance > 0. || _MinBackfaceDistance > 0.) {
      ResolveSurfaceCollisions(dx, _HardNonSelfIntersection, _MinFrontfaceDistance, _MinBackfaceDistance);
    }
  }
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void DeformableSurfaceModel::WriteDataSets(const char *prefix, const char *suffix, bool all) const
{
  const int sz = 1024;
  char      fname[sz];

  snprintf(fname, sz, "%soutput%s%s", prefix, suffix, _PointSet.DefaultExtension());
  _PointSet.Write(fname);

  if (all) {
    for (int i = 0; i < _NumberOfTerms; ++i) {
      const EnergyTerm *term = Term(i);
      if (term->Weight() != .0) {
        term->WriteDataSets(prefix, suffix, all);
      }
    }
  }
}

// -----------------------------------------------------------------------------
void DeformableSurfaceModel::WriteGradient(const char *prefix, const char *suffix) const
{
  for (int i = 0; i < _NumberOfTerms; ++i) {
    const EnergyTerm *term = Term(i);
    if (term->Weight() != .0) {
      term->WriteGradient(prefix, suffix);
    }
  }
}


} // namespace mirtk
