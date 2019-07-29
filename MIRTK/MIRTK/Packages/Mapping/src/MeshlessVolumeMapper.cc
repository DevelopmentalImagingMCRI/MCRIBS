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

#include "mirtk/MeshlessVolumeMapper.h"

#include "mirtk/CommonExport.h"

#include "mirtk/Math.h"
#include "mirtk/Assert.h"
#include "mirtk/Parallel.h"
#include "mirtk/MeshSmoothing.h"
#include "mirtk/PointSetIO.h"

#include "mirtk/Vtk.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkImplicitModeller.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkContourFilter.h"
#include "vtkLinearSubdivisionFilter.h"
#include "vtkQuadricDecimation.h"
#include "vtkAbstractCellLocator.h"
#include "vtkCellLocator.h"
#include "vtkGenericCell.h"

#include "mirtk/Eigen.h"
#include "Eigen/LU"
#include "Eigen/SVD"


namespace mirtk {


// Global flags (cf. mirtk/Options.h)
MIRTK_Common_EXPORT extern int verbose;
MIRTK_Common_EXPORT extern int debug;


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace MeshlessVolumeMapperUtils {


// -----------------------------------------------------------------------------
/// Compute df = f - \sum f_i
struct ComputeResidualMap
{
  vtkPoints     *_BoundarySet;
  vtkDataArray  *_BoundaryMap;
  Mapping       *_OutputMap;
  vtkDataArray  *_ResidualMap;
  int            _OutputDimension;
  double         _SquaredError;
  double         _SquaredError2;
  double         _MinSquaredError;
  double         _MaxSquaredError;

  ComputeResidualMap()
  :
    _SquaredError(.0),
    _SquaredError2(.0),
    _MinSquaredError(numeric_limits<double>::infinity()),
    _MaxSquaredError(.0)
  {}

  ComputeResidualMap(const ComputeResidualMap &other, split)
  :
    _BoundarySet    (other._BoundarySet),
    _BoundaryMap    (other._BoundaryMap),
    _OutputMap      (other._OutputMap),
    _ResidualMap    (other._ResidualMap),
    _OutputDimension(other._OutputDimension),
    _SquaredError   (.0),
    _SquaredError2  (.0),
    _MinSquaredError(other._MinSquaredError),
    _MaxSquaredError(other._MaxSquaredError)
  {}

  void join(const ComputeResidualMap &other)
  {
    _SquaredError   += other._SquaredError;
    _SquaredError2  += other._SquaredError2;
    _MinSquaredError = min(_MinSquaredError, other._MinSquaredError);
    _MaxSquaredError = max(_MaxSquaredError, other._MaxSquaredError);
  }

  inline double Dot(double *v, double *w)
  {
    double dot = .0;
    for (int i = 0; i < _OutputDimension; ++i) {
      dot += v[i] * w[i];
    }
    return dot;
  }

  void operator ()(const blocked_range<vtkIdType> &re)
  {
    double p[3], dist2;
    double *f  = new double[_OutputDimension];
    double *df = new double[_OutputDimension];
    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      _BoundarySet->GetPoint(ptId, p);
      _OutputMap->Evaluate(f, p);
      for (int i = 0; i < _OutputDimension; ++i) {
        df[i] = _BoundaryMap->GetComponent(ptId, i) - f[i];
      }
      _ResidualMap->SetTuple(ptId, df);
      dist2 = Dot(df, df);
      _SquaredError  += dist2;
      _SquaredError2 += dist2 * dist2;
      if (dist2 < _MinSquaredError) _MinSquaredError = dist2;
      if (dist2 > _MaxSquaredError) _MaxSquaredError = dist2;
    }
    delete[] f;
    delete[] df;
  }
};


} // namespace MeshlessVolumeMapperUtils
using namespace MeshlessVolumeMapperUtils;

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void MeshlessVolumeMapper
::CopyAttributes(const MeshlessVolumeMapper &other)
{
  _BoundaryPointsRatio         = other._BoundaryPointsRatio;
  _SourcePointsRatio           = other._SourcePointsRatio;
  _MaximumNumberOfSourcePoints = other._MaximumNumberOfSourcePoints;
  _NumberOfIterations          = other._NumberOfIterations;
  _ImplicitSurfaceSize         = other._ImplicitSurfaceSize;
  _ImplicitSurfaceSpacing      = other._ImplicitSurfaceSpacing;
  _DistanceOffset              = other._DistanceOffset;
  _MaximumConditionNumber      = other._MaximumConditionNumber;
  _OffsetSurface               = other._OffsetSurface;
  _SourcePartition             = other._SourcePartition;
}

// -----------------------------------------------------------------------------
MeshlessVolumeMapper::MeshlessVolumeMapper()
:
  _BoundaryPointsRatio(.1),
  _SourcePointsRatio(.1),
  _MaximumNumberOfSourcePoints(500),
  _NumberOfIterations(10),
  _ImplicitSurfaceSize(0),
  _ImplicitSurfaceSpacing(.0),
  _DistanceOffset(-.1),
  _MaximumConditionNumber(1.0e6)
{
}

// -----------------------------------------------------------------------------
MeshlessVolumeMapper
::MeshlessVolumeMapper(const MeshlessVolumeMapper &other)
:
  VolumeMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
MeshlessVolumeMapper &MeshlessVolumeMapper
::operator =(const MeshlessVolumeMapper &other)
{
  if (this != &other) {
    VolumeMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
MeshlessVolumeMapper::~MeshlessVolumeMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void MeshlessVolumeMapper
::GetClosestPointOnOffsetSurface(double x[3], double p[3])
{
  vtkIdType cellId;
  int       subId;
  double    dist2;
  _OffsetPointLocator->FindClosestPoint(x, p, cellId, subId, dist2);
}

// -----------------------------------------------------------------------------
void MeshlessVolumeMapper::UpdateBoundary(vtkPolyData *boundary)
{
  // Initialize new boundary map
  vtkSmartPointer<vtkDataArray> boundary_map;
  boundary_map.TakeReference(_BoundaryMap->NewInstance());
  boundary_map->SetName(_BoundaryMap->GetName());
  boundary_map->SetNumberOfComponents(_BoundaryMap->GetNumberOfComponents());
  boundary_map->SetNumberOfTuples(boundary->GetNumberOfPoints());

  boundary->GetPointData()->Initialize();
  boundary->GetPointData()->AddArray(boundary_map);

  // Initialize old boundary surface cell locator
  vtkSmartPointer<vtkAbstractCellLocator> locator;
  locator = vtkSmartPointer<vtkCellLocator>::New();
  locator->SetDataSet(_Boundary);
  locator->BuildLocator();

  // Project new boundary points onto old boundary surface and interpolate
  // boundary map value at projected point
  double    x[3], p[3], pcoords[3], dist2;
  double   *weights = new double[_Boundary->GetMaxCellSize()];
  double   *pmap    = new double[_BoundaryMap->GetNumberOfComponents()];
  vtkIdType cellId, cellPtId;
  int       subId;

  vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
  for (vtkIdType ptId = 0; ptId < boundary->GetNumberOfPoints(); ++ptId) {
    boundary->GetPoint(ptId, x);
    locator->FindClosestPoint(x, p, cell, cellId, subId, dist2);
    cell->EvaluatePosition(p, NULL, subId, pcoords, dist2, weights);
    memset(pmap, 0, _BoundaryMap->GetNumberOfComponents() * sizeof(double));
    for (vtkIdType i = 0; i < cell->GetNumberOfPoints(); ++i) {
      cellPtId = cell->GetPointId(i);
      for (int j = 0; j < _BoundaryMap->GetNumberOfComponents(); ++j) {
        pmap[j] += weights[i] * _BoundaryMap->GetComponent(cellPtId, j);
      }
    }
    boundary_map->SetTuple(ptId, pmap);
  }
  delete[] weights;
  delete[] pmap;

  // Replace boundary surface and boundary map
  _Boundary    = boundary;
  _BoundaryMap = boundary_map;
}

// -----------------------------------------------------------------------------
void MeshlessVolumeMapper::PlaceBoundaryPoints()
{
  mirtkAssert(_Boundary != NULL, "input boundary surface is initialized");

  if (verbose) {
    cout << "Place boundary points...", cout.flush();
  }

  // Discard all data arrays
  vtkSmartPointer<vtkPolyData> boundary;
  boundary.TakeReference(_Boundary->NewInstance());
  boundary->ShallowCopy(_Boundary);
  boundary->GetPointData()->Initialize();
  boundary->GetCellData ()->Initialize();
  boundary->GetFieldData()->Initialize();

  if (_BoundaryMap->GetNumberOfComponents() == 3) {
    boundary->GetPointData()->SetTCoords(_BoundaryMap);
  }

  // Subdivide input surface to generate dense surface mesh
  vtkSmartPointer<vtkLinearSubdivisionFilter> subdiv;
  subdiv = vtkSmartPointer<vtkLinearSubdivisionFilter>::New();
  subdiv->SetNumberOfSubdivisions(1);
  SetVTKInput(subdiv, boundary);

  // Calculate target reduction of dense surface mesh
  subdiv->Update();
  double n     = _Boundary->GetNumberOfPoints();
  double m     = subdiv->GetOutput()->GetNumberOfPoints();
  double ratio = _BoundaryPointsRatio * n / m;

  // Decimate dense surface mesh while minimizing the quadric error metric
  vtkSmartPointer<vtkQuadricDecimation> decimate;
  decimate = vtkSmartPointer<vtkQuadricDecimation>::New();
  decimate->SetTargetReduction(1.0 - ratio);
  if (boundary->GetPointData()->GetTCoords()) {
    decimate->AttributeErrorMetricOn();
    decimate->TCoordsAttributeOn();
    decimate->SetTCoordsWeight(.1);
  }
  SetVTKConnection(decimate, subdiv);

  decimate->Update();

  this->UpdateBoundary(decimate->GetOutput());

  if (debug) WritePolyData("boundary_surface.vtp", _Boundary);

  if (verbose) {
    cout << " done: N_c = " << _Boundary->GetNumberOfPoints() << endl;
  }
}

// -----------------------------------------------------------------------------
void MeshlessVolumeMapper::PlaceSourcePoints()
{
  mirtkAssert(_Boundary != NULL, "input boundary surface is initialized");

  if (verbose) {
    cout << "Place source points...", cout.flush();
  }

  // Compute input surface bounds
  double bounds[6];
  _Boundary->GetBounds(bounds);

  // Offset surface distance
  double offset = _DistanceOffset;
  if (offset < .0) {
    offset *= -sqrt(pow(bounds[1] - bounds[0], 2) +
                    pow(bounds[3] - bounds[2], 2) +
                    pow(bounds[5] - bounds[4], 2));
  }

  // Adjust implicit surface model bounds
  const double margin = 1.15 * offset;
  bounds[0] -= margin;
  bounds[1] += margin;
  bounds[2] -= margin;
  bounds[3] += margin;
  bounds[4] -= margin;
  bounds[5] += margin;

  // Sampling of implicit surface model
  int nx = _ImplicitSurfaceSize;
  int ny = _ImplicitSurfaceSize;
  int nz = _ImplicitSurfaceSize;

  if (nx <= 0 || ny <= 0 || nz <= 0) {
    double dx = _ImplicitSurfaceSpacing;
    double dy = _ImplicitSurfaceSpacing;
    double dz = _ImplicitSurfaceSpacing;
    if (dx <= .0 && nx > 0) dx = (bounds[1] - bounds[0]) / nx;
    if (dy <= .0 && ny > 0) dy = (bounds[3] - bounds[2]) / ny;
    if (dz <= .0 && nz > 0) dz = (bounds[5] - bounds[4]) / nz;
    double ds = min(min(dx, dy), dz);
    if (ds <= .0) {
      ds = sqrt(pow(bounds[1] - bounds[0], 2) +
                pow(bounds[3] - bounds[2], 2) +
                pow(bounds[5] - bounds[4], 2)) / 128;
    }
    if (dx <= .0) dx = ds;
    if (dy <= .0) dy = ds;
    if (dz <= .0) dz = ds;
    if (nx <=  0) nx = int(ceil((bounds[1] - bounds[0]) / dx));
    if (ny <=  0) ny = int(ceil((bounds[3] - bounds[2]) / dy));
    if (nz <=  0) nz = int(ceil((bounds[5] - bounds[4]) / dz));
  }

  // Create implicit surface model from input surface
  vtkSmartPointer<vtkImplicitModeller> model;
  model = vtkSmartPointer<vtkImplicitModeller>::New();
  model->SetSampleDimensions(nx, ny, nz);
  model->SetMaximumDistance(1.1 * offset);
  model->SetModelBounds(bounds);
  model->AdjustBoundsOff();
  model->ScaleToMaximumDistanceOff();
  model->SetOutputScalarTypeToFloat();
  model->CappingOn();
  model->SetCapValue(margin);
  model->SetProcessModeToPerVoxel();
  SetVTKInput(model, _Boundary);

  // Extract inside/outside offset surfaces
  // Note: Unfortunately, vtkImplicitModeller creates an unsigned distance field.
  vtkSmartPointer<vtkContourFilter> contours;
  contours = vtkSmartPointer<vtkContourFilter>::New();
  contours->UseScalarTreeOn();
  contours->SetNumberOfContours(1);
  contours->SetValue(0, offset);
  SetVTKConnection(contours,  model);

  // Only keep offset surface closest to bounding box corner (i.e., outside)
  vtkSmartPointer<vtkPolyDataConnectivityFilter> outside;
  outside = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
  outside->SetClosestPoint(bounds[0], bounds[2], bounds[4]);
  outside->SetExtractionModeToClosestPointRegion();
  SetVTKConnection(outside, contours);

  // Subdivide offset surface to generate dense surface mesh
  vtkSmartPointer<vtkLinearSubdivisionFilter> subdiv;
  subdiv = vtkSmartPointer<vtkLinearSubdivisionFilter>::New();
  subdiv->SetNumberOfSubdivisions(1);
  SetVTKConnection(subdiv, outside);

  // Execute dense offset surface mesh generation
  subdiv->Update();

  // Smooth offset surface to reduce sampling artifacts
  MeshSmoothing smoother;
  smoother.Input(outside->GetOutput());
  smoother.NumberOfIterations(5);
  smoother.Weighting(MeshSmoothing::Gaussian);
  smoother.SmoothPoints(true);
  smoother.AdjacentValuesOnlyOn();
  smoother.Run();

  // Calculate target reduction of dense surface mesh
  double n     = _Boundary->GetNumberOfPoints();
  double m     = smoother.Output()->GetNumberOfPoints();
  double ratio = _SourcePointsRatio * n / m;

  // Decimate dense surface mesh while minimizing the quadric error metric
  vtkSmartPointer<vtkQuadricDecimation> decimate;
  decimate = vtkSmartPointer<vtkQuadricDecimation>::New();
  decimate->SetTargetReduction(1.0 - ratio);
  decimate->AttributeErrorMetricOff();
  SetVTKInput(decimate, smoother.Output());

  decimate->Update();
  _OffsetSurface = decimate->GetOutput();

  // Initialize offset surface point locator
  _OffsetPointLocator = vtkSmartPointer<vtkCellLocator>::New();
  _OffsetPointLocator->SetDataSet(_OffsetSurface);
  _OffsetPointLocator->BuildLocator();

  if (debug) WritePolyData("offset_surface.vtp", _OffsetSurface);

  if (verbose) {
    cout << " done: N_s = " << _OffsetSurface->GetNumberOfPoints() << endl;
  }
}

// -----------------------------------------------------------------------------
bool MeshlessVolumeMapper::AddSourcePoint(double q[3])
{
  MeshlessMap *map = dynamic_cast<MeshlessMap *>(_Output.get());
  return map->AddSourcePoint(q, 1e-9);
}

// -----------------------------------------------------------------------------
void MeshlessVolumeMapper::PartitionSourcePoints()
{
  int nsubsets = 1;
  if (_MaximumNumberOfSourcePoints > 0) {
    nsubsets = int(ceil(double(NumberOfSourcePoints()) / _MaximumNumberOfSourcePoints));
  }
  _SourcePartition.resize(nsubsets);
  if (nsubsets > 0) {
    int npoints = NumberOfSourcePoints() / nsubsets;
    for (int k = 0; k < nsubsets; ++k) {
      _SourcePartition[k].clear();
      _SourcePartition[k].reserve(npoints);
    }
    int i = 0;
    do {
      for (int k = 0; k < nsubsets && i < NumberOfSourcePoints(); ++k, ++i) {
        _SourcePartition[k].push_back(i);
      }
    } while (i < NumberOfSourcePoints());
  }
}

// -----------------------------------------------------------------------------
void MeshlessVolumeMapper::InitializeResidualMap()
{
  _ResidualMap.TakeReference(_BoundaryMap->NewInstance());
  _ResidualMap->DeepCopy(_BoundaryMap);
  _ResidualMap->SetName("ResidualMap");
  _Boundary->GetPointData()->AddArray(_ResidualMap);
}

// -----------------------------------------------------------------------------
double MeshlessVolumeMapper::UpdateResidualMap(double *min, double *max, double *std)
{
  ComputeResidualMap eval;
  eval._BoundarySet     = _Boundary->GetPoints();
  eval._BoundaryMap     = _BoundaryMap;
  eval._OutputMap       = _Output.get();
  eval._ResidualMap     = _ResidualMap;
  eval._OutputDimension = _Output->NumberOfComponents();
  parallel_reduce(blocked_range<vtkIdType>(0, _Boundary->GetNumberOfPoints()), eval);
  double avg = eval._SquaredError / NumberOfBoundaryPoints();
  if (min) *min = eval._MinSquaredError;
  if (max) *max = eval._MaxSquaredError;
  if (std) *std = sqrt(eval._SquaredError2 / NumberOfBoundaryPoints() - avg * avg);
  return avg;
}

// -----------------------------------------------------------------------------
void MeshlessVolumeMapper::Initialize()
{
  if (verbose) cout << endl;

  // Initialize base class
  VolumeMapper::Initialize();

  // Initialize boundary surface
  this->InitializeBoundary(_InputSet, _InputMap);

  // Sample source points from offset surface
  this->PlaceSourcePoints();

  // Sample boundary points from boundary surface
  this->PlaceBoundaryPoints();

  // Initialize residual boundary map
  this->InitializeResidualMap();

  if (_MaximumNumberOfSourcePoints <= 0) {
    _MaximumNumberOfSourcePoints = NumberOfBoundaryPoints() + 1;
  }
}

// -----------------------------------------------------------------------------
void MeshlessVolumeMapper::Solve()
{
  typedef Eigen::JacobiSVD<Eigen::MatrixXd> SVD;
  typedef SVD::SingularValuesType           SingularValues;

  Matrix         A, b, x; // coefficients matrix, right-hand side, solution of Ax = b
  SingularValues sigma;   // singular values
  double         alpha;   // weight of regularization term
  double         error;   // error of boundary map approximation
  double         min_error, max_error, std_error;
  double         p[3], q[3], *df;

  const int d = NumberOfComponents();

  df    = new double[d];
  alpha = .015;

  // Compute initial error
  if (verbose) {
    cout << "Initialize residual boundary map...";
    cout.flush();
  }
  error = this->UpdateResidualMap(&min_error, &max_error, &std_error);
  if (verbose) {
    cout << " done\nBoundary fitting error (MSE) = " << error
         << " (+/-" << std_error << "), range = ["
         << min_error << ", " << max_error << "]" << endl;
  }

  // Iteratively approximate volumetric map
  for (int iter = 0; iter < _NumberOfIterations; ++iter) {

    if (verbose) cout << "\nIteration " << (iter+1) << endl;

    // Evenly partition set of source points
    this->PartitionSourcePoints();

    // Perform boundary fitting for each subset
    for (int k = 0; k < NumberOfSourcePointSets(); ++k) {

      if (verbose) {
        cout << "Source points subset " << (k+1) << " out of " << NumberOfSourcePointSets() << endl;
      }

      // Get linear system of unregularized boundary fitting problem
      this->GetCoefficients(k, A);
      this->GetConstraints (k, b);

      // This generic implementation minimizes an energy function
      //
      //   E = w^T A w - b^T w + c
      //
      // where A = K^T K is a square matrix of size (k N_s) x (k N_s), where
      // k N_s is an integer multiple k of the number of source points, N_s,
      // and K is the matrix containing the sum of the kernel function weights
      // for each pair of source and boundary points. In particular, in case
      // of the harmonic map, K_ij = H(q_i, p_j), and in case of the biharmonic
      // map, K_ij = H(q_i, p_j) + dH(q_i, p_j) for 0 <= j < N_s and
      // K_ij = B(q_i, p_j) + dB(q_i, p_j) for N_s <= j < 2 N_s,
      // and i is the boundary point / constraint index.
      //
      // In case of the harmonic map, a different linear system with A = K
      // can be solved instead, using the (truncated or randomized) SVD as in
      // (Li et al., 2010). This alternative (slower!) method is implemented by
      // MeshlessHarmonicVolumeMapper::Parameterize for comparison.
      mirtkAssert(A.Rows() == A.Cols(), "coefficients matrix is square");
      mirtkAssert(b.Rows() == A.Rows(), "right-hand side has required number of rows");

      // Compute singular values of coefficients matrix
      if (alpha == .0) {
        if (verbose) cout << "Compute singular values...", cout.flush();
        SVD svd(MatrixToEigen(A));
        sigma = svd.singularValues();
        alpha = sigma(0) / (_MaximumConditionNumber - 1.0);
        if (verbose) {
          cout << " done\nmax(sigma) = " << sigma(0)
               << ", alpha = " << alpha
               << ", cond(A) = " << ((alpha + sigma(0)) / alpha) << endl;
        }
      }

      // Add regularization weight to diagonal matrix elements
      for (int i = 0; i < A.Rows(); ++i) {
        A(i, i) += alpha;
      }

      // Solve linear system using LU decomposition
      if (verbose) cout << "Solve linear system using LU decomposition...", cout.flush();
      x = EigenToMatrix(MatrixToEigen(A).partialPivLu().solve(MatrixToEigen(b)));
      if (verbose) cout << " done" << endl;

      // Add solution to volumetric map
      if (verbose) cout << "Add solution to harmonic map...", cout.flush();
      this->AddWeights(k, x);
      if (verbose) cout << " done" << endl;

      // Update residual boundary map
      if (verbose) cout << "Update residual boundary map...", cout.flush();
      error = this->UpdateResidualMap(&min_error, &max_error, &std_error);
      if (verbose) {
        cout << " done" << endl;
        cout << "Boundary fitting error (MSE) = " << error
             << " (+/-" << std_error << "), range = ["
             << min_error << ", " << max_error << "]" << endl;
      }

      // TODO: Remove source points with insignificant contribution
      //       if possible as done in (Li et al., 2010) and also mentioned
      //       in (Xu et al., 2013). This, however, is based on the singular
      //       values associated with each source point and thus requires
      //       an expensive SVD computation.
    }

    // Insert new source points by projecting boundary points with
    // high residual error onto the offset surface (cf. Xu et al., 2013)
    if (verbose) cout << "Insert new source points...", cout.flush();
    const int n = NumberOfSourcePoints();
    for (vtkIdType ptId = 0; ptId < _Boundary->GetNumberOfPoints(); ++ptId) {
      _ResidualMap->GetTuple(ptId, df);
      if ((df[0]*df[0] + df[1]*df[1] + df[2]*df[2]) > (error + 1.5 * std_error)) {
        _Boundary->GetPoint(ptId, p);
        GetClosestPointOnOffsetSurface(p, q);
        this->AddSourcePoint(q);
      }
    }
    if (verbose) {
      cout << " done: #points = " << NumberOfSourcePoints()
           << " (+" << (NumberOfSourcePoints() - n) << ")" << endl;
    }
  }

  delete[] df;

  if (debug) WritePolyData("boundary_surface.vtp", _Boundary);
}


} // namespace mirtk
