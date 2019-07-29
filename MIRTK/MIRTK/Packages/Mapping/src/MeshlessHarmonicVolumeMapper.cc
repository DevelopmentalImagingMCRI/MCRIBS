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

#include "mirtk/MeshlessHarmonicVolumeMapper.h"
#include "mirtk/MeshlessHarmonicMap.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Array.h"
#include "mirtk/PointSet.h"
#include "mirtk/Matrix.h"
#include "mirtk/Parallel.h"
#include "mirtk/VtkMath.h"

#include "mirtk/Eigen.h"
#include "Eigen/SVD"


namespace mirtk {


// Global flags (cf. mirtk/Options.h)
MIRTK_Common_EXPORT extern int verbose;


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace MeshlessHarmonicVolumeMapperUtils {

// -----------------------------------------------------------------------------
/// Compute A = K^T K
struct ComputeCoefficients
{
  const Matrix     *_Kernel;
  const Array<int> *_ColIdx;
  Matrix           *_Coeffs;

  void operator ()(const blocked_range2d<int> &re) const
  {
    double sum;
    const double *ci, *cj;
    for (int j = re.cols().begin(); j < re.cols().end(); ++j)
    for (int i = re.rows().begin(); i < re.rows().end(); ++i) {
      sum = .0;
      ci = _Kernel->RawPointer(0, (*_ColIdx)[i]);
      cj = _Kernel->RawPointer(0, (*_ColIdx)[j]);
      for (int r = 0; r < _Kernel->Rows(); ++r, ++ci, ++cj) {
        sum += (*ci) * (*cj);
      }
      _Coeffs->Put(i, j, sum);
    }
  }
};

// -----------------------------------------------------------------------------
/// Compute b = K^T f
struct ComputeConstraints
{
  const Matrix     *_Kernel;
  const Array<int> *_ColIdx;
  vtkDataArray     *_BoundaryMap;
  Matrix           *_b;

  void operator ()(const blocked_range2d<int> &re) const
  {
    double sum;
    const double *ci;
    for (int j = re.cols().begin(); j < re.cols().end(); ++j)
    for (int i = re.rows().begin(); i < re.rows().end(); ++i) {
      sum = .0;
      ci = _Kernel->RawPointer(0, (*_ColIdx)[i]);
      for (int r = 0; r < _Kernel->Rows(); ++r, ++ci) {
        sum += (*ci) * _BoundaryMap->GetComponent(r, j);
      }
      _b->Put(i, j, sum);
    }
  }
};


} // namespace MeshlessHarmonicVolumeMapperUtils
using namespace MeshlessHarmonicVolumeMapperUtils;

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void MeshlessHarmonicVolumeMapper
::CopyAttributes(const MeshlessHarmonicVolumeMapper &other)
{
  _Kernel = other._Kernel;
  _UseSVD = other._UseSVD;
}

// -----------------------------------------------------------------------------
MeshlessHarmonicVolumeMapper::MeshlessHarmonicVolumeMapper()
:
  _UseSVD(false)
{
}

// -----------------------------------------------------------------------------
MeshlessHarmonicVolumeMapper
::MeshlessHarmonicVolumeMapper(const MeshlessHarmonicVolumeMapper &other)
:
  MeshlessVolumeMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
MeshlessHarmonicVolumeMapper &MeshlessHarmonicVolumeMapper
::operator =(const MeshlessHarmonicVolumeMapper &other)
{
  if (this != &other) {
    MeshlessVolumeMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
MeshlessHarmonicVolumeMapper::~MeshlessHarmonicVolumeMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void MeshlessHarmonicVolumeMapper::Initialize()
{
  // Initialize base class
  MeshlessVolumeMapper::Initialize();

  const int m = NumberOfBoundaryPoints();
  const int n = NumberOfSourcePoints();
  const int d = NumberOfComponents();

  // Initialize harmonic map and precompute kernel function values
  double p[3], q[3], dist;

  SharedPtr<MeshlessHarmonicMap> map = NewShared<MeshlessHarmonicMap>();

  PointSet &points  = map->SourcePoints();
  Matrix   &weights = map->Coefficients();

  points.Resize(n);
  weights.Initialize(n, d);
  _Kernel.Initialize(m, n);

  double *c = _Kernel.RawPointer();
  for (int j = 0; j < n; ++j) {
    _OffsetSurface->GetPoint(j, q);
    points.SetPoint(j, q);
    for (int i = 0; i < m; ++i, ++c) {
      _Boundary->GetPoint(i, p);
      dist = sqrt(vtkMath::Distance2BetweenPoints(p, q));
      *c = MeshlessHarmonicMap::H(dist);
    }
  }

  // Set output map
  _Output = map;
}

// -----------------------------------------------------------------------------
bool MeshlessHarmonicVolumeMapper::AddSourcePoint(double q[3])
{
  if (!MeshlessVolumeMapper::AddSourcePoint(q)) return false;

  const int m = NumberOfBoundaryPoints();
  const int n = NumberOfSourcePoints();

  _Kernel.Resize(m, n);

  double p[3], dist, *c = _Kernel.RawPointer(0, n - 1);
  for (int i = 0; i < m; ++i, ++c) {
    _Boundary->GetPoint(i, p);
    dist = sqrt(vtkMath::Distance2BetweenPoints(p, q));
    *c = MeshlessHarmonicMap::H(dist);
  }

  return true;
}

// -----------------------------------------------------------------------------
void MeshlessHarmonicVolumeMapper::Solve()
{
  // Use base class implementation if possible
  if (!_UseSVD) {
    MeshlessVolumeMapper::Solve();
    return;
  }

  typedef Eigen::JacobiSVD<Eigen::MatrixXd> SVD;
  typedef SVD::SingularValuesType           SingularValues;

  const int m = NumberOfBoundaryPoints();
  const int d = NumberOfComponents();

  Matrix         A, b(m, d), x; // coefficients matrix, right-hand side, solution of Ax = b
  SingularValues sigma;         // singular values of coefficients matrix
  double         error;         // error of boundary map approximation
  double         min_error, max_error, std_error;
  double         p[3], q[3];

  // Compute initial error
  if (verbose) {
    cout << "Initialize residual boundary map...", cout.flush();
  }
  error = this->UpdateResidualMap(&min_error, &max_error, &std_error);
  if (verbose) {
    cout << " done\nBoundary fitting error (MSE) = " << error
         << " (+/-" << std_error << "), range = ["
         << min_error << ", " << max_error << "]" << endl;
  }

  // Iteratively approximate volumetric map
  double *df = new double[d];
  for (int iter = 0; iter < _NumberOfIterations; ++iter) {

    if (verbose) cout << "\nIteration " << (iter+1) << endl;

    // Evenly partition set of source points
    this->PartitionSourcePoints();

    // Perform boundary fitting for each subset
    for (int k = 0; k < NumberOfSourcePointSets(); ++k) {
      const int n = NumberOfSourcePoints(k);

      if (verbose) {
        cout << "Use source points subset " << (k+1) << " out of "
             << NumberOfSourcePointSets() << " with " << n << " points" << endl;
      }

      // Get coefficients
      A.Initialize(m, n);
      for (int j = 0, c; j < n; ++j) {
        c = SourcePointIndex(k, j);
        for (int i = 0; i < m; ++i) {
          A(i, j) = _Kernel(i, c);
        }
      }

      // Get right-hand side
      double *c = b.RawPointer();
      for (int j = 0; j < d; ++j)
      for (int i = 0; i < m; ++i, ++c) {
        *c = _ResidualMap->GetComponent(i, j);
      }

      // Solve linear system using SVD
      if (verbose) cout << "Solve linear system using SVD...", cout.flush();
      SVD svd(MatrixToEigen(A), Eigen::ComputeThinU | Eigen::ComputeThinV);
      x = EigenToMatrix(svd.solve(MatrixToEigen(b)));
      sigma = svd.singularValues();
      if (verbose) {
        cout << " done\nmax(sigma) = " << sigma(0)
             << ", min(sigma) = " << sigma(sigma.size() - 1)
             << ", cond(A) = " << (sigma(0) / sigma(sigma.size() - 1)) << endl;
      }

      // Add solution to volumetric map
      if (verbose) cout << "Add solution to harmonic map...", cout.flush();
      this->AddWeights(k, x);
      if (verbose) cout << " done" << endl;

      // Update residual boundary map
      if (verbose) cout << "Update residual boundary map...", cout.flush();
      error = this->UpdateResidualMap(&min_error, &max_error, &std_error);
      if (verbose) {
        cout << " done\n";
        cout << "Boundary fitting error (MSE) = " << error
             << " (+/-" << std_error << "), range = ["
             << min_error << ", " << max_error << "]" << endl;
      }

      // TODO: Adapt set of source points (Xu et al., 2013)
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
}

// =============================================================================
// Linear system
// =============================================================================

// -----------------------------------------------------------------------------
void MeshlessHarmonicVolumeMapper
::GetCoefficients(int k, Matrix &coeffs) const
{
  const int n = NumberOfSourcePoints(k);
  coeffs.Initialize(n, n);
  ComputeCoefficients eval;
  eval._Kernel = &_Kernel;
  eval._ColIdx = &_SourcePartition[k];
  eval._Coeffs = &coeffs;
  parallel_for(blocked_range2d<int>(0, n, 0, n), eval);
}

// -----------------------------------------------------------------------------
void MeshlessHarmonicVolumeMapper
::GetConstraints(int k, Matrix &b) const
{
  const int n = NumberOfSourcePoints(k);
  const int d = NumberOfComponents();
  b.Initialize(n, d);
  ComputeConstraints eval;
  eval._Kernel      = &_Kernel;
  eval._ColIdx      = &_SourcePartition[k];
  eval._BoundaryMap = _ResidualMap;
  eval._b           = &b;
  parallel_for(blocked_range2d<int>(0, n, 0, d), eval);
}

// -----------------------------------------------------------------------------
void MeshlessHarmonicVolumeMapper
::AddWeights(int k, const Matrix &w)
{
  const int d = NumberOfComponents();
  MeshlessHarmonicMap *map = dynamic_cast<MeshlessHarmonicMap *>(_Output.get());
  Matrix &weights = map->Coefficients();
  for (int i = 0; i < NumberOfSourcePoints(k); ++i) {
    const int r = SourcePointIndex(k, i);
    for (int j = 0; j < d; ++j) {
      weights(r, j) += w(i, j);
    }
  }
}


} // namespace mirtk
