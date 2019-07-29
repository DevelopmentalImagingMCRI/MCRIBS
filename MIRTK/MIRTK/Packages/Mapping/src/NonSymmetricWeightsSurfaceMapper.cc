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

#include "mirtk/NonSymmetricWeightsSurfaceMapper.h"

#include "Eigen/SparseCore"
#include "Eigen/SparseLU"
#include "Eigen/OrderingMethods"
#include "Eigen/IterativeLinearSolvers"


namespace mirtk {


// Global flags (cf. mirtk/Options.h)
MIRTK_Common_EXPORT extern int verbose;


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void NonSymmetricWeightsSurfaceMapper
::CopyAttributes(const NonSymmetricWeightsSurfaceMapper &other)
{
}

// -----------------------------------------------------------------------------
NonSymmetricWeightsSurfaceMapper::NonSymmetricWeightsSurfaceMapper()
{
}

// -----------------------------------------------------------------------------
NonSymmetricWeightsSurfaceMapper
::NonSymmetricWeightsSurfaceMapper(const NonSymmetricWeightsSurfaceMapper &other)
:
  LinearFixedBoundarySurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
NonSymmetricWeightsSurfaceMapper &NonSymmetricWeightsSurfaceMapper
::operator =(const NonSymmetricWeightsSurfaceMapper &other)
{
  if (this != &other) {
    LinearFixedBoundarySurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
NonSymmetricWeightsSurfaceMapper::~NonSymmetricWeightsSurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void NonSymmetricWeightsSurfaceMapper
::Weights(int i, const int *j, double *w_i, int d_i) const
{
  for (int k = 0; k < d_i; ++k) {
    w_i[k] = this->Weight(i, j[k]);
  }
}

// -----------------------------------------------------------------------------
double NonSymmetricWeightsSurfaceMapper::Weight(int i, int j) const
{
  cerr << this->NameOfClass() << "::Weight: Either this or Weights function must be implemented by subclass" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void NonSymmetricWeightsSurfaceMapper::ComputeMap()
{
  const bool use_direct_solver = (_NumberOfIterations < 0 || _NumberOfIterations == 1);

  typedef Eigen::MatrixXd             Values;
  typedef Eigen::SparseMatrix<double> Matrix;
  typedef Eigen::Triplet<double>      NZEntry;

  const int N = NumberOfPoints();
  const int n = NumberOfFreePoints();
  const int m = NumberOfComponents();

  int i, k, l, r, c;
  const int *j;

  Matrix A(n, n);
  Values b(n, m);
  {
    Array<NZEntry> w;
    Array<double>  w_ii(n, .0);

    b.setZero();
    w.reserve(2 * _EdgeTable->NumberOfEdges() + n);

    int     d_i;
    double *w_i = new double[_EdgeTable->MaxNumberOfAdjacentPoints()];

    EdgeIterator edgeIt(*_EdgeTable);
    for (i = 0; i < N; ++i) {
      r = FreePointIndex(i);
      if (r >= 0) {
        _EdgeTable->GetAdjacentPoints(i, d_i, j);
        this->Weights(i, j, w_i, d_i);
        for (k = 0; k < d_i; ++k) {
          c = FreePointIndex(j[k]);
          if (c >= 0) {
            w.push_back(NZEntry(r, c, -w_i[k]));
          } else {
            for (l = 0; l < m; ++l) {
              b(r, l) += w_i[k] * GetValue(j[k], l);
            }
          }
          w_ii[r] += w_i[k];
        }
      }
    }
    for (r = 0; r < n; ++r) {
      w.push_back(NZEntry(r, r, w_ii[r]));
    }
    delete[] w_i;

    A.setFromTriplets(w.begin(), w.end());
  }

  if (verbose) {
    cout << "\n";
    cout << "  No. of surface points        = " << NumberOfPoints() << "\n";
    cout << "  No. of free points           = " << n << "\n";
    cout << "  No. of non-zero coefficients = " << A.nonZeros() << "\n";
    cout << "  Dimension of map codomain    = " << m << "\n";
    cout.flush();
  }

  Values x(n, m);
  int    niter = 0;
  double error = nan;

  if (use_direct_solver) {
    Eigen::SparseLU<Matrix> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    x = solver.solve(b);
  } else {
    Eigen::BiCGSTAB<Matrix> solver(A);
    if (_NumberOfIterations > 0 ) solver.setMaxIterations(_NumberOfIterations);
    if (_Tolerance          > 0.) solver.setTolerance(_Tolerance);
    for (r = 0; r < n; ++r) {
      i = FreePointId(r);
      for (l = 0; l < m; ++l) {
        x(r, l) = GetValue(i, l);
      }
    }
    x = solver.solveWithGuess(b, x);
    niter = solver.iterations();
    error = solver.error();
  }

  for (r = 0; r < n; ++r) {
    i = FreePointId(r);
    for (l = 0; l < m; ++l) {
      SetValue(i, l, x(r, l));
    }
  }

  if (verbose && !use_direct_solver) {
    cout << "  No. of iterations            = " << niter << "\n";
    cout << "  Estimated error              = " << error << "\n";
    cout.flush();
  }
}


} // namespace mirtk
