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

#include "mirtk/SymmetricWeightsSurfaceMapper.h"

#include "mirtk/EdgeTable.h"

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
void SymmetricWeightsSurfaceMapper
::CopyAttributes(const SymmetricWeightsSurfaceMapper &other)
{
}

// -----------------------------------------------------------------------------
SymmetricWeightsSurfaceMapper::SymmetricWeightsSurfaceMapper()
{
}

// -----------------------------------------------------------------------------
SymmetricWeightsSurfaceMapper
::SymmetricWeightsSurfaceMapper(const SymmetricWeightsSurfaceMapper &other)
:
  LinearFixedBoundarySurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
SymmetricWeightsSurfaceMapper &SymmetricWeightsSurfaceMapper
::operator =(const SymmetricWeightsSurfaceMapper &other)
{
  if (this != &other) {
    LinearFixedBoundarySurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
SymmetricWeightsSurfaceMapper::~SymmetricWeightsSurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void SymmetricWeightsSurfaceMapper::ComputeMap()
{
  const bool use_direct_solver = (_NumberOfIterations < 0 || _NumberOfIterations == 1);

  typedef Eigen::MatrixXd             Values;
  typedef Eigen::SparseMatrix<double> Matrix;
  typedef Eigen::Triplet<double>      NZEntry;

  const int n = NumberOfFreePoints();
  const int m = NumberOfComponents();

  int i, j, r, c, l;

  Matrix A(n, n);
  Values b(n, m);
  {
    Array<NZEntry> w;
    Array<double>  w_ii(n, .0);
    double         w_ij; // = w_ji

    w.reserve(2 * _EdgeTable->NumberOfEdges() + n);
    b.setZero();

    EdgeIterator edgeIt(*_EdgeTable);
    for (edgeIt.InitTraversal(); edgeIt.GetNextEdge(i, j) != -1;) {
      r = FreePointIndex(i);
      c = FreePointIndex(j);
      if (r >= 0 || c >= 0) {
        w_ij = this->Weight(i, j);
        if (r >= 0 && c >= 0) {
          w.push_back(NZEntry(r, c, -w_ij));
          w.push_back(NZEntry(c, r, -w_ij));
        } else if (r >= 0) {
          for (l = 0; l < m; ++l) {
            b(r, l) += w_ij * GetValue(j, l);
          }
        } else if (c >= 0) {
          for (l = 0; l < m; ++l) {
            b(c, l) += w_ij * GetValue(i, l);
          }
        }
        if (r >= 0) {
          w_ii[r] += w_ij;
        }
        if (c >= 0) {
          w_ii[c] += w_ij;
        }
      }
    }
    for (r = 0; r < n; ++r) {
      w.push_back(NZEntry(r, r, w_ii[r]));
    }

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
    Eigen::ConjugateGradient<Matrix> solver(A);
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
