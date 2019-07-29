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

#include "mirtk/IntrinsicLeastAreaDistortionSurfaceMapper.h"

#include "mirtk/Math.h"
#include "mirtk/VtkMath.h"
#include "mirtk/Triangle.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/PolynomialSolvers.h"

#include "boost/math/tools/polynomial.hpp"
using namespace boost::math::tools;


namespace mirtk {


// Global flags (cf. mirtk/Options.h)
MIRTK_Common_EXPORT extern int verbose;


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void IntrinsicLeastAreaDistortionSurfaceMapper
::CopyAttributes(const IntrinsicLeastAreaDistortionSurfaceMapper &other)
{
  _MinStepLength = other._MinStepLength;
  _MaxStepLength = other._MaxStepLength;
}

// -----------------------------------------------------------------------------
IntrinsicLeastAreaDistortionSurfaceMapper
::IntrinsicLeastAreaDistortionSurfaceMapper()
:
  _MinStepLength(.1),
  _MaxStepLength(100.)
{
}

// -----------------------------------------------------------------------------
IntrinsicLeastAreaDistortionSurfaceMapper
::IntrinsicLeastAreaDistortionSurfaceMapper(
  const IntrinsicLeastAreaDistortionSurfaceMapper &other
) :
  NearOptimalIntrinsicSurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
IntrinsicLeastAreaDistortionSurfaceMapper &
IntrinsicLeastAreaDistortionSurfaceMapper
::operator =(const IntrinsicLeastAreaDistortionSurfaceMapper &other)
{
  if (this != &other) {
    NearOptimalIntrinsicSurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
IntrinsicLeastAreaDistortionSurfaceMapper
::~IntrinsicLeastAreaDistortionSurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
double IntrinsicLeastAreaDistortionSurfaceMapper
::ComputeLambda(vtkDataArray *u0, vtkDataArray *u1) const
{
  const auto digits   = cout.precision(5);
  const auto fmtflags = cout.flags(ios::fixed);

  const double eps = 1e-12; // Small number added to areas to avoid division by
                            // (close to) zero area in case of degenerate triangles

  // Determine range of coordinates used to normalize them such that area measures
  // in the original mesh and the mapped mesh have comparable order of magnitude
  const double scale = pow(Scale(u0, u1) / Scale(_Surface), 2);

  // Compute coefficients of area distortion polynomial
  // (except of the constant coefficient e which is not needed to find minimum)
  const int    energy_degree         = 8;
  const double zero[energy_degree+1] = {0.};

  vtkIdType          npts, *pts;
  double             a[3], b[3], c[3];
  polynomial<double> u_i(zero, 1), v_i(zero, 1), u_j(zero, 1), v_j(zero, 1);
  polynomial<double> area[3], ratio(zero, 2), distortion(zero, energy_degree);
  polynomial<double> energy(zero, energy_degree);

  for (int i = 0; i < 3; ++i) {
    area[i] = polynomial<double>(zero, 2);
  }
  for (vtkIdType cellId = 0; cellId < _Surface->GetNumberOfCells(); ++cellId) {
    _Surface->GetCellPoints(cellId, npts, pts);
    if (npts != 3) {
      cerr << this->NameOfType() << "::ComputeLambda: Surface mesh must be triangulated" << endl;
      exit(1);
    }

    // Get surface triangle corner points
    _Surface->GetPoint(pts[0], a);
    _Surface->GetPoint(pts[1], b);
    _Surface->GetPoint(pts[2], c);

    // Coefficients of first edge start point polynomial
    // u = (u1 - u0) lambda + u0
    u_i[0] = u0->GetComponent(pts[0], 0);
    u_i[1] = u1->GetComponent(pts[0], 0) - u0->GetComponent(pts[0], 0);
    v_i[0] = u0->GetComponent(pts[0], 1);
    v_i[1] = u1->GetComponent(pts[0], 1) - u0->GetComponent(pts[0], 1);
    for (int i = 0; i < 3; ++i) {
      // Coefficients of edge end point polynomial
      auto const&ptId = pts[(i + 1) % 3];
      u_j[0] = u0->GetComponent(ptId, 0);
      u_j[1] = u1->GetComponent(ptId, 0) - u0->GetComponent(ptId, 0);
      v_j[0] = u0->GetComponent(ptId, 1);
      v_j[1] = u1->GetComponent(ptId, 1) - u0->GetComponent(ptId, 1);
      // Polynomial of edge area contribution
      area[i] = u_i * v_j - v_i * u_j;
      // Set start point of next edge
      u_i = u_j, v_i = v_j;
    }

    // Compute parametric double area polynomial
    // and divide by twice the original area
    ratio  = area[0] + area[1] + area[2];
    ratio += eps;
    ratio *= scale / (Triangle::DoubleArea(a, b, c) + eps);

    // Calculate polynomial of distortion energy
    distortion  = ratio * ratio;
    distortion -= 1.;

    energy += distortion * distortion;
  }
  energy *= 1. / _Surface->GetNumberOfCells();

  // Find lambda value with minimum area distortion, starting with lambda=0
  // which corresponds to the map which minimizes the authalic energy
  double lambda  = .0;
  double step    = _MinStepLength + .5 * (_MaxStepLength - _MinStepLength);
  double epsilon = 1e-9;
  double value   = energy.evaluate(lambda);
  double x, dx, f, df;

  polynomial<double> derivative(zero, static_cast<unsigned int>(energy.degree() - 1));
  for (polynomial<double>::size_type d = 1; d <= energy.degree(); ++d) {
    derivative[d-1] = d * energy[d];
  }

  if (verbose > 1) {
    cout << "  Minimizing area distortion using gradient descent..." << endl;
  }
  for (int iter = 0, reject = 0; iter < 100 && reject < 5; ++iter) {
    dx = derivative.evaluate(lambda);
    x  = lambda - step * dx;
    f  = energy.evaluate(x);
    df = value - f;
    if (abs(df) < epsilon) {
      break;
    } else if (df > epsilon) {
      if (verbose > 1) {
        cout << "    Accepted step: step=" << step << ", lambda=" << x << ", value=" << f << endl;
      }
      lambda = x;
      value  = f;
      step   = min(_MaxStepLength, 1.1 * step);
      reject = 0;
    } else {
      if (verbose > 1) {
        cout << "    Rejected step: step=" << step << ", lambda=" << x << ", value=" << f << endl;
      }
      step = max(_MinStepLength, .5 * step);
      ++reject;
    }
  }
  if (verbose > 0) {
    if (verbose > 1) {
      cout << "  Minimizing area distortion using gradient descent... done\n\n";
    }
    cout << "  RMS area distortion          = " << sqrt(value) << endl;
  }

  cout.precision(digits);
  cout.flags(fmtflags);
  return lambda;
}


} // namespace mirtk
