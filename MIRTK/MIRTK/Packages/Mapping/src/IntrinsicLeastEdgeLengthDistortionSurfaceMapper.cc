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

#include "mirtk/IntrinsicLeastEdgeLengthDistortionSurfaceMapper.h"

#include "mirtk/Math.h"
#include "mirtk/VtkMath.h"
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
void IntrinsicLeastEdgeLengthDistortionSurfaceMapper
::CopyAttributes(const IntrinsicLeastEdgeLengthDistortionSurfaceMapper &other)
{
}

// -----------------------------------------------------------------------------
IntrinsicLeastEdgeLengthDistortionSurfaceMapper
::IntrinsicLeastEdgeLengthDistortionSurfaceMapper()
{
}

// -----------------------------------------------------------------------------
IntrinsicLeastEdgeLengthDistortionSurfaceMapper
::IntrinsicLeastEdgeLengthDistortionSurfaceMapper(
  const IntrinsicLeastEdgeLengthDistortionSurfaceMapper &other
) :
  NearOptimalIntrinsicSurfaceMapper(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
IntrinsicLeastEdgeLengthDistortionSurfaceMapper &
IntrinsicLeastEdgeLengthDistortionSurfaceMapper
::operator =(const IntrinsicLeastEdgeLengthDistortionSurfaceMapper &other)
{
  if (this != &other) {
    NearOptimalIntrinsicSurfaceMapper::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
IntrinsicLeastEdgeLengthDistortionSurfaceMapper
::~IntrinsicLeastEdgeLengthDistortionSurfaceMapper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
double IntrinsicLeastEdgeLengthDistortionSurfaceMapper
::ComputeLambda(vtkDataArray *u0, vtkDataArray *u1) const
{
  const double eps = 1e-12; // Small number added to edge length measures to avoid
                            // division by zero in case of degenerated edge

  // Determine range of coordinates used to normalize them such that length measures
  // in the original mesh and the mapped mesh have comparable order of magnitude
  const double scale = pow(Scale(u0, u1) / Scale(_Surface), 2);

  // Compute coefficients of edge-length distortion polynomial
  // (except of the constant coefficient e which is not needed to find minimum)
  const int    energy_degree       = 4;
  const double zero[energy_degree] = {0.};

  polynomial<double> dist2(zero, 2);
  polynomial<double> distortion(zero, energy_degree);
  polynomial<double> energy(zero, energy_degree);

  int    i, j;
  double a, b, p_i[3], p_j[3], u0_i[2], u0_j[2], u1_i[2], u1_j[2];

  EdgeIterator it(*_EdgeTable);
  for (it.InitTraversal(); it.GetNextEdge(i, j) != -1;) {
    _Surface->GetPoint(i, p_i);
    _Surface->GetPoint(j, p_j);

    u0->GetTuple(i, u0_i);
    u0->GetTuple(j, u0_j);
    u1->GetTuple(i, u1_i);
    u1->GetTuple(j, u1_j);

    a = (u1_i[0] - u0_i[0]) - (u1_j[0] - u0_j[0]);
    b = (u0_i[0] - u0_j[0]);

    dist2[0] = b * b;
    dist2[1] = a * b;
    dist2[2] = a * a;

    a = (u1_i[1] - u0_i[1]) - (u1_j[1] - u0_j[1]);
    b = (u0_i[1] - u0_j[1]);

    dist2[0] += b * b;
    dist2[1] += a * b;
    dist2[2] += a * a;

    dist2[0] += eps;
    dist2[1] *= 2.;

    distortion  = dist2;
    distortion *= scale / (vtkMath::Distance2BetweenPoints(p_i, p_j) + eps);
    distortion -= 1.;

    energy += distortion * distortion;
  }

  // Find roots of derivative using cubic equation formula
  double lambda = MinimumOf4thDegreePolynomial(energy);
  if (verbose) {
    energy *= 1. / _EdgeTable->NumberOfEdges();
    cout << "  RMS edge-length distortion   = " << sqrt(energy.evaluate(lambda)) << endl;
  }
  return lambda;
}


} // namespace mirtk
