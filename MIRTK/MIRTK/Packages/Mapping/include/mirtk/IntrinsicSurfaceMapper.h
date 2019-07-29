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

#ifndef MIRTK_IntrinsicSurfaceMapper_H
#define MIRTK_IntrinsicSurfaceMapper_H

#include "mirtk/NonSymmetricWeightsSurfaceMapper.h"


namespace mirtk {


/**
 * Surface maps based on a class of intrinsic surface parameterizations
 *
 * Surface mapping methods based on this class use a linear combination of edge
 * weights which result from the discrete formulations of Dirichlet (angle preserving)
 * and Chi (Euler characteristic, area preserving) energy functionals (Desbrun et al., 2002)
 * obtaind using the finite element method (FEM). A discrete harmonic map, also
 * referred to as Discrete Conformal Parameterization (DCP), is obtained by
 * minimizing the Dirichlet energy functional. A convex combination map referred
 * to as Discrete Authalic Parameterization (DAP) is obtained by minimizing only
 * the Chi energy. The normalized weights are a generalization of Barycentric
 * coordinates (Meyer, 2002). The discrete Dirichlet energy functional used by
 * this mapper is the piecewise linear finite element approximation of the Laplace
 * equation (Eck et al., 1995; Wardetzky et al., 2007). These cotangent weights
 * were previously used by Pinkall and Polthier (1993) to compute a conformal surface map.
 *
 * - Wachspress (1975). A Rational Finite Element Basis.
 * - Pinkall and Polthier (1993). Computing Discrete Minimal Surfaces and Their Conjugates.
 *   Experiment. Math., 2(1), 15–36.
 * - Eck et al. (1995). Multiresolution analysis of arbitrary meshes. SIGGRAPH.
 * - Desbrun, Meyer, and Alliez (2002). Intrinsic parameterizations of surface meshes.
 *   Computer Graphics Forum, 21(3), 209–218.
 * - Meyer et al. (2002). Generalized Barycentric Coordinates on Irregular Polygons.
 *   Journal of Graphics Tools, 7(1), 13–22.
 * - Wardetzky et al. (2007). Discrete quadratic curvature energies.
 *   Computer Aided Geometric Design, 24(8–9), 499–518.
 */
class IntrinsicSurfaceMapper : public NonSymmetricWeightsSurfaceMapper
{
  mirtkObjectMacro(IntrinsicSurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Weight of conformal energy in [0, 1], weight of authalic energy is 1 - Lambda.
  mirtkPublicAttributeMacro(double, Lambda);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const IntrinsicSurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  ///
  /// \param[in] lambda Weight of conformal energy in [0, 1].
  ///                   The weight of the authalic energy is 1 - \p lambda.
  IntrinsicSurfaceMapper(double lambda = .5);

  /// Copy constructor
  IntrinsicSurfaceMapper(const IntrinsicSurfaceMapper &);

  /// Assignment operator
  IntrinsicSurfaceMapper &operator =(const IntrinsicSurfaceMapper &);

  /// Destructor
  virtual ~IntrinsicSurfaceMapper();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Weight of directed edge (i, j)
  ///
  /// \param[in] i Index of start point.
  /// \param[in] j Index of end point.
  ///
  /// \returns Weight of directed edge (i, j).
  virtual double Weight(int i, int j) const;

};


} // namespace mirtk

#endif // MIRTK_IntrinsicSurfaceMapper_H
