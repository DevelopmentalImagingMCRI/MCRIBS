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

#ifndef MIRTK_HarmonicSurfaceMapper_H
#define MIRTK_HarmonicSurfaceMapper_H

#include "mirtk/SymmetricWeightsSurfaceMapper.h"


namespace mirtk {


/**
 * Computes piecewise linear harmonic surface map with fixed boundary
 *
 * The harmonic surface map minimizes the discrete Dirichlet energy. It is the
 * solution of the discrete Laplace equation. Moreover, because the boundary
 * points are fixed and therefore the area of the parameteric domain is constant,
 * this minimum is also a critical point of the conformal energy (Mullen et al., 2008).
 * But because of these Dirichlet boundary conditions, the resulting map can only
 * be as conformal as possible in a least squares sense. For a free boundary
 * surface map which minimizes the conformal energy, see the NaturalConformalSurfaceMapper.
 *
 * The surface map is the solution of a sparse linear system of equations, where
 * the coefficients are the well-known cotangent weights of the discretized
 * Laplace-Beltrami operator.
 *
 * - Pinkall and Polthier (1993). Computing Discrete Minimal Surfaces and Their Conjugates.
 *   Experiment. Math., 2(1), 15–36.
 * - Eck et al. (1995). Multiresolution analysis of arbitrary meshes. SIGGRAPH.
 * - Meyer et al. (2002). Generalized Barycentric Coordinates on Irregular Polygons.
 *   Journal of Graphics Tools, 7(1), 13–22.
 * - Desbrun, Meyer, and Alliez (2002). Intrinsic parameterizations of surface meshes.
 *   Computer Graphics Forum, 21(3), 209–218.
 * - Mullen et al. (2008). Spectral conformal parameterization.
 *   Eurographics Symposium on Geometry Processing, 27(5), 1487–1494.
 */
class HarmonicSurfaceMapper : public SymmetricWeightsSurfaceMapper
{
  mirtkObjectMacro(HarmonicSurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Copy attributes of this class from another instance
  void CopyAttributes(const HarmonicSurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  HarmonicSurfaceMapper();

  /// Copy constructor
  HarmonicSurfaceMapper(const HarmonicSurfaceMapper &);

  /// Assignment operator
  HarmonicSurfaceMapper &operator =(const HarmonicSurfaceMapper &);

public:

  /// Destructor
  virtual ~HarmonicSurfaceMapper();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Weight of undirected edge (i, j)
  ///
  /// \param[in] i First end point.
  /// \param[in] j Second end point.
  ///
  /// \returns Weight of undirected edge (i, j).
  virtual double Weight(int i, int j) const;

};


} // namespace mirtk

#endif // MIRTK_HarmonicSurfaceMapper_H
