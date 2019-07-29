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

#ifndef MIRTK_AuthalicSurfaceMapper_H
#define MIRTK_AuthalicSurfaceMapper_H

#include "mirtk/NonSymmetricWeightsSurfaceMapper.h"


namespace mirtk {


/**
 * Computes piecewise linear authalic surface map with fixed boundary
 *
 * This surface map minimizes the discrete authalic energy outlined in
 * Desbrun et al. (2002). The coefficients of the sparse linear system of
 * equations solved correspond to the generalized Barycentric coordinates
 * in Meyer et al. (2002).
 *
 * - Desbrun, Meyer, and Alliez (2002). Intrinsic parameterizations of surface meshes.
 *   Computer Graphics Forum, 21(3), 209–218.
 * - Meyer et al. (2002). Generalized Barycentric Coordinates on Irregular Polygons.
 *   Journal of Graphics Tools, 7(1), 13–22.
 */
class AuthalicSurfaceMapper : public NonSymmetricWeightsSurfaceMapper
{
  mirtkObjectMacro(AuthalicSurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Copy attributes of this class from another instance
  void CopyAttributes(const AuthalicSurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  AuthalicSurfaceMapper();

  /// Copy constructor
  AuthalicSurfaceMapper(const AuthalicSurfaceMapper &);

  /// Assignment operator
  AuthalicSurfaceMapper &operator =(const AuthalicSurfaceMapper &);

  /// Destructor
  virtual ~AuthalicSurfaceMapper();

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

#endif // MIRTK_AuthalicSurfaceMapper_H
