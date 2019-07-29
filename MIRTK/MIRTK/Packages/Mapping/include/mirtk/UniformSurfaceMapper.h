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

#ifndef MIRTK_UniformSurfaceMapper_H
#define MIRTK_UniformSurfaceMapper_H

#include "mirtk/SymmetricWeightsSurfaceMapper.h"


namespace mirtk {


/**
 * Piecewise linear surface mapper with uniform convex combination weights
 *
 * - Tutte (1963). How to draw a graph. Proc. London Math. Soc, 8(May 1962), 743–767.
 * - Taubin (1995). A Signal Processing Approach To Fair Surface Design. SIGGRAPH, 351–358.
 * - Floater (1997). Parametrization and smooth approximation of surface triangulations.
 *   Computer Aided Geometric Design, 14(3), 231–250.
 */
class UniformSurfaceMapper : public SymmetricWeightsSurfaceMapper
{
  mirtkObjectMacro(UniformSurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Copy attributes of this class from another instance
  void CopyAttributes(const UniformSurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  UniformSurfaceMapper();

  /// Copy constructor
  UniformSurfaceMapper(const UniformSurfaceMapper &);

  /// Assignment operator
  UniformSurfaceMapper &operator =(const UniformSurfaceMapper &);

  /// Destructor
  virtual ~UniformSurfaceMapper();

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

#endif // MIRTK_UniformSurfaceMapper_H
