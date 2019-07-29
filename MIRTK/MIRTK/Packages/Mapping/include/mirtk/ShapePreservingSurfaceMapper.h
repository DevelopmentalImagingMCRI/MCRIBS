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

#ifndef MIRTK_ShapePreservingSurfaceMapper_H
#define MIRTK_ShapePreservingSurfaceMapper_H

#include "mirtk/NonSymmetricWeightsSurfaceMapper.h"


namespace mirtk {


/**
 * Linear surface mapper with Floater's shape preserving convex combination weights
 *
 * The shape preserving convex combination weights are the average Barycentric coordinates
 * of overlapping triangles containing the respective center node, where the Barycentric
 * coordinates are the unique coordinates which satisfy the convex combination interpolation
 * criteria for the case where the center node is adjacent to exactly three nodes.
 *
 * - Floater (1997). Parametrization and smooth approximation of surface triangulations.
 *   Computer Aided Geometric Design, 14(3), 231–250.
 * - Guskov (2004). An anisotropic mesh parameterization scheme. 
 *   Engineering with Computers, 20(2), 129–135.
 */
class ShapePreservingSurfaceMapper : public NonSymmetricWeightsSurfaceMapper
{
  mirtkObjectMacro(ShapePreservingSurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Copy attributes of this class from another instance
  void CopyAttributes(const ShapePreservingSurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  ShapePreservingSurfaceMapper();

  /// Copy constructor
  ShapePreservingSurfaceMapper(const ShapePreservingSurfaceMapper &);

  /// Assignment operator
  ShapePreservingSurfaceMapper &operator =(const ShapePreservingSurfaceMapper &);

  /// Destructor
  virtual ~ShapePreservingSurfaceMapper();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Weights of edges adjacent to node i
  ///
  /// \param[in]  i Index of central node.
  /// \param[in]  j Indices of adjacent nodes.
  /// \param[out] w Weights of \p d_i edges (i, j).
  /// \param[in]  d Number of adjacent nodes, i.e., node degree.
  virtual void Weights(int i, const int *j, double *w, int d) const;

};


} // namespace mirtk

#endif // MIRTK_ShapePreservingSurfaceMapper_H
