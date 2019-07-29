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

#ifndef MIRTK_MeanValueSurfaceMapper_H
#define MIRTK_MeanValueSurfaceMapper_H

#include "mirtk/NonSymmetricWeightsSurfaceMapper.h"


namespace mirtk {


/**
 * Piecewise linear surface map using Floater's mean value convex combination weights
 *
 * - Floater (2003). Mean value coordinates. Computer Aided Geometric Design, 20(1):19–37.
 */
class MeanValueSurfaceMapper : public NonSymmetricWeightsSurfaceMapper
{
  mirtkObjectMacro(MeanValueSurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Copy attributes of this class from another instance
  void CopyAttributes(const MeanValueSurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  MeanValueSurfaceMapper();

  /// Copy constructor
  MeanValueSurfaceMapper(const MeanValueSurfaceMapper &);

  /// Assignment operator
  MeanValueSurfaceMapper &operator =(const MeanValueSurfaceMapper &);

  /// Destructor
  virtual ~MeanValueSurfaceMapper();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Weight of directed edge (i, j)
  ///
  /// \param[in] i Index of start node.
  /// \param[in] j Index of end node.
  ///
  /// \returns Weight of directed edge (i, j).
  virtual double Weight(int i, int j) const;

};


} // namespace mirtk

#endif // MIRTK_MeanValueSurfaceMapper_H
