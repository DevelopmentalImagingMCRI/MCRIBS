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

#ifndef MIRTK_SymmetricWeightsSurfaceMapper_H
#define MIRTK_SymmetricWeightsSurfaceMapper_H

#include "mirtk/LinearFixedBoundarySurfaceMapper.h"


namespace mirtk {


/**
 * Piecewise linear fixed boundary surface mapping with symmetric edge weights
 */
class SymmetricWeightsSurfaceMapper : public LinearFixedBoundarySurfaceMapper
{
  mirtkAbstractMacro(SymmetricWeightsSurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Copy attributes of this class from another instance
  void CopyAttributes(const SymmetricWeightsSurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  SymmetricWeightsSurfaceMapper();

  /// Copy constructor
  SymmetricWeightsSurfaceMapper(const SymmetricWeightsSurfaceMapper &);

  /// Assignment operator
  SymmetricWeightsSurfaceMapper &operator =(const SymmetricWeightsSurfaceMapper &);

public:

  /// Destructor
  virtual ~SymmetricWeightsSurfaceMapper();

  // ---------------------------------------------------------------------------
  // Execution

  /// Construct and solve symmetric system of linear equations
  virtual void ComputeMap();

protected:

  /// Weight of undirected edge (i, j)
  ///
  /// \param[in] i First end point.
  /// \param[in] j Second end point.
  ///
  /// \returns Weight of undirected edge (i, j).
  virtual double Weight(int i, int j) const = 0;

};


} // namespace mirtk

#endif // MIRTK_SymmetricWeightsSurfaceMapper_H
