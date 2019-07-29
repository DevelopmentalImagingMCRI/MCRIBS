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

#ifndef MIRTK_ChordLengthSurfaceMapper_H
#define MIRTK_ChordLengthSurfaceMapper_H

#include "mirtk/SymmetricWeightsSurfaceMapper.h"


namespace mirtk {


/**
 * Piecewise linear surface mapper with convex combination weights inverse proportional to edge length
 *
 * - Floater (1997). Parametrization and smooth approximation of surface triangulations.
 *   Computer Aided Geometric Design, 14(3), 231–250.
 */
class ChordLengthSurfaceMapper : public SymmetricWeightsSurfaceMapper
{
  mirtkObjectMacro(ChordLengthSurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Number of times the edge length is exponentiated
  mirtkPublicAttributeMacro(int, Exponent);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const ChordLengthSurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  ///
  /// \param[in] Edge length exponent.
  ChordLengthSurfaceMapper(int p = 2);

  /// Copy constructor
  ChordLengthSurfaceMapper(const ChordLengthSurfaceMapper &);

  /// Assignment operator
  ChordLengthSurfaceMapper &operator =(const ChordLengthSurfaceMapper &);

  /// Destructor
  virtual ~ChordLengthSurfaceMapper();

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

#endif // MIRTK_ChordLengthSurfaceMapper_H
