/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2015-2016 Imperial College London
 * Copyright 2015-2016 Andreas Schuh
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

#ifndef MIRTK_MeshlessHarmonicMap_H
#define MIRTK_MeshlessHarmonicMap_H

#include "mirtk/MeshlessMap.h"

#include "mirtk/Math.h"
#include "mirtk/PointSet.h"
#include "mirtk/Vector.h"


namespace mirtk {


/**
 * Harmonic map computed with the method of fundamental solutions (MFS)
 *
 * Li et al. (2009). Meshless harmonic volumetric mapping using fundamental solution methods.
 * IEEE Transactions on Automation Science and Engineering, 6(3), 409–422.
 */
class MeshlessHarmonicMap : public MeshlessMap
{
  mirtkObjectMacro(MeshlessHarmonicMap);

public:

  // ---------------------------------------------------------------------------
  // Auxiliaries

  /// Harmonic kernel function
  ///
  /// \param[in] d Distance of point to source point.
  ///
  /// \returns Harmonic kernel function value.
  static double H(double d);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  MeshlessHarmonicMap();

  /// Copy constructor
  MeshlessHarmonicMap(const MeshlessHarmonicMap &);

  /// Assignment operator
  MeshlessHarmonicMap &operator =(const MeshlessHarmonicMap &);

  /// Make deep copy of this volumetric map
  virtual Mapping *NewCopy() const;

  /// Destructor
  virtual ~MeshlessHarmonicMap();

  // ---------------------------------------------------------------------------
  // Evaluation

  // Import other overloads
  using MeshlessMap::Evaluate;

  /// Evaluate map at a given point
  ///
  /// \param[out] v Map value.
  /// \param[in]  x Coordinate of point along x axis at which to evaluate map.
  /// \param[in]  y Coordinate of point along y axis at which to evaluate map.
  /// \param[in]  z Coordinate of point along z axis at which to evaluate map.
  ///
  /// \returns Whether input point is inside map domain.
  virtual bool Evaluate(double *v, double x, double y, double z = .0) const;

  /// Evaluate map at a given point
  ///
  /// \param[in] x Coordinate of point along x axis at which to evaluate map.
  /// \param[in] y Coordinate of point along y axis at which to evaluate map.
  /// \param[in] z Coordinate of point along z axis at which to evaluate map.
  /// \param[in] l Index of map value component.
  ///
  /// \returns The l-th component of the map value evaluated at the given point
  ///          or the \c OutsideValue when input point is outside the map domain.
  virtual double Evaluate(double x, double y, double z = .0, int l = 0) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline double MeshlessHarmonicMap::H(double d)
{
  return .25 / (d * pi);
}


} // namespace mirtk

#endif // MIRTK_MeshlessHarmonicMap_H
