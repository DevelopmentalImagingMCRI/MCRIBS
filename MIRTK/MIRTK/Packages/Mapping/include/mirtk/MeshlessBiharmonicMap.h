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

#ifndef MIRTK_MeshlessBiharmonicMap_H
#define MIRTK_MeshlessBiharmonicMap_H

#include "mirtk/MeshlessHarmonicMap.h"


namespace mirtk {


/**
 * Biharmonic map computed using the method of fundamental solutions (MFS)
 *
 * Xu et al. (2013). Biharmonic volumetric mapping using fundamental solutions. 
 * IEEE Transactions on Visualization and Computer Graphics, 19(5), 787–798.
 */
class MeshlessBiharmonicMap : public MeshlessHarmonicMap
{
  mirtkObjectMacro(MeshlessBiharmonicMap);

public:

  // ---------------------------------------------------------------------------
  // Auxiliaries

  /// Biharmonic kernel function
  ///
  /// \param[in] d Distance of point to source point.
  ///
  /// \returns Biharmonic kernel function value.
  static double B(double d);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  MeshlessBiharmonicMap();

  /// Copy constructor
  MeshlessBiharmonicMap(const MeshlessBiharmonicMap &);

  /// Assignment operator
  MeshlessBiharmonicMap &operator =(const MeshlessBiharmonicMap &);

  /// Initialize map after inputs and parameters are set
  virtual void Initialize();

  /// Make deep copy of this volumetric map
  virtual Mapping *NewCopy() const;

  /// Destructor
  virtual ~MeshlessBiharmonicMap();

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
  virtual bool Evaluate(double *v, double x, double y, double z = 0) const;

  /// Evaluate map at a given point
  ///
  /// \param[in] x Coordinate of point along x axis at which to evaluate map.
  /// \param[in] y Coordinate of point along y axis at which to evaluate map.
  /// \param[in] z Coordinate of point along z axis at which to evaluate map.
  /// \param[in] l Index of map value component.
  ///
  /// \returns The l-th component of the map value evaluated at the given point
  ///          or the \c OutsideValue when input point is outside the map domain.
  virtual double Evaluate(double x, double y, double z = 0, int l = 0) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline double MeshlessBiharmonicMap::B(double d)
{
  return d / (8.0 * pi);
}


} // namespace mirtk

#endif // MIRTK_MeshlessBiharmonicMap_H
