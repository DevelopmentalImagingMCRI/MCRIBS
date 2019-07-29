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

#ifndef MIRTK_BoundarySegmentParameterizer_H
#define MIRTK_BoundarySegmentParameterizer_H

#include "mirtk/Object.h"

#include "mirtk/Memory.h"
#include "mirtk/Array.h"

#include "mirtk/BoundarySegment.h"


namespace mirtk {


/**
 * Base class of filters which parameterize a closed boundary segment
 *
 * A boundary segment parameterizer is used to assign each point of a closed
 * surface boundary segment a value t in [0, 1), with strictly increasing
 * (or decreasing) values along the boundary curve. This parameterization is
 * used in particular by a boundary mapper to assign fixed map values to the
 * boundary points.
 *
 * For example, the BoundaryToSquareMapper maps points with a t value in
 * [0, .25) to the first side of the square, points with a t value in [.25, .5)
 * to the second side, and so on. Note that the t value is proportional to the
 * distance of the point from the curve point at t=0 along the mapped curve.
 * Note further that the boundary segment point with index 0 need not be the
 * point with parameter value t=0. It is also not required that there exists
 * even a discrete curve point with parameter t=0.
 *
 * Specialized boundary parameterizers may assign parameter values to the
 * boundary segment of a given surface mesh not only based on this surface
 * and possibly related imaging data, but also other surface meshes with
 * known or to be determined correspondences between surface boundaries.
 * Such parameterizer will assign equal t values to corresponding boundary
 * points such that these points are mapped to the same value by the boundary
 * mapper. Other parameterizers may take an ordered/labeled list of manually
 * selected boundary points as input and assign parameter values based on
 * this selection.
 *
 * In a nutshell, a boundary parameterizer indirectly establishes
 * correspondences between boundary points of 2D manifolds, while the
 * boundary mapper itself only defines the shape of the parametric domain,
 * i.e., the codomain of the boundary map.
 *
 * Subclasses may further consider the ordered (not sorted!) list of user
 * selected boundary segment points (see BoundarySegment::IsSelected).
 * 
 * How many points need to be selected depends on the specific parameterizer.
 * In general, no points have to be selected, but selecting boundary points
 * can be used to orient the curve, i.e., to decide in which direction along
 * the closed boundary the t values are increasing. The curve orientation is
 * generally defined by the order of the selected points. The first point in
 * the selection usually defines the boundary point with parameter value t=0.
 *
 * \sa BoundaryMapper, BoundarySegment
 */
class BoundarySegmentParameterizer : public Object
{
  mirtkAbstractMacro(BoundarySegmentParameterizer);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Surface boundary segment
  mirtkPublicAttributeMacro(BoundarySegment, Boundary);

  /// Parameter value at boundary segment points
  mirtkReadOnlyAttributeMacro(Array<double>, Values);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const BoundarySegmentParameterizer &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  BoundarySegmentParameterizer();

  /// Copy constructor
  BoundarySegmentParameterizer(const BoundarySegmentParameterizer &);

  /// Assignment operator
  BoundarySegmentParameterizer &operator =(const BoundarySegmentParameterizer &);

public:

  /// Destructor
  virtual ~BoundarySegmentParameterizer();

  /// New copy of this parameterizer
  virtual BoundarySegmentParameterizer *NewCopy() const = 0;

  // ---------------------------------------------------------------------------
  // Execution

public:

  /// Parameterize boundary segment
  virtual void Run();

  /// Get parameter value of i-th boundary segment point
  ///
  /// \param[in] i Index of boundary segment point.
  ///
  /// \returns Parameter value of i-th boundary segment point
  double Value(int i) const;

protected:

  /// Initialize parameterizer after input and parameters are set
  virtual void Initialize();

  /// Parameterize boundary curve
  virtual void Parameterize() = 0;

  /// Finalize parameterization
  virtual void Finalize();

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline double BoundarySegmentParameterizer::Value(int i) const
{
  return _Values[i];
}


} // namespace mirtk

#endif // MIRTK_BoundarySegmentParameterizer_H
