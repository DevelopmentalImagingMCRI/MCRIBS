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

#ifndef MIRTK_BoundarySegmentMapper_H
#define MIRTK_BoundarySegmentMapper_H

#include "mirtk/BoundaryMapper.h"

#include "mirtk/Memory.h"
#include "mirtk/BoundarySegmentParameterizer.h"


namespace mirtk {


/**
 * Base class of filters which assign map values to one surface boundary segment
 */
class BoundarySegmentMapper : public BoundaryMapper
{
  mirtkAbstractMacro(BoundarySegmentMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Closed boundary segment parameterizer
  mirtkPublicAttributeMacro(SharedPtr<BoundarySegmentParameterizer>, Parameterizer);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const BoundarySegmentMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  BoundarySegmentMapper();

  /// Copy constructor
  BoundarySegmentMapper(const BoundarySegmentMapper &);

  /// Assignment operator
  BoundarySegmentMapper &operator =(const BoundarySegmentMapper &);

public:

  /// Destructor
  virtual ~BoundarySegmentMapper();

  // ---------------------------------------------------------------------------
  // Execution

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Process longest boundary segment
  void MapLongest();

  /// Process boundary segment with the most number of points
  void MapLargest();

  /// Process specified boundary segment
  ///
  /// \param[in] n Index of boundary segment.
  void MapSegment(int n);

  /// Finalize boundary map
  virtual void Finalize();

protected:

  /// Assign map values to longest boundary of input surface
  ///
  /// When the surface has more than one boundary segment and a different boundary
  /// segment should be processed or when multiple boundary segments should be
  /// processed with different boundary map settings, perform the following steps
  /// instead of calling the Run() function.
  ///
  /// 1. Initialize() filter
  /// 2. Repeat the following steps for each boundary segment to be mapped:
  ///    a) Set selected boundary segment points and map parameters.
  ///    b) MapSegment() to assign map values to current segment.
  /// 3. Finalize() boundary map
  virtual void ComputeMap();

  /// Map boundary segment
  ///
  /// \param[in] n         Index of boundary segment.
  /// \param[in] indices   Indices of boundary points forming a closed line strip
  ///                      that discretizes the current surface boundary segment.
  /// \param[in] tvalues   Curve parameter in [0, 1) for each boundary segment point,
  ///                      where the first point has value t=0 and the parameter value
  ///                      for consecutive points is proportional to the distance of
  ///                      the point from the first point along the boundary curve.
  /// \param[in] selection Indices in \p i and \p t arrays corresponding to
  ///                      selected boundary points.
  virtual void MapBoundarySegment(int n, const Array<int>    &indices,
                                         const Array<double> &tvalues,
                                         const Array<int>    &selection) = 0;

};


} // namespace mirtk

#endif // MIRTK_BoundarySegmentMapper_H
