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

#ifndef MIRTK_ChordLengthBoundarySegmentParameterizer_H
#define MIRTK_ChordLengthBoundarySegmentParameterizer_H

#include "mirtk/BoundarySegmentParameterizer.h"


namespace mirtk {


/**
 * Boundary curve parameterization which preserves relative point distances
 *
 * This class parameterizes the boundary segment as a curve with uniform
 * speed, i.e., without local distortion of the boundary segment. The curve
 * parameter value is linearly proportional to the distance of the boundary
 * points from the point at t=0. The optionally selected points are used to
 * define the point with parameter value t=0 and in which direction the t
 * values are increasing.
 *
 * This boundary curve parameterization is referred to as chord length
 * parameterization in Floater (1997).
 *
 * - Floater (1997). Parametrization and smooth approximation of surface triangulations.
 *   Computer Aided Geometric Design, 14(3), 231–250.
 */
class ChordLengthBoundarySegmentParameterizer : public BoundarySegmentParameterizer
{
  mirtkObjectMacro(ChordLengthBoundarySegmentParameterizer);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Copy attributes of this class from another instance
  void CopyAttributes(const ChordLengthBoundarySegmentParameterizer &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  ChordLengthBoundarySegmentParameterizer();

  /// Copy constructor
  ChordLengthBoundarySegmentParameterizer(const ChordLengthBoundarySegmentParameterizer &);

  /// Assignment operator
  ChordLengthBoundarySegmentParameterizer &operator =(const ChordLengthBoundarySegmentParameterizer &);

  /// Destructor
  virtual ~ChordLengthBoundarySegmentParameterizer();

  /// New copy of this parameterizer
  virtual BoundarySegmentParameterizer *NewCopy() const;

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Parameterize boundary curve
  virtual void Parameterize();

};


} // namespace mirtk

#endif // MIRTK_ChordLengthBoundarySegmentParameterizer_H
