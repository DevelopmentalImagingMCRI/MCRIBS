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

#ifndef MIRTK_SubdividedBoundarySegmentParameterizer_H
#define MIRTK_SubdividedBoundarySegmentParameterizer_H

#include "mirtk/BoundarySegmentParameterizer.h"


namespace mirtk {


/**
 * Hybrid of uniform and chord length parameterization
 *
 * This boundary curve parameterization is equivalent to the chord length
 * parameterization when less than two boundary segment points are selected.
 * When all points are selected, the resulting parameterization corresponds
 * to the uniform parameterization, instead. Otherwise, when more than one
 * point is selected, each curve segment between selected points is mapped
 * to the same fraction of the total boundary curve length, with chord length
 * parameterization within each sub-segment.
 *
 * \sa UniformBoundarySegmentParameterizer, ChordLengthBoundarySegmentParameterizer
 */
class SubdividedBoundarySegmentParameterizer : public BoundarySegmentParameterizer
{
  mirtkObjectMacro(SubdividedBoundarySegmentParameterizer);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Copy attributes of this class from another instance
  void CopyAttributes(const SubdividedBoundarySegmentParameterizer &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  SubdividedBoundarySegmentParameterizer();

  /// Copy constructor
  SubdividedBoundarySegmentParameterizer(const SubdividedBoundarySegmentParameterizer &);

  /// Assignment operator
  SubdividedBoundarySegmentParameterizer &operator =(const SubdividedBoundarySegmentParameterizer &);

  /// Destructor
  virtual ~SubdividedBoundarySegmentParameterizer();

  /// New copy of this parameterizer
  virtual BoundarySegmentParameterizer *NewCopy() const;

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Parameterize boundary curve
  virtual void Parameterize();

};


} // namespace mirtk

#endif // MIRTK_SubdividedBoundarySegmentParameterizer_H
