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

#ifndef MIRTK_MaximumCurvatureConstraint_H
#define MIRTK_MaximumCurvatureConstraint_H

#include "mirtk/SurfaceConstraint.h"


namespace mirtk {


/**
 * Maximum principle curvature constraint term
 */
class MaximumCurvatureConstraint : public SurfaceConstraint
{
  mirtkEnergyTermMacro(MaximumCurvatureConstraint, EM_MaximumCurvature);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Maximum curvature threshold
  mirtkPublicAttributeMacro(double, Threshold);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const MaximumCurvatureConstraint &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  MaximumCurvatureConstraint(const char * = "", double = 1.0);

  /// Copy constructor
  MaximumCurvatureConstraint(const MaximumCurvatureConstraint &);

  /// Assignment operator
  MaximumCurvatureConstraint &operator =(const MaximumCurvatureConstraint &);

  /// Destructor
  virtual ~MaximumCurvatureConstraint();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Initialize force term once input and parameters have been set
  virtual void Initialize();

  /// Update internal force data structures
  virtual void Update(bool);

protected:

  /// Evaluate energy of internal force term
  virtual double Evaluate();

  /// Evaluate internal force w.r.t. transformation parameters or surface nodes
  virtual void EvaluateGradient(double *, double, double);

};


} // namespace mirtk

#endif // MIRTK_MaximumCurvatureConstraint_H
