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

#ifndef MIRTK_MeanCurvatureConstraint_H
#define MIRTK_MeanCurvatureConstraint_H

#include "mirtk/SurfaceConstraint.h"

#include "vtkSmartPointer.h"
#include "vtkDataArray.h"


namespace mirtk {


/**
 * Mean curvature smoothness term
 */
class MeanCurvatureConstraint : public SurfaceConstraint
{
  mirtkEnergyTermMacro(MeanCurvatureConstraint, EM_MeanCurvature);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Maximum absolute mean curvature
  mirtkAttributeMacro(double, MaxMeanCurvature);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const MeanCurvatureConstraint &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  MeanCurvatureConstraint(const char * = "", double = 1.0);

  /// Copy constructor
  MeanCurvatureConstraint(const MeanCurvatureConstraint &);

  /// Assignment operator
  MeanCurvatureConstraint &operator =(const MeanCurvatureConstraint &);

  /// Destructor
  virtual ~MeanCurvatureConstraint();

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

#endif // MIRTK_MeanCurvatureConstraint_H
