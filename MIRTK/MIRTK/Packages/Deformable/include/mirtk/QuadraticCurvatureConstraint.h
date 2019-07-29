/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Andreas Schuh
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

#ifndef MIRTK_QuadraticCurvatureConstraint_H
#define MIRTK_QuadraticCurvatureConstraint_H

#include "mirtk/SurfaceConstraint.h"

#include "mirtk/Array.h"


namespace mirtk {


/**
 * Surface curvature constraint
 *
 * This surface constraint fits a 1D quadratic to the surface locally 
 * and moves the vertex in the normal direction to improve the fit.
 * It is based on FreeSurfer's equivalent mrisComputeQuadraticCurvatureTerm
 * function defined in the mrisurf.c open source file.
 */
class QuadraticCurvatureConstraint : public SurfaceConstraint
{
  mirtkEnergyTermMacro(QuadraticCurvatureConstraint, EM_QuadraticCurvature);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Copy attributes of this class from another instance
  void CopyAttributes(const QuadraticCurvatureConstraint &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  QuadraticCurvatureConstraint(const char * = "", double = 1.0);

  /// Copy constructor
  QuadraticCurvatureConstraint(const QuadraticCurvatureConstraint &);

  /// Assignment operator
  QuadraticCurvatureConstraint &operator =(const QuadraticCurvatureConstraint &);

  /// Destructor
  virtual ~QuadraticCurvatureConstraint();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Initialize force term once input and parameters have been set
  virtual void Initialize();

  /// Update internal force data structures
  virtual void Update(bool);

protected:

  /// Compute penalty for current transformation estimate
  virtual double Evaluate();

  /// Compute internal force w.r.t. transformation parameters
  virtual void EvaluateGradient(double *, double, double);

};


} // namespace mirtk

#endif // MIRTK_QuadraticCurvatureConstraint_H
