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

#ifndef MIRTK_NormalForce_H
#define MIRTK_NormalForce_H

#include "mirtk/InternalForce.h"


namespace mirtk {


/**
 * Constant internal force in normal direction
 */
class NormalForce : public InternalForce
{
  mirtkEnergyTermMacro(NormalForce, EM_NormalForce);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  NormalForce(const char * = "", double = 1.0);

  /// Copy constructor
  NormalForce(const NormalForce &);

  /// Assignment operator
  NormalForce &operator =(const NormalForce &);

  /// Destructor
  virtual ~NormalForce();

  // ---------------------------------------------------------------------------
  // Evaluation

protected:

  /// Evaluate external force term
  virtual double Evaluate();

  /// Compute internal force w.r.t. transformation parameters
  virtual void EvaluateGradient(double *, double, double);

};


} // namespace mirtk

#endif // MIRTK_NormalForce_H
