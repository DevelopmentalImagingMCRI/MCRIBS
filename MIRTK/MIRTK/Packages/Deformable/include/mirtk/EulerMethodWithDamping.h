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

#ifndef MIRTK_EulerMethodWithDamping_H
#define MIRTK_EulerMethodWithDamping_H

#include "mirtk/EulerMethod.h"


namespace mirtk {


/**
 * Minimizes deformable surface model using Euler's method with momentum
 *
 * This method has for example been used in the following paper on which this
 * implementation is also based on:
 *
 *   Park et al. (2001), A non-self-intersecting adaptive deformable surface for
 *   complex boundary extraction from volumetric images,
 *   Computer & Graphics, 25, 421–440
 *
 * Another variant of an Euler method with momentum used by the FreeSurfer tools
 * is implemented by EulerMethodWithMomentum.
 */
class EulerMethodWithDamping : public EulerMethod
{
  mirtkOptimizerMacro(EulerMethodWithDamping, OM_EulerMethodWithDamping);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Common mass of all nodes
  mirtkPublicAttributeMacro(double, BodyMass);

  /// Damping coefficient
  mirtkPublicAttributeMacro(double, DampingFactor);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const EulerMethodWithDamping &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  EulerMethodWithDamping(ObjectiveFunction * = NULL);

  /// Copy constructor
  EulerMethodWithDamping(const EulerMethodWithDamping &);

  /// Assignment operator
  EulerMethodWithDamping &operator =(const EulerMethodWithDamping &);

  /// Destructor
  virtual ~EulerMethodWithDamping();

  // ---------------------------------------------------------------------------
  // Parameters
  using LocalOptimizer::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameters as key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Execution

  /// Initialize optimization
  ///
  /// This member funtion is implicitly called by Run.
  virtual void Initialize();

  /// Update node displacements
  virtual void UpdateDisplacement();

};


} // namespace mirtk

#endif // MIRTK_EulerMethodWithDamping_H
