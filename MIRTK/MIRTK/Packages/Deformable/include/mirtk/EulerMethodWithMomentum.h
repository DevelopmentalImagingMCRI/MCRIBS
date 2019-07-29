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

#ifndef MIRTK_EulerMethodWithMomentum_H
#define MIRTK_EulerMethodWithMomentum_H

#include "mirtk/EulerMethod.h"


namespace mirtk {


/**
 * Minimizes deformable surface model using Euler's method with momentum
 *
 * This method is used by the FreeSurfer tools:
 *
 *   Fischl et al., Cortical Surface-Based Analysis II: Inflation, Flattening,
 *   and a Surface-Based Coordinate System. NeuroImage, 9(2), 195–207 (1999).
 *
 * Another variant of an Euler method with momentum is implemented by
 * EulerMethodWithDamping.
 */
class EulerMethodWithMomentum : public EulerMethod
{
  mirtkOptimizerMacro(EulerMethodWithMomentum, OM_EulerMethodWithMomentum);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Momentum factor
  mirtkPublicAttributeMacro(double, Momentum);

  /// Whether to exclude the momentum component from the total displacements
  /// in normal direction. By default, the _Displacement in normal direction
  /// by which a node is actually moved is integrated.
  mirtkPublicAttributeMacro(bool, ExcludeMomentumFromNormalDisplacement);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const EulerMethodWithMomentum &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  EulerMethodWithMomentum(ObjectiveFunction * = NULL);

  /// Copy constructor
  EulerMethodWithMomentum(const EulerMethodWithMomentum &);

  /// Assignment operator
  EulerMethodWithMomentum &operator =(const EulerMethodWithMomentum &);

  /// Destructor
  virtual ~EulerMethodWithMomentum();

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

protected:

  /// Update node displacements
  virtual void UpdateDisplacement();

  /// Update recorded node displacement in normal direction
  virtual void UpdateNormalDisplacement();

};


} // namespace mirtk

#endif // MIRTK_EulerMethodWithMomentum_H
