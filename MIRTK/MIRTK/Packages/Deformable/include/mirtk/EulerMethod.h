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

#ifndef MIRTK_EulerMethod_H
#define MIRTK_EulerMethod_H

#include "mirtk/LocalOptimizer.h"

#include "vtkSmartPointer.h"
#include "vtkDataArray.h"


namespace mirtk {


class DeformableSurfaceModel;


/**
 * Minimizes deformable surface model using Euler's method
 */
class EulerMethod : public LocalOptimizer
{
  mirtkOptimizerMacro(EulerMethod, OM_EulerMethod);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Deformable surface model
  mirtkReadOnlyAggregateMacro(DeformableSurfaceModel, Model);

  /// Length of each integration step (\f$\delta t\f$)
  ///
  /// When _NormalizeStepLength is \c true, the step length is divided by the
  /// maximum magnitude of the computed node forces at a given iteration.
  /// Otherwise it is divided by the number of deformable surface points to
  /// remove the normalization of the gradient of the deformable surface terms
  /// and the resulting node forces are multplied by the unnormalized step length.
  /// Use _MaximumDisplacement in this case to limit the maximum allowed
  /// displacement of each node in this case.
  mirtkPublicAttributeMacro(double, StepLength);

  /// Whether to normalize the step length using the magnitude of the maximum force
  mirtkPublicAttributeMacro(bool, NormalizeStepLength);

  /// Maximum node displacement
  ///
  /// Maximum allowed magnitude of a node displacement, i.e., node force times
  /// current time step. When _NormalizeStepLength is \c true (default), the
  /// maximum displacement is naturally limited to 1 [mm]. When set to a value
  /// less than one in this case, the node displacements are further restricted.
  /// This parameter is in particular useful to limit the maximum node displacement
  /// when _NormalizeStepLength is \c false. The default value is 1 [mm].
  mirtkPublicAttributeMacro(double, MaximumDisplacement);

  /// Point data array of current node displacement
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkDataArray>, Displacement);

  /// Point data array used to track node displacement in normal direction (optional)
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkDataArray>, NormalDisplacement);

  /// Current negated force acting on each node
  mirtkReadOnlyComponentMacro(double, Gradient);

  /// Last maximum node displacement
  mirtkReadOnlyAttributeMacro(double, LastDelta);

private:

  /// Size of allocated vectors, may be larger than actual number of model DoFs!
  int _NumberOfDOFs;

  /// Copy attributes of this class from another instance
  void CopyAttributes(const EulerMethod &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  EulerMethod(ObjectiveFunction * = NULL);

  /// Copy constructor
  EulerMethod(const EulerMethod &);

  /// Assignment operator
  EulerMethod &operator =(const EulerMethod &);

  /// Destructor
  virtual ~EulerMethod();

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

  /// Integrate deformable surface model
  virtual double Run();

protected:

  /// Perform local adaptive remeshing (optional)
  virtual bool RemeshModel();

  /// Get node force normalization factor
  virtual double GradientNorm() const;

  /// Update node displacements
  virtual void UpdateDisplacement();

  /// Truncate magnitude of node displacements to the range [0, max]
  ///
  /// \param[in] force Force truncation of displacements. If \c false, the
  ///                  function does nothing if the step length was normalized
  ///                  and the maximum displacement is greater or equal one.
  virtual void TruncateDisplacement(bool force = false);

  /// Update recorded node displacement in normal direction
  virtual void UpdateNormalDisplacement();

  /// Finalize optimization
  virtual void Finalize();

};


} // namespace mirtk

#endif // MIRTK_EulerMethod_H
