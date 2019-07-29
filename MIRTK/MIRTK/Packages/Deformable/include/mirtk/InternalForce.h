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

#ifndef MIRTK_InternalForce_H
#define MIRTK_InternalForce_H

#include "mirtk/PointSetForce.h"
#include "mirtk/InternalForceTerm.h"


namespace mirtk {


/**
 * Base class for a penalty term imposed on a transformed point set
 *
 * Subclasses represent in particular internal forces on deformable simplicial
 * complexes such as elasticity, strain, curvature, non-self-intersection (i.e.
 * triangle repulsion), and node repulsion.
 *
 * The penalty is minimized by the registration using an instance of
 * RegistrationEnergy with set data similarity and regularization terms.
 * Higher penalties lead to stronger enforcement of the constraint.
 *
 * The same force terms can be used in a "non-parametric" DeformableSurfaceModel
 * as internal force terms to regularize the evolution of the surface model.
 */
class InternalForce : public PointSetForce
{
  mirtkAbstractMacro(InternalForce);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Name of point data array with external force magnitude
  ///
  /// When specified, the force vectors of the internal force should be
  /// multiplied by the absolute value of the external force magnitude.
  /// The internal force is therefore proportional to the external force
  /// and vanishes at the same time as the external force.
  ///
  /// The deformable model has to ensure that the external forces are updated
  /// before the internal forces such that the magnitude array is up-to-date.
  mirtkPublicAttributeMacro(string, ExternalMagnitudeArrayName);

  /// Weight of this internal force when "inside" the object,
  /// i.e., when the external force magnitude is negative.
  /// The force _Weight is multiplied by this positive weight.
  mirtkPublicAttributeMacro(double, WeightInside);

  /// Weight of this internal force when "outside" the object,
  /// i.e., when the external force magnitude is positive.
  /// The force _Weight is multiplied by this positive weight.
  mirtkPublicAttributeMacro(double, WeightOutside);

  /// Minimum weight of this internal force when external force magnitude
  /// is applied to scale the force. It can be used to preserve some internal
  /// force even when the external force vanishes.
  mirtkPublicAttributeMacro(double, WeightMinimum);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const InternalForce &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Constructor
  InternalForce(const char * = "", double = 1.0);

  /// Copy constructor
  InternalForce(const InternalForce &);

  /// Assignment operator
  InternalForce &operator =(const InternalForce &);

public:

  /// Instantiate new constraint representing specified internal forces
  static InternalForce *New(InternalForceTerm, const char * = "", double = 1.0);

  /// Destructor
  virtual ~InternalForce();

  // ---------------------------------------------------------------------------
  // Configuration

protected:

  /// Set parameter value from string
  virtual bool SetWithoutPrefix(const char *, const char *);

public:

  // Import other overloads
  using PointSetForce::Parameter;

  /// Get parameter key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Auxiliaries

protected:

  /// Get magnitude array of external force term
  vtkDataArray *ExternalMagnitude() const;

  // ---------------------------------------------------------------------------
  // Evaluation

protected:

  /// Evaluate gradient of force term
  ///
  /// \param[in,out] gradient Gradient to which the computed gradient of the
  ///                         force term is added after multiplying by \p weight.
  /// \param[in]     step     Step length for finite differences (unused).
  /// \param[in]     weight   Weight of force term.
  virtual void EvaluateGradient(double *gradient, double step, double weight);

};


} // namespace mirtk

#endif // MIRTK_InternalForce_H
