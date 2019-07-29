/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
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

#ifndef MIRTK_RepulsiveForce_H
#define MIRTK_RepulsiveForce_H

#include "mirtk/SurfaceConstraint.h"

#include "vtkSmartPointer.h"
#include "vtkAbstractPointLocator.h"


namespace mirtk {


/**
 * Internal force which makes non-neighboring nodes repell each other
 *
 * This force term is based on the respective mrisComputeRepulsiveTerm
 * found in mrisurf.c of the FreeSurfer open source implementation.
 */
class RepulsiveForce : public SurfaceConstraint
{
  mirtkEnergyTermMacro(RepulsiveForce, EM_RepulsiveForce);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Radius within which repelling force is active
  ///
  /// Marks the intersection point of the quadratic distance weighting function
  /// with the x axis. The intersection point with the y axis is equal to the
  /// weight of this force term.
  mirtkPublicAttributeMacro(double, FrontfaceRadius);

  /// Radius within which repelling force is active
  ///
  /// Marks the intersection point of the quadratic distance weighting function
  /// with the x axis. The intersection point with the y axis is equal to the
  /// weight of this force term.
  mirtkPublicAttributeMacro(double, BackfaceRadius);

  /// Point locator
  mirtkAttributeMacro(vtkSmartPointer<vtkAbstractPointLocator>, Locator);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const RepulsiveForce &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  RepulsiveForce(const char * = "", double = 1.0);

  /// Copy constructor
  RepulsiveForce(const RepulsiveForce &);

  /// Assignment operator
  RepulsiveForce &operator =(const RepulsiveForce &);

  /// Destructor
  virtual ~RepulsiveForce();

  /// Initialize force term once input and parameters have been set
  virtual void Initialize();

  /// Reinitialize force term after change of input topology
  ///
  /// This function is called in particular when an input surface has been
  /// reparameterized, e.g., by a local remeshing filter.
  virtual void Reinitialize();

  // ---------------------------------------------------------------------------
  // Configuration

protected:

  /// Set parameter value from string
  virtual bool SetWithoutPrefix(const char *, const char *);

public:

  // Import other overloads
  using SurfaceConstraint::Parameter;

  /// Get parameter name/value pairs
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Update moving input points and internal state of force term
  virtual void Update(bool = true);

protected:

  /// Evaluate energy of internal force term
  virtual double Evaluate();

  /// Evaluate internal force w.r.t. transformation parameters or surface nodes
  virtual void EvaluateGradient(double *, double, double);

};


} // namespace mirtk

#endif // MIRTK_RepulsiveForce_H
