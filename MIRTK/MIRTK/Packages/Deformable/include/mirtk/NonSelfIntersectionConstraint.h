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

#ifndef MIRTK_NonSelfIntersectionConstraint_H
#define MIRTK_NonSelfIntersectionConstraint_H

#include "mirtk/SurfaceConstraint.h"
#include "mirtk/SurfaceCollisions.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkDataArray.h"


namespace mirtk {


/**
 * Non-self intersecting surface forces
 *
 * This internal force makes nearby non-adjacent faces to repel each other
 * in order to avoid self-intersection of the surface.
 *
 * Park et al., A non-self-intersecting adaptive deformable surface for
 * complex boundary extraction from volumetric images, 25, 421–440 (2001).
 */
class NonSelfIntersectionConstraint : public SurfaceConstraint
{
  mirtkEnergyTermMacro(NonSelfIntersectionConstraint, EM_NonSelfIntersection);

  // ---------------------------------------------------------------------------
  // Attributes

protected:

  /// Use fast approximate surface triangle collision test
  mirtkPublicAttributeMacro(bool, FastCollisionTest);

  /// Minimum distance
  /// Set to non-positive value to use average edge length.
  mirtkPublicAttributeMacro(double, MinDistance);

  /// Maximum angle between face normal and center to closest point vector
  /// required for collision to be detected
  mirtkPublicAttributeMacro(double, MaxAngle);

  /// Detected surface collisions
  mirtkAttributeMacro(SurfaceCollisions::CollisionsArray, Collisions);

  /// Number of found surface collisions
  mirtkReadOnlyAttributeMacro(int, NumberOfCollisions);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const NonSelfIntersectionConstraint &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  NonSelfIntersectionConstraint(const char * = "", double = 1.0);

  /// Copy constructor
  NonSelfIntersectionConstraint(const NonSelfIntersectionConstraint &);

  /// Assignment operator
  NonSelfIntersectionConstraint &operator =(const NonSelfIntersectionConstraint &);

  /// Destructor
  virtual ~NonSelfIntersectionConstraint();

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

  /// Initialize internal force term once input and parameters have been set
  virtual void Initialize();

  /// Reinitialize internal force term after change of input topology
  virtual void Reinitialize();

  /// Update internal state upon change of input
  virtual void Update(bool = true);

protected:

  /// Common (re-)initialization steps of this class (non-virtual function!)
  void Init();

  /// Compute penalty for current transformation estimate
  virtual double Evaluate();

  /// Compute internal force w.r.t. transformation parameters
  virtual void EvaluateGradient(double *, double, double);

  // ---------------------------------------------------------------------------
  // Debugging
public:

  /// Write input of force term
  virtual void WriteDataSets(const char *, const char *, bool = true) const;

};


} // namespace mirtk

#endif // MIRTK_NonSelfIntersectionConstraint_H
