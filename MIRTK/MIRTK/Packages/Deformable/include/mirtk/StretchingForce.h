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

#ifndef MIRTK_StretchingForce_H
#define MIRTK_StretchingForce_H

#include "mirtk/InternalForce.h"


namespace mirtk {


/**
 * A spring force which favors a given rest edge length
 *
 * This spring force keeps the edge lengths of a simplicial complex either
 * similar to a specified average edge length or a predefined rest edge lenght.
 * It is usually used in conjunction a with local remeshing filter.
 *
 * Lachaud and Montanvert, Deformable meshes with automated topology changes
 * for coarse-to-fine three-dimensional surface extraction. Medical Image Analysis,
 * 3(2), 187–207, doi:10.1016/S1361-8415(99)80012-7 (1999).
 *
 * Park et al., A non-self-intersecting adaptive deformable surface for
 * complex boundary extraction from volumetric images, 25, 421–440 (2001).
 */
class StretchingForce : public InternalForce
{
  mirtkEnergyTermMacro(StretchingForce, EM_Stretching);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Rest edge length
  mirtkPublicAttributeMacro(double, RestLength);

  /// Average edge length, used only when _RestLength < 0
  mirtkAttributeMacro(double, AverageLength);

  /// Whether to use current average edge length at each iteration
  mirtkPublicAttributeMacro(bool, UseCurrentAverageLength);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const StretchingForce &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  StretchingForce(const char * = "", double = 1.0);

  /// Copy constructor
  StretchingForce(const StretchingForce &);

  /// Assignment operator
  StretchingForce &operator =(const StretchingForce &);

  /// Destructor
  virtual ~StretchingForce();

  // ---------------------------------------------------------------------------
  // Configuration

protected:

  /// Set parameter value from string
  virtual bool SetWithoutPrefix(const char *, const char *);

public:

  // Import other overloads
  using InternalForce::Parameter;

  /// Get parameter name/value pairs
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Initialize internal force term once input and parameters have been set
  virtual void Initialize();

  /// Reinitialize internal force term after change of input topology
  virtual void Reinitialize();

  /// Update internal force data structures
  virtual void Update(bool);

protected:

  /// Common (re-)initialization steps of this internal force term
  /// \note Must be a non-virtual function!
  void Init();

  /// Compute penalty for current transformation estimate
  virtual double Evaluate();

  /// Compute internal force w.r.t. transformation parameters
  virtual void EvaluateGradient(double *, double, double);

};


} // namespace mirtk

#endif // MIRTK_StretchingForce_H
