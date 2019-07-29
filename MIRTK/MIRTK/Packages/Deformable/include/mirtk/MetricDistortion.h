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

#ifndef MIRTK_MetricDistortion_H
#define MIRTK_MetricDistortion_H

#include "mirtk/SurfaceConstraint.h"

#include "mirtk/Array.h"


namespace mirtk {


/**
 * Metric distortion constraint
 *
 * Fischl et al.. Cortical Surface-Based Analysis II: Inflation, Flattening,
 * and a Surface-Based Coordinate System. NeuroImage, 9(2), 195–207 (1999).
 *
 * Additional to the formulation in the NeuroImage paper of Fischl et al.,
 * this implementation contains a slight modification found in the implementation
 * of FreeSurfer's mris_inflate tool. The normal component of the force term
 * is removed from the gradient vector.
 */
class MetricDistortion : public SurfaceConstraint
{
  mirtkEnergyTermMacro(MetricDistortion, EM_MetricDistortion);

  // ---------------------------------------------------------------------------
  // Types

public:

  struct NodeDistances
  {
    double _Distance0; ///< Initial node distance
    double _Distance;  ///< Current node distance
  };

  /// Array of computed distances to neighboring nodes
  typedef Array<NodeDistances> DistancesArray;

  // ---------------------------------------------------------------------------
  // Attributes

  /// Radius of node neighborhood (i.e., maximum edge-connectivity)
  mirtkPublicAttributeMacro(int, Radius);

  /// Area of initial surface mesh
  mirtkAttributeMacro(double, InitialArea);

  /// For each node, the initial and current distances to their neighboring nodes
  mirtkAttributeMacro(Array<DistancesArray>, Distances);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const MetricDistortion &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  MetricDistortion(const char * = "", double = 1.0);

  /// Copy constructor
  MetricDistortion(const MetricDistortion &);

  /// Assignment operator
  MetricDistortion &operator =(const MetricDistortion &);

  /// Destructor
  virtual ~MetricDistortion();

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

#endif // MIRTK_MetricDistortion_H
