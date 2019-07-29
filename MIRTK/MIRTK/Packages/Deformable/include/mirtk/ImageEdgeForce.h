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

#ifndef MIRKT_ImageEdgeForce_H
#define MIRKT_ImageEdgeForce_H

#include "mirtk/ExternalForce.h"

#include "mirtk/GenericImage.h"


namespace mirtk {


/**
 * External surface force which attracts the surface to image edges
 */
class ImageEdgeForce : public ExternalForce
{
  mirtkEnergyTermMacro(ImageEdgeForce, EM_ImageEdgeForce);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Standard deviation of Gaussian smoothing kernel
  mirtkPublicAttributeMacro(double, Sigma);

  /// Whether to project edge field gradient onto face normal
  mirtkPublicAttributeMacro(bool, InNormalDirection);

  /// Edge field
  mirtkAttributeMacro(RealImage, EdgeField);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Copy attributes of this class from another instance
  void CopyAttributes(const ImageEdgeForce &);

public:

  /// Constructor
  ImageEdgeForce(const char * = "", double = 1.0);

  /// Copy constructor
  ImageEdgeForce(const ImageEdgeForce &);

  /// Assignment operator
  ImageEdgeForce &operator =(const ImageEdgeForce &);

  /// Destructor
  virtual ~ImageEdgeForce();

  // ---------------------------------------------------------------------------
  // Configuration

protected:

  /// Set parameter value from string
  virtual bool SetWithoutPrefix(const char *, const char *);

public:

  // Import other overloads
  using ExternalForce::Parameter;

  /// Get parameter name/value pairs
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Initialization

  /// Initialize external force once input and parameters have been set
  virtual void Initialize();

  // ---------------------------------------------------------------------------
  // Evaluation

protected:

  /// Evaluate external force
  virtual void EvaluateGradient(double *, double, double);

};


} // namespace mirtk

#endif // MIRKT_ImageEdgeForce_H
