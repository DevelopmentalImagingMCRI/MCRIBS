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

#ifndef MIRTK_ExternalForce_H
#define MIRTK_ExternalForce_H

#include "mirtk/PointSetForce.h"

#include "mirtk/ExternalForceTerm.h"
#include "mirtk/RegisteredImage.h"


namespace mirtk {


/**
 * Base class for external point set force terms
 *
 * Subclasses implement in particular external forces for deformable surface
 * models such as inflation/balloon forces and intensity edge forces.
 */
class ExternalForce : public PointSetForce
{
  mirtkAbstractMacro(ExternalForce);

  // ---------------------------------------------------------------------------
  // Types
public:

  /// Type of input image voxel values
  typedef RegisteredImage::VoxelType VoxelType;

  /// Non-abstract base type of input image
  typedef GenericImage<VoxelType> ImageType;

  // ---------------------------------------------------------------------------
  // Attributes
private:

  /// (Transformed) reference image
  mirtkPublicAggregateMacro(RegisteredImage, Image);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  ExternalForce(const char * = "", double = 1.0);

  /// Copy constructor
  ExternalForce(const ExternalForce &);

  /// Assignment operator
  ExternalForce &operator =(const ExternalForce &);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const ExternalForce &);

public:

  /// Instantiate specified external force
  static ExternalForce *New(ExternalForceTerm, const char * = "", double = 1.0);

  /// Destructor
  virtual ~ExternalForce();

  // ---------------------------------------------------------------------------
  // Initialization
public:

  /// Initialize external force once input and parameters have been set
  virtual void Initialize();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Update moving input points and internal state of external force term
  virtual void Update(bool = true);

protected:

  /// Evaluate external force term
  virtual double Evaluate();

};


} // namespace mirtk

#endif // MIRTK_ExternalForce_H
