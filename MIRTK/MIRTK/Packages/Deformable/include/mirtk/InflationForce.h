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

#ifndef MIRTK_InflationForce_H
#define MIRTK_InflationForce_H

#include "mirtk/SurfaceConstraint.h"

#include "vtkSmartPointer.h"
#include "vtkDataArray.h"


namespace mirtk {


/**
 * Surface inflation force
 *
 * The inflation is driven by the average convexity or concavity of a region.
 * That is, points which lie in concave regions move outwards over time,
 * while points in convex regions move inwards.
 *
 *   Fischl et al., Cortical Surface-Based Analysis II: Inflation, Flattening,
 *   and a Surface-Based Coordinate System. NeuroImage, 9(2), 195–207 (1999).
 *
 * This force term is similar to the CurvatureConstraint. Both reduce the
 * curvature of the surface mesh by moving points towards the centroid of the
 * adjacent nodes. The inflation force is defined as the sum of the average
 * squared adjacent edge lengths, while the curvature constraint is defined as
 * the squared distance to the explicitly computed centroid position.
 *
 * Moreover, the average normal component is subtracted from the node force.
 * This modification is not presented in the NeuroImage paper by Fischl et al,
 * but can be seen in the implementation of FreeSurfer's mris_inflate tool
 * (cf. mrisComputeNormalizedSpringTerm in mrisurf.c).
 *
 * @sa CurvatureConstraint
 */
class InflationForce : public SurfaceConstraint
{
  mirtkEnergyTermMacro(InflationForce, EM_InflationForce);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  InflationForce(const char * = "", double = 1.0);

  /// Copy constructor
  InflationForce(const InflationForce &);

  /// Assignment operator
  InflationForce &operator =(const InflationForce &);

  /// Destructor
  virtual ~InflationForce();

  // ---------------------------------------------------------------------------
  // Evaluation

protected:

  /// Compute penalty for current transformation estimate
  virtual double Evaluate();

  /// Compute internal force w.r.t. transformation parameters
  virtual void EvaluateGradient(double *, double, double);

};


} // namespace mirtk

#endif // MIRTK_InflationForce_H
