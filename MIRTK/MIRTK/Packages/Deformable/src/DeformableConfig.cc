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

#include "mirtk/DeformableConfig.h"
#include "mirtk/ObjectFactory.h"

#ifndef MIRTK_AUTO_REGISTER
  // Optimizers
  #include "mirtk/EulerMethod.h"
  #include "mirtk/EulerMethodWithDamping.h"
  #include "mirtk/EulerMethodWithMomentum.h"
  // External forces
  #include "mirtk/BalloonForce.h"
  #include "mirtk/ImageEdgeForce.h"
  #include "mirtk/ImageEdgeDistance.h"
  #include "mirtk/ImplicitSurfaceDistance.h"
  // Internal forces
  #include "mirtk/CurvatureConstraint.h"
  #include "mirtk/InflationForce.h"
  #include "mirtk/NormalForce.h"
  #include "mirtk/MetricDistortion.h"
  #include "mirtk/NonSelfIntersectionConstraint.h"
  #include "mirtk/QuadraticCurvatureConstraint.h"
  #include "mirtk/GaussCurvatureConstraint.h"
  #include "mirtk/MeanCurvatureConstraint.h"
  #include "mirtk/MaximumCurvatureConstraint.h"
  #include "mirtk/RepulsiveForce.h"
  #include "mirtk/SpringForce.h"
  #include "mirtk/StretchingForce.h"
#endif


namespace mirtk {


// -----------------------------------------------------------------------------
static void RegisterOptimizers()
{
  #ifndef MIRTK_AUTO_REGISTER
    mirtkRegisterOptimizerMacro(EulerMethod);
    mirtkRegisterOptimizerMacro(EulerMethodWithDamping);
    mirtkRegisterOptimizerMacro(EulerMethodWithMomentum);
  #endif
}

// -----------------------------------------------------------------------------
static void RegisterExternalForces()
{
  #ifndef MIRTK_AUTO_REGISTER
    mirtkRegisterEnergyTermMacro(BalloonForce);
    mirtkRegisterEnergyTermMacro(ImageEdgeForce);
    mirtkRegisterEnergyTermMacro(ImageEdgeDistance);
    mirtkRegisterEnergyTermMacro(ImplicitSurfaceDistance);
  #endif
}

// -----------------------------------------------------------------------------
static void RegisterInternalForces()
{
  #ifndef MIRTK_AUTO_REGISTER
    mirtkRegisterEnergyTermMacro(CurvatureConstraint);
    mirtkRegisterEnergyTermMacro(InflationForce);
    mirtkRegisterEnergyTermMacro(MetricDistortion);
    mirtkRegisterEnergyTermMacro(NonSelfIntersectionConstraint);
    mirtkRegisterEnergyTermMacro(QuadraticCurvatureConstraint);
    mirtkRegisterEnergyTermMacro(GaussCurvatureConstraint);
    mirtkRegisterEnergyTermMacro(MeanCurvatureConstraint);
    mirtkRegisterEnergyTermMacro(MaximumCurvatureConstraint);
    mirtkRegisterEnergyTermMacro(RepulsiveForce);
    mirtkRegisterEnergyTermMacro(SpringForce);
    mirtkRegisterEnergyTermMacro(StretchingForce);
    mirtkRegisterEnergyTermMacro(NormalForce);
  #endif
}

// -----------------------------------------------------------------------------
void InitializeDeformableLibrary()
{
  static bool initialized = false;
  if (!initialized) {
    RegisterOptimizers();
    RegisterExternalForces();
    RegisterInternalForces();
    initialized = true;
  }
}


} // namespace mirtk
