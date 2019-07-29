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

#ifndef MIRTK_InternalForceTerm_H
#define MIRTK_InternalForceTerm_H

#include "mirtk/EnergyMeasure.h"


namespace mirtk {


// -----------------------------------------------------------------------------
enum InternalForceTerm
{
  IFT_Unknown             = IFT_Begin,
  IFT_Distortion          = EM_MetricDistortion,
  IFT_Stretching          = EM_Stretching,
  IFT_Curvature           = EM_Curvature,
  IFT_QuadraticCurvature  = EM_QuadraticCurvature,
  IFT_GaussCurvature      = EM_GaussCurvature,
  IFT_MeanCurvature       = EM_MeanCurvature,
  IFT_MaximumCurvature    = EM_MaximumCurvature,
  IFT_NonSelfIntersection = EM_NonSelfIntersection,
  IFT_Repulsion           = EM_RepulsiveForce,
  IFT_Inflation           = EM_InflationForce,
  IFT_Spring              = EM_SpringForce,
  IFT_Normal              = EM_NormalForce
};

// -----------------------------------------------------------------------------
template <>
inline string ToString(const InternalForceTerm &ift, int w, char c, bool left)
{
  EnergyMeasure em = static_cast<EnergyMeasure>(ift);
  if (em <= IFT_Begin || em >= IFT_End) return ToString("Unknown", w, c, left);
  return ToString(em, w, c, left);
}

// -----------------------------------------------------------------------------
template <>
inline bool FromString(const char *str, InternalForceTerm &ift)
{
  EnergyMeasure em = EM_Unknown;
  if (FromString(str, em) && IFT_Begin < em && em < IFT_End) {
    ift = static_cast<InternalForceTerm>(em);
    return true;
  } else {
    ift = IFT_Unknown;
    return false;
  }
}


} // namespace mirtk

#endif // MIRTK_InternalForceTerm_H
