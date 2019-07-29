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

#ifndef MIRTK_ExternalForceTerm_H
#define MIRTK_ExternalForceTerm_H

#include "mirtk/EnergyMeasure.h"


namespace mirtk {


// -----------------------------------------------------------------------------
/// Enumeration of available external point set force terms
///
/// \note This enumeration constains only a subset of all EnergyMeasure
///       enumeration values, whereby the integer value of corresponding
///       enumeration values is equal.
///
/// \sa EnergyMeasure
enum ExternalForceTerm
{
  EFT_Unknown                 = EFT_Begin,
  EFT_BalloonForce            = EM_BalloonForce,
  EFT_ImageEdgeForce          = EM_ImageEdgeForce,
  EFT_ImageEdgeDistance       = EM_ImageEdgeDistance,
  EFT_ImplicitSurfaceDistance = EM_ImplicitSurfaceDistance
};

// -----------------------------------------------------------------------------
template <>
inline string ToString(const ExternalForceTerm &eft, int w, char c, bool left)
{
  EnergyMeasure em = static_cast<EnergyMeasure>(eft);
  if (em <= EFT_Begin || em >= EFT_End) return ToString("Unknown", w, c, left);
  return ToString(em, w, c, left);
}

// -----------------------------------------------------------------------------
template <>
inline bool FromString(const char *str, ExternalForceTerm &eft)
{
  EnergyMeasure em = EM_Unknown;
  if (FromString(str, em) && EFT_Begin < em && em < EFT_End) {
    eft = static_cast<ExternalForceTerm>(em);
    return true;
  } else {
    eft = EFT_Unknown;
    return false;
  }
}


} // namespace mirtk

#endif // MIRTK_ExternalForceTerm_H
