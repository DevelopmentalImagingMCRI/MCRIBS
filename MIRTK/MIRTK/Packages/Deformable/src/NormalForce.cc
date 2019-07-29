/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
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

#include "mirtk/NormalForce.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(NormalForce);


// =============================================================================
// Auxiliaries
// =============================================================================

namespace NormalForceUtils {


// -----------------------------------------------------------------------------
/// Evaluate gradient of bending penalty (i.e., negative internal bending force)
struct EvaluateGradient
{
  typedef NormalForce::GradientType Force;

  vtkDataArray *_Status;
  vtkDataArray *_Normals;
  Force        *_Gradient;

  void operator ()(const blocked_range<int> &re) const
  {
    double n[3];
    for (int ptId = re.begin(); ptId != re.end(); ++ptId) {
      if (_Status && _Status->GetComponent(ptId, 0) == .0) {
        n[0] = n[1] = n[2] = 0.;
      } else {
        _Normals->GetTuple(ptId, n);
      }
      _Gradient[ptId] = -Force(n);
    }
  }
};


} // namespace NormalForceUtils

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
NormalForce::NormalForce(const char *name, double weight)
:
  InternalForce(name, weight)
{
  _SurfaceForce = true;
}

// -----------------------------------------------------------------------------
NormalForce::NormalForce(const NormalForce &other)
:
  InternalForce(other)
{
}

// -----------------------------------------------------------------------------
NormalForce &NormalForce::operator =(const NormalForce &other)
{
  InternalForce::operator =(other);
  return *this;
}

// -----------------------------------------------------------------------------
NormalForce::~NormalForce()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double NormalForce::Evaluate()
{
  return numeric_limits<double>::infinity();
}

// -----------------------------------------------------------------------------
void NormalForce::EvaluateGradient(double *gradient, double step, double weight)
{
  if (_NumberOfPoints == 0) return;

  MIRTK_START_TIMING();
  NormalForceUtils::EvaluateGradient eval;
  eval._Status   = Status();
  eval._Normals  = Normals();
  eval._Gradient = _Gradient;
  parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);
  InternalForce::EvaluateGradient(gradient, step, weight / _NumberOfPoints);
  MIRTK_DEBUG_TIMING(3, "evaluation of normal force");
}


} // namespace mirtk
