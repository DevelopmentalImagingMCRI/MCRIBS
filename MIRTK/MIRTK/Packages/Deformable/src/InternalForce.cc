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

#include "mirtk/InternalForce.h"


namespace mirtk {


// =============================================================================
// Factory method
// =============================================================================

// -----------------------------------------------------------------------------
InternalForce *InternalForce::New(InternalForceTerm ift, const char *name, double w)
{
  enum EnergyMeasure em = static_cast<enum EnergyMeasure>(ift);
  if (IFT_Begin < em && em < IFT_End) {
    EnergyTerm *term = EnergyTerm::TryNew(em, name, w);
    if (term) return dynamic_cast<InternalForce *>(term);
    ThrowStatic(ERR_RuntimeError, NameOfType(), __FUNCTION__,
                "Internal point set force not available: ", em, " (", int(em), ")");
  } else {
    ThrowStatic(ERR_RuntimeError, NameOfType(), __FUNCTION__,
                "Energy term is not an internal point set force: ", em, " (", int(em), ")");
  }
  return nullptr;
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void InternalForce::CopyAttributes(const InternalForce &other)
{
  _ExternalMagnitudeArrayName = other._ExternalMagnitudeArrayName;
  _WeightInside               = other._WeightInside;
  _WeightOutside              = other._WeightOutside;
  _WeightMinimum              = other._WeightMinimum;
}

// -----------------------------------------------------------------------------
InternalForce::InternalForce(const char *name, double weight)
:
  PointSetForce(name, weight),
  _WeightInside (1.),
  _WeightOutside(1.),
  _WeightMinimum(0.)
{
}

// -----------------------------------------------------------------------------
InternalForce::InternalForce(const InternalForce &other)
:
  PointSetForce(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
InternalForce &InternalForce::operator =(const InternalForce &other)
{
  if (this != &other) {
    PointSetForce::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
InternalForce::~InternalForce()
{
}

// =============================================================================
// Configuration
// =============================================================================

// -----------------------------------------------------------------------------
bool InternalForce::SetWithoutPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Weight inside") == 0 ||
      strcmp(param, "Weight factor inside") == 0) {
    double weight;
    if (!FromString(value, weight) || weight < 0.) return false;
    _WeightInside = weight;
    return true;
  }
  if (strcmp(param, "Weight outside") == 0 ||
      strcmp(param, "Weight factor outside") == 0) {
    double weight;
    if (!FromString(value, weight) || weight < 0.) return false;
    _WeightOutside = weight;
    return true;
  }
  if (strcmp(param, "Minimum weight") == 0 ||
      strcmp(param, "Weight minimum") == 0 ||
      strcmp(param, "Minimum weight factor") == 0 ||
      strcmp(param, "Weight factor minimum") == 0) {
    double weight;
    if (!FromString(value, weight) || weight < 0.) return false;
    _WeightMinimum = weight;
    return true;
  }
  return PointSetForce::SetWithoutPrefix(param, value);
}

// -----------------------------------------------------------------------------
ParameterList InternalForce::Parameter() const
{
  ParameterList params = PointSetForce::Parameter();
  InsertWithPrefix(params, "Weight factor inside",  _WeightInside);
  InsertWithPrefix(params, "Weight factor outside", _WeightOutside);
  InsertWithPrefix(params, "Minimum weight factor", _WeightMinimum);
  return params;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
vtkDataArray *InternalForce::ExternalMagnitude() const
{
  vtkDataArray *mag = nullptr;
  if (!_ExternalMagnitudeArrayName.empty()) {
    const bool optional = false;
    mag = PointData(_ExternalMagnitudeArrayName.c_str(), optional);
    if (mag->GetNumberOfComponents() != 1) {
      Throw(ERR_LogicError, __FUNCTION__, "External force magnitude array must"
                                          " have only one scalar component!");
    }
  }
  return mag;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void InternalForce::EvaluateGradient(double *gradient, double step, double weight)
{
  // Scale by external force magnitude
  vtkDataArray * const scale = ExternalMagnitude();

  if (scale != nullptr) {
    double s;
    GradientType *grad = _Gradient;
    for (int i = 0; i < _NumberOfPoints; ++i, ++grad) {
      s  = scale->GetComponent(i, 0);
      s *= (s < 0. ? _WeightOutside : _WeightInside);
      (*grad) *= _WeightMinimum + abs(s);
    }
  }

  // Compute parametric gradient
  PointSetForce::EvaluateGradient(gradient, step, weight);
}


} // namespace mirtk
