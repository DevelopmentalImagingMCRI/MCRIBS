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

#include "mirtk/EulerMethodWithDamping.h"

#include "mirtk/Parallel.h"
#include "mirtk/DeformableSurfaceModel.h"
#include "mirtk/ObjectFactory.h"

#include "vtkPointData.h"
#include "vtkFloatArray.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterOptimizerMacro(EulerMethodWithDamping);


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace EulerMethodWithDampingUtils {


// -----------------------------------------------------------------------------
/// Compute node displacements given previous velocity, negated force and time step
class ComputeDisplacements
{
  const double *_Gradient;
  vtkDataArray *_Velocity;
  vtkDataArray *_Displacement;
  double        _DampingFactor;
  double        _BodyMass;
  double        _StepLength;
  double        _GradientNorm;

public:

  ComputeDisplacements(vtkDataArray *dx, vtkDataArray *v, const double *gradient,
                       double damping, double mass, double dt, double norm)
  :
    _Gradient(gradient),
    _Velocity(v),
    _Displacement(dx),
    _DampingFactor(damping),
    _BodyMass(mass),
    _StepLength(dt),
    _GradientNorm(norm)
  {}

  void operator ()(const blocked_range<int> &ptIds) const
  {
    double v[3], d[3];
    const double s1 = -_StepLength / (_BodyMass * _GradientNorm);
    const double s2 =  _StepLength * _DampingFactor / _BodyMass;
    const double *g = _Gradient     + 3 * ptIds.begin();
    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId, g += 3) {
      _Velocity->GetTuple(ptId, v);
      v[0] += s1 * g[0] - s2 * v[0];
      v[1] += s1 * g[1] - s2 * v[1];
      v[2] += s1 * g[2] - s2 * v[2];
      _Velocity->SetTuple(ptId, v);
      d[0] = v[0] * _StepLength;
      d[1] = v[1] * _StepLength;
      d[2] = v[2] * _StepLength;
      _Displacement->SetTuple(ptId, d);
    }
  }
};


} // namespace EulerMethodWithDampingUtils
using namespace EulerMethodWithDampingUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
EulerMethodWithDamping::EulerMethodWithDamping(ObjectiveFunction *f)
:
  EulerMethod(f),
  _BodyMass(1.0),
  _DampingFactor(.1)
{
}

// -----------------------------------------------------------------------------
void EulerMethodWithDamping::CopyAttributes(const EulerMethodWithDamping &other)
{
  _DampingFactor = other._DampingFactor;
  _BodyMass      = other._BodyMass;
}

// -----------------------------------------------------------------------------
EulerMethodWithDamping::EulerMethodWithDamping(const EulerMethodWithDamping &other)
:
  EulerMethod(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
EulerMethodWithDamping &EulerMethodWithDamping::operator =(const EulerMethodWithDamping &other)
{
  if (this != &other) {
    EulerMethod::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
EulerMethodWithDamping::~EulerMethodWithDamping()
{
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool EulerMethodWithDamping::Set(const char *name, const char *value)
{
  if (strcmp(name, "Deformable surface mass") == 0) {
    return FromString(value, _BodyMass);
  }
  if (strcmp(name, "Deformable surface damping") == 0) {
    return FromString(value, _DampingFactor);
  }
  if (strcmp(name, "Deformable surface momentum") == 0) {
    double momentum;
    if (!FromString(value, momentum) || momentum < .0 || momentum > 1.0) return false;
    _DampingFactor = 1.0 - momentum;
    return true;
  }
  return EulerMethod::Set(name, value);
}

// -----------------------------------------------------------------------------
ParameterList EulerMethodWithDamping::Parameter() const
{
  ParameterList params = EulerMethod::Parameter();
  Insert(params, "Deformable surface mass",    _BodyMass);
  Insert(params, "Deformable surface damping", _DampingFactor);
  return params;
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void EulerMethodWithDamping::Initialize()
{
  // Initialize base class
  EulerMethod::Initialize();

  // Get model point data
  vtkPointData *modelPD = _Model->Output()->GetPointData();

  // Add point data array with initial node velocities such that these
  // are interpolated at new node positions during the remeshing
  //
  // An initial node "Velocity" can also be provided as input
  // (e.g., from a previous Euler integration with different parameters).
  vtkSmartPointer<vtkDataArray> velocity;
  velocity = modelPD->GetArray("Velocity");
  if (!velocity) {
    velocity = vtkSmartPointer<vtkFloatArray>::New();
    velocity->SetName("Velocity");
    velocity->SetNumberOfComponents(3);
    velocity->SetNumberOfTuples(_Model->NumberOfPoints());
    velocity->FillComponent(0, .0);
    velocity->FillComponent(1, .0);
    velocity->FillComponent(2, .0);
    modelPD->AddArray(velocity);
  }

  // Limit damping factor to the interval [0, 1]
  _DampingFactor = max(.0, min(_DampingFactor, 1.0));
}

// -----------------------------------------------------------------------------
void EulerMethodWithDamping::UpdateDisplacement()
{
  vtkPointData *modelPD  = _Model->Output()->GetPointData();
  vtkDataArray *velocity = modelPD->GetArray("Velocity");
  ComputeDisplacements eval(_Displacement, velocity, _Gradient, _DampingFactor,
                            _BodyMass, _StepLength, this->GradientNorm());
  parallel_for(blocked_range<int>(0, _Model->NumberOfPoints()), eval);
  this->TruncateDisplacement();
}


} // namespace mirtk
