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

#include "mirtk/EulerMethodWithMomentum.h"

#include "mirtk/Math.h"
#include "mirtk/Parallel.h"
#include "mirtk/DeformableSurfaceModel.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/ObjectFactory.h"

#include "vtkPointData.h"
#include "vtkFloatArray.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterOptimizerMacro(EulerMethodWithMomentum);


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace EulerMethodWithMomentumUtils {


// -----------------------------------------------------------------------------
/// Compute node displacements given previous velocity, negated force and time step
class ComputeDisplacements
{
  const double *_Gradient;
  vtkDataArray *_Displacement;
  double        _Momentum;
  double        _Maximum;
  double        _StepLength;

public:

  ComputeDisplacements(vtkDataArray *dx, const double *gradient,
                       double momentum, double max_dx, double dt, double norm)
  :
    _Gradient(gradient),
    _Displacement(dx),
    _Momentum(momentum),
    _Maximum(max_dx * max_dx),
    _StepLength(-dt / norm)
  {}

  void operator ()(const blocked_range<int> &ptIds) const
  {
    double norm, d[3];
    const double *g = _Gradient + 3 * ptIds.begin();
    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId, g += 3) {
      _Displacement->GetTuple(ptId, d);
      d[0] = _StepLength * g[0] + _Momentum * d[0];
      d[1] = _StepLength * g[1] + _Momentum * d[1];
      d[2] = _StepLength * g[2] + _Momentum * d[2];
      norm = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
      if (norm > _Maximum) {
        norm = sqrt(_Maximum / norm);
        d[0] *= norm, d[1] *= norm, d[2] *= norm;
      }
      _Displacement->SetTuple(ptId, d);
    }
  }
};


} // namespace EulerMethodWithMomentumUtils
using namespace EulerMethodWithMomentumUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
EulerMethodWithMomentum::EulerMethodWithMomentum(ObjectiveFunction *f)
:
  EulerMethod(f),
  _Momentum(.9),
  _ExcludeMomentumFromNormalDisplacement(false)
{
}

// -----------------------------------------------------------------------------
void EulerMethodWithMomentum::CopyAttributes(const EulerMethodWithMomentum &other)
{
  _Momentum                              = other._Momentum;
  _ExcludeMomentumFromNormalDisplacement = other._ExcludeMomentumFromNormalDisplacement;
}

// -----------------------------------------------------------------------------
EulerMethodWithMomentum::EulerMethodWithMomentum(const EulerMethodWithMomentum &other)
:
  EulerMethod(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
EulerMethodWithMomentum &EulerMethodWithMomentum::operator =(const EulerMethodWithMomentum &other)
{
  if (this != &other) {
    EulerMethod::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
EulerMethodWithMomentum::~EulerMethodWithMomentum()
{
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool EulerMethodWithMomentum::Set(const char *name, const char *value)
{
  if (strcmp(name, "Deformable surface momentum") == 0) {
    return FromString(value, _Momentum) && .0 <= _Momentum && _Momentum <= 1.0;
  }
  if (strcmp(name, "Deformable surface damping") == 0) {
    double damping;
    if (!FromString(value, damping) || damping < .0 || damping > 1.0) return false;
    _Momentum = 1.0 - damping;
    return true;
  }
  if (strcmp(name, "Exclude momentum from tracked normal displacement") == 0) {
    return FromString(value, _ExcludeMomentumFromNormalDisplacement);
  }
  return EulerMethod::Set(name, value);
}

// -----------------------------------------------------------------------------
ParameterList EulerMethodWithMomentum::Parameter() const
{
  ParameterList params = EulerMethod::Parameter();
  Insert(params, "Deformable surface momentum", _Momentum);
  Insert(params, "Exclude momentum from tracked normal displacement", _ExcludeMomentumFromNormalDisplacement);
  return params;
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void EulerMethodWithMomentum::Initialize()
{
  // Initialize base class
  EulerMethod::Initialize();

  // Limit momentum factor to the interval [0, 1]
  _Momentum = max(0., min(_Momentum, 1.));
}

// -----------------------------------------------------------------------------
void EulerMethodWithMomentum::UpdateDisplacement()
{
  double norm   = this->GradientNorm();
  double max_dx = _MaximumDisplacement;
  if (max_dx <= .0) max_dx = _NormalizeStepLength ? _StepLength : 1.0;
  ComputeDisplacements eval(_Displacement, _Gradient, _Momentum, max_dx, _StepLength, norm);
  parallel_for(blocked_range<int>(0, _Model->NumberOfPoints()), eval);
}

// -----------------------------------------------------------------------------
void EulerMethodWithMomentum::UpdateNormalDisplacement()
{
  if (_ExcludeMomentumFromNormalDisplacement) {
    // Note: The actual displacement of a node is given by _Displacement.
    //       However, the tracking of the normal displacements is in particular
    //       used to measure sulcal depth during cortical surface inflation.
    //       The FreeSurfer function MRISinflateBrain only sums up the current
    //       displacements (scaled forces) excluding the momentum.
    if (_NormalDisplacement && IsSurfaceMesh(_Model->Output())) {
      double d, n[3];
      const double dt = - _StepLength / this->GradientNorm();
      const double *g = _Gradient;
      vtkDataArray *normals = _Model->PointSet().SurfaceNormals();
      for (int ptId = 0; ptId < _Model->NumberOfPoints(); ++ptId, g += 3) {
        normals->GetTuple(ptId, n);
        d  = _NormalDisplacement->GetComponent(ptId, 0);
        d += dt*g[0]*n[0] + dt*g[1]*n[1] + dt*g[2]*n[2];
        _NormalDisplacement->SetComponent(ptId, 0, d);
      }
    }
  } else {
    EulerMethod::UpdateNormalDisplacement();
  }
}


} // namespace mirtk
