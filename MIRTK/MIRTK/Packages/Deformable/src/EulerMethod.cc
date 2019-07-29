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

#include "mirtk/EulerMethod.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/DeformableSurfaceModel.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/ObjectFactory.h"

#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkUnsignedCharArray.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterOptimizerMacro(EulerMethod);


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace EulerMethodUtils {


// -----------------------------------------------------------------------------
/// Compute node displacements given negated force and time step
class ComputeDisplacements
{
  const double *_Gradient;
  vtkDataArray *_Displacement;
  double        _StepLength;

public:

  ComputeDisplacements(vtkDataArray *dx, const double *gradient, double dt, double norm)
  :
    _Gradient(gradient), _Displacement(dx), _StepLength(-dt / norm)
  {}

  void operator ()(const blocked_range<int> &ptIds) const
  {
    double d[3];
    const double *g = _Gradient + 3 * ptIds.begin();
    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId, g += 3) {
      d[0] = _StepLength * g[0];
      d[1] = _StepLength * g[1];
      d[2] = _StepLength * g[2];
      _Displacement->SetTuple(ptId, d);
    }
  }
};

// -----------------------------------------------------------------------------
/// Clamp magnitudes of node displacements to [0, max]
class ClampDisplacements
{
  vtkDataArray *_Displacement;
  double        _Maximum;

public:

  ClampDisplacements(vtkDataArray *dx, double max_dx)
  :
    _Displacement(dx), _Maximum(max_dx * max_dx)
  {}

  void operator ()(const blocked_range<int> &ptIds) const
  {
    double norm, d[3];
    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      _Displacement->GetTuple(ptId, d);
      norm = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
      if (norm > _Maximum) {
        norm = sqrt(_Maximum / norm);
        d[0] *= norm, d[1] *= norm, d[2] *= norm;
        _Displacement->SetTuple(ptId, d);
      }
    }
  }
};


} // namespace EulerMethodUtils
using namespace EulerMethodUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
EulerMethod::EulerMethod(ObjectiveFunction *f)
:
  LocalOptimizer(f),
  _StepLength(1.0),
  _NormalizeStepLength(true),
  _MaximumDisplacement(.0),
  _Gradient(nullptr),
  _NumberOfDOFs(0)
{
  _Epsilon = 1e-9;
}

// -----------------------------------------------------------------------------
void EulerMethod::CopyAttributes(const EulerMethod &other)
{
  _StepLength          = other._StepLength;
  _NormalizeStepLength = other._NormalizeStepLength;
  _MaximumDisplacement = other._MaximumDisplacement;
  _Displacement        = other._Displacement;
  _NormalDisplacement  = other._NormalDisplacement;
  _NumberOfDOFs        = other._NumberOfDOFs;

  if (_NumberOfDOFs != other._NumberOfDOFs || other._NumberOfDOFs == 0) {
    Deallocate(_Gradient);
    if (_NumberOfDOFs > 0) {
      _Gradient = Allocate<double>(_NumberOfDOFs);
    }
  }
  if (_NumberOfDOFs > 0) {
    memcpy(_Gradient, other._Gradient, _NumberOfDOFs * sizeof(double));
  }
}

// -----------------------------------------------------------------------------
EulerMethod::EulerMethod(const EulerMethod &other)
:
  LocalOptimizer(other),
  _Gradient(nullptr)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
EulerMethod &EulerMethod::operator =(const EulerMethod &other)
{
  if (this != &other) {
    LocalOptimizer::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
EulerMethod::~EulerMethod()
{
  Deallocate(_Gradient);
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool EulerMethod::Set(const char *name, const char *value)
{
  if (strcmp(name, "Deformable surface step length") == 0 ||
      strcmp(name, "Length of steps")                == 0 ||
      strcmp(name, "Maximum length of steps")        == 0) {
    return FromString(value, _StepLength);
  }
  if (strcmp(name, "Normalize length of steps") == 0 ||
      strcmp(name, "Normalise length of steps") == 0 ||
      strcmp(name, "Normalize maximum length of steps") == 0 ||
      strcmp(name, "Normalise maximum length of steps") == 0 ||
      strcmp(name, "Normalize deformable surface step length") == 0 ||
      strcmp(name, "Normalise deformable surface step length") == 0) {
    return FromString(value, _NormalizeStepLength);
  }
  if (strcmp(name, "Maximum deformable surface displacement") == 0 ||
      strcmp(name, "Maximum node displacement")               == 0) {
    return FromString(value, _MaximumDisplacement);
  }
  return LocalOptimizer::Set(name, value);
}

// -----------------------------------------------------------------------------
ParameterList EulerMethod::Parameter() const
{
  ParameterList params = LocalOptimizer::Parameter();
  Insert(params, "Length of steps",           _StepLength);
  Insert(params, "Normalize length of steps", _NormalizeStepLength);
  Insert(params, "Maximum node displacement", _MaximumDisplacement);
  return params;
}

// =============================================================================
// Optimization
// =============================================================================

// -----------------------------------------------------------------------------
void EulerMethod::Initialize()
{
  // Initialize base class
  LocalOptimizer::Initialize();

  // Cast objective function to deformable surface model
  _Model = dynamic_cast<DeformableSurfaceModel *>(_Function);
  if (_Model == nullptr) {
    cerr << "EulerMethod::Initialize: Objective function must be a deformable surface model" << endl;
    exit(1);
  }
  if (_Model->Transformation()) {
    cerr << "EulerMethod::Initialize: Optimizer can only be used for non-parametric deformable surface models" << endl;
    exit(1);
  }

  // Allocate memory for node vectors if not done before
  if (_Model->NumberOfDOFs() > _NumberOfDOFs) {
    Deallocate(_Gradient);
    _NumberOfDOFs = _Model->NumberOfDOFs();
    Allocate(_Gradient, _NumberOfDOFs);
  }

  // Get model point data
  vtkPointData * const modelPD = _Model->Output()->GetPointData();

  // Add point data array used to keep track of active/passive nodes
  //
  // Note that the status of a node may be set to passive (i.e., status = .0)
  // during the evaluation of a stopping criterion such as in particular the
  // MinActiveStoppingCriterion. It will then be excluded from further
  // changes (and force evaluation) of the deformable surface model.
  //
  // A "Status" array can also be provided as input mask to restrict the
  // movement of model nodes to a subset only (e.g., cortical vertices).
  vtkSmartPointer<vtkDataArray> status;
  status = modelPD->GetArray("Status");
  if (!status) {
    status = vtkSmartPointer<vtkUnsignedCharArray>::New();
    status->SetName("Status");
    status->SetNumberOfComponents(1);
    status->SetNumberOfTuples(_Model->NumberOfPoints());
    status->FillComponent(0, 1.0);
    modelPD->AddArray(status);
  }

  // Add point data array for node update such that remesher can modify it
  // when necessary before it is passed on to the convergence test functions
  //
  // Moreover, this displacement array is needed by the EulerMethodWithMomentum.
  // An initial node "Displacement" can also be provided as input for this method
  // (e.g., from a previous Euler integration with different parameters).
  _Displacement = modelPD->GetArray("Displacement");
  if (!_Displacement || _Displacement->GetNumberOfComponents() != 3) {
    _Displacement = vtkSmartPointer<vtkDoubleArray>::New();
    _Displacement->SetName("Displacement");
    _Displacement->SetNumberOfComponents(3);
    _Displacement->SetNumberOfTuples(_Model->NumberOfPoints());
    _Displacement->FillComponent(0, 0.);
    _Displacement->FillComponent(1, 0.);
    _Displacement->FillComponent(2, 0.);
    modelPD->RemoveArray(_Displacement->GetName());
    modelPD->AddArray(_Displacement);
  } else if (_Displacement->GetDataType() != VTK_DOUBLE) {
    vtkSmartPointer<vtkDataArray> displacement = _Displacement;
    _Displacement = vtkSmartPointer<vtkDoubleArray>::New();
    _Displacement->SetName("Displacement");
    _Displacement->SetNumberOfComponents(3);
    _Displacement->SetNumberOfTuples(_Model->NumberOfPoints());
    _Displacement->CopyComponent(0, displacement, 0);
    _Displacement->CopyComponent(1, displacement, 1);
    _Displacement->CopyComponent(2, displacement, 2);
    modelPD->RemoveArray(displacement->GetName());
    modelPD->AddArray(_Displacement);
  }

  // Add point data array used to keep track of node displacement in
  // normal direction such that the remesher can modify it when necessary
  if (_NormalDisplacement) {
    if (_NormalDisplacement->GetName() == nullptr) {
      _NormalDisplacement->SetName("NormalDisplacement");
    }
    int i;
    for (i = 0; i < modelPD->GetNumberOfArrays(); ++i) {
      if (modelPD->GetArray(i) == _NormalDisplacement) break;
    }
    if (i == modelPD->GetNumberOfArrays()) {
      modelPD->RemoveArray(_NormalDisplacement->GetName());
      modelPD->AddArray(_NormalDisplacement);
    }
  }
}

// -----------------------------------------------------------------------------
double EulerMethod::Run()
{
  double *dx;

  // Initialize
  this->Initialize();

  // Initial update of deformable surface model before start event because
  // the update can trigger some lazy initialization which in turn may
  // broadcast some log events for verbose command output
  _Model->Update(true);

  // Notify observers about start of optimization
  Broadcast(StartEvent);

  // Get initial energy value
  double value = _Model->Value();
  _LastValues.clear();
  _LastValues.push_back(value);

  // Perform explicit integration steps
  _Converged = false;
  Iteration step(0, _NumberOfSteps);
  while (!_Converged && step.Next()) {

    // Notify observers about start of iteration
    Broadcast(IterationStartEvent, &step);

    // Compute (negative) node forces
    _Model->Gradient(_Gradient);

    // Update current node displacements
    this->UpdateDisplacement();
    dx = static_cast<double *>(_Displacement->GetVoidPointer(0));

    // Perform time step
    _LastDelta = _Model->Step(dx);
    if (_LastDelta <= _Delta) break;

    // Track node displacement in normal direction
    // (e.g., sulcal depth measure during cortical surface inflation)
    this->UpdateNormalDisplacement();

    // Perform local adaptive remeshing
    if (this->RemeshModel()) {
      dx = static_cast<double *>(_Displacement->GetVoidPointer(0));
    }

    // Update model terms
    _Model->Update(true);

    // Test stopping criteria
    //
    // For deformable surface models which are not defined as a local minimum
    // of a well-defined energy function, but only via the equilibrium of
    // internal and external forces, the energy values corresponding to the
    // external forces are infinite and hence the total energy value.
    if (!IsInf(value)) value = _Model->Value();
    _Converged = Converged(step.Iter(), value, dx);

    // Notify observers about end of iteration
    Broadcast(IterationEndEvent, &step);
  }

  // Notify observers about end of optimization
  Broadcast(EndEvent, &value);

  // Finalize
  this->Finalize();

  return value;
}

// -----------------------------------------------------------------------------
bool EulerMethod::RemeshModel()
{
  if (_Model->Remesh()) {
    if (_Model->NumberOfDOFs() > _NumberOfDOFs) {
      Deallocate(_Gradient);
      _NumberOfDOFs = _Model->NumberOfDOFs();
      Allocate(_Gradient, _NumberOfDOFs);
    }
    vtkPointData * const modelPD = _Model->Output()->GetPointData();
    _Displacement = modelPD->GetArray(_Displacement->GetName());
    if (_NormalDisplacement) {
      _NormalDisplacement = modelPD->GetArray(_NormalDisplacement->GetName());
    }
    return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
double EulerMethod::GradientNorm() const
{
  if (_NormalizeStepLength) {
    const double norm = _Model->GradientNorm(_Gradient);
    return (norm > .0 ? norm : 1.0);
  } else {
    // Maximum displacement limited per node using _MaximumDisplacement instead.
    // Node forces are normalized by the number of mesh nodes such that when used
    // within an image + surface registration, their magnitude is compatible with
    // the image similarity gradient. In case of a deformable surface model,
    // however, we want to revert this normalization of the force magnitude here.
    return 1.0 / _Model->NumberOfPoints();
  }
}

// -----------------------------------------------------------------------------
void EulerMethod::UpdateDisplacement()
{
  ComputeDisplacements eval(_Displacement, _Gradient, _StepLength, this->GradientNorm());
  parallel_for(blocked_range<int>(0, _Model->NumberOfPoints()), eval);
  this->TruncateDisplacement();
}

// -----------------------------------------------------------------------------
void EulerMethod::TruncateDisplacement(bool force)
{
  double max_dx = _MaximumDisplacement;
  if (max_dx <= .0) max_dx = _NormalizeStepLength ? _StepLength : 1.0;
  if (force || !_NormalizeStepLength || max_dx < _StepLength) {
    ClampDisplacements clamp(_Displacement, max_dx);
    parallel_for(blocked_range<int>(0, _Model->NumberOfPoints()), clamp);
  }
}

// -----------------------------------------------------------------------------
void EulerMethod::UpdateNormalDisplacement()
{
  if (_NormalDisplacement && IsSurfaceMesh(_Model->Output())) {
    double m, n[3], dx[3];
    vtkDataArray * const normals = _Model->PointSet().SurfaceNormals();
    for (int ptId = 0; ptId < _Model->NumberOfPoints(); ++ptId) {
      normals->GetTuple(ptId, n);
      _Displacement->GetTuple(ptId, dx);
      m  = _NormalDisplacement->GetComponent(ptId, 0);
      m += dx[0]*n[0] + dx[1]*n[1] + dx[2]*n[2];
      _NormalDisplacement->SetComponent(ptId, 0, m);
    }
  }
}

// -----------------------------------------------------------------------------
void EulerMethod::Finalize()
{
}


} // namespace mirtk
