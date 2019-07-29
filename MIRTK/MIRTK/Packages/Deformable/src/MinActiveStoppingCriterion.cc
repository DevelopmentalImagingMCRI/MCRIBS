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

#include "mirtk/MinActiveStoppingCriterion.h"

#include "mirtk/LocalOptimizer.h"
#include "mirtk/DeformableSurfaceModel.h"

#include "vtkPointSet.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
MinActiveStoppingCriterion
::MinActiveStoppingCriterion(const ObjectiveFunction *f)
:
  StoppingCriterion(f),
  _Threshold(.01),
  _StreakOfPassiveIterations(5),
  _Active(1.0)
{
}

// -----------------------------------------------------------------------------
void MinActiveStoppingCriterion
::CopyAttributes(const MinActiveStoppingCriterion &other)
{
  _Threshold                 = other._Threshold;
  _StreakOfPassiveIterations = other._StreakOfPassiveIterations;
  _Active                    = other._Active;
}

// -----------------------------------------------------------------------------
MinActiveStoppingCriterion
::MinActiveStoppingCriterion(const MinActiveStoppingCriterion &other)
:
  StoppingCriterion(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
MinActiveStoppingCriterion &MinActiveStoppingCriterion
::operator =(const MinActiveStoppingCriterion &other)
{
  if (this != &other) {
    StoppingCriterion::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
StoppingCriterion *MinActiveStoppingCriterion::New() const
{
  return new MinActiveStoppingCriterion(*this);
}

// -----------------------------------------------------------------------------
MinActiveStoppingCriterion::~MinActiveStoppingCriterion()
{
}

// -----------------------------------------------------------------------------
void MinActiveStoppingCriterion::Initialize()
{
  // Initialize base class
  StoppingCriterion::Initialize();

  // Add point data array to keep track of when a node was last modified
  const DeformableSurfaceModel *model;
  model = dynamic_cast<const DeformableSurfaceModel *>(_Function);
  if (model && model->NumberOfPoints() > 0 && !model->Transformation()) {
    vtkSmartPointer<vtkDataArray> modified;
    vtkPointData *modelPD = model->Output()->GetPointData();
    modified = modelPD->GetArray("LastModified");

    if (!modified) {
      modified = vtkSmartPointer<vtkIntArray>::New();
      modified->SetName("LastModified");
      modified->SetNumberOfComponents(1);
      modelPD->AddArray(modified);
    }
    modified->SetNumberOfTuples(model->NumberOfPoints());
    modified->FillComponent(0, .0);
    // Preserve node status from a previous optimization
    // (cf. EulerMethod and deformmesh tool)
    vtkDataArray *status = modelPD->GetArray("Status");
    if (status) {
      for (vtkIdType ptId = 0; ptId < model->Output()->GetNumberOfPoints(); ++ptId) {
        if (status->GetComponent(ptId, 0) == .0) {
          modified->SetComponent(ptId, 0, -_StreakOfPassiveIterations);
        }
      }
    }
  }
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
bool MinActiveStoppingCriterion::Fulfilled(int iter, double, const double *dx)
{
  const DeformableSurfaceModel *model;
  model = dynamic_cast<const DeformableSurfaceModel *>(_Function);
  if (!model) return false;

  if (model->NumberOfPoints() == 0) return true;

  vtkPointData *modelPD  = model->Output()->GetPointData();
  vtkDataArray *modified = modelPD->GetArray("LastModified");
  vtkDataArray *status   = modelPD->GetArray("Status");

  if (!modified) return false;

  // Minimum _LastModified iteration for a node to be considered still active
  const int    min_iter   = iter - _StreakOfPassiveIterations;
  const double min_delta2 = _Optimizer->Delta() * _Optimizer->Delta();

  // Record iteration when a node was last modified
  int nactive = model->NumberOfPoints();
  for (int ptId = 0; ptId < model->NumberOfPoints(); ++ptId, dx += 3) {
    if (static_cast<int>(modified->GetComponent(ptId, 0)) < min_iter) {
      --nactive;
      if (status) status->SetComponent(ptId, 0, .0);
    } else if ((dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]) >= min_delta2) {
      modified->SetComponent(ptId, 0, static_cast<double>(iter));
    }
  }

  // Ratio of active nodes
  _Active = double(nactive) / model->NumberOfPoints();

  return _Active <= _Threshold;
}

// =============================================================================
// Logging
// =============================================================================

// -----------------------------------------------------------------------------
void MinActiveStoppingCriterion::Print(ostream &os) const
{
  const ios::fmtflags fmt = os.flags();
  const size_t ndigits = os.precision(1);
  os << "active = " << fixed << setw(5) << (100.0 * _Active) << "%";
  os.precision(ndigits);
  os.flags(fmt);
}


} // namespace mirtk
