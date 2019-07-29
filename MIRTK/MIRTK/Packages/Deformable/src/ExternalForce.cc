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

#include "mirtk/ExternalForce.h"

#include "mirtk/Math.h"


namespace mirtk {


// =============================================================================
// Factory
// =============================================================================

// -----------------------------------------------------------------------------
ExternalForce *ExternalForce::New(ExternalForceTerm eft, const char *name, double w)
{
  enum EnergyMeasure em = static_cast<enum EnergyMeasure>(eft);
  if (EFT_Begin < em && em < EFT_End) {
    EnergyTerm *term = EnergyTerm::TryNew(em, name, w);
    if (term) return dynamic_cast<ExternalForce *>(term);
    cerr << NameOfType() << "::New: External point set force not available: ";
  } else {
    cerr << NameOfType() << "::New: Energy term is not an external point set force: ";
  }
  cerr << ToString(em) << " (" << em << ")" << endl;
  exit(1);
  return NULL;
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
ExternalForce::ExternalForce(const char *name, double weight)
:
  PointSetForce(name, weight),
  _Image(NULL)
{
}

// -----------------------------------------------------------------------------
void ExternalForce::CopyAttributes(const ExternalForce &other)
{
  _Image = other._Image;
}

// -----------------------------------------------------------------------------
ExternalForce::ExternalForce(const ExternalForce &other)
:
  PointSetForce(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
ExternalForce &ExternalForce::operator =(const ExternalForce &other)
{
  if (this != &other) {
    PointSetForce::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
ExternalForce::~ExternalForce()
{
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
void ExternalForce::Initialize()
{
  // Initialize base class
  PointSetForce::Initialize();
  if (_NumberOfPoints == 0) return;

  // Check input
  if (_Image == NULL) {
    cerr << "ExternalForce::Initialize: Reference image is NULL" << endl;
    exit(1);
  }
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void ExternalForce::Update(bool gradient)
{
  // Before PointSetForce::Update which sets _InitialUpdate = false
  if (_InitialUpdate || _Image->Transformation()) {
    _Image->Update(_InitialUpdate && _Image->SelfUpdate());
  }
  PointSetForce::Update(gradient);
}

// -----------------------------------------------------------------------------
double ExternalForce::Evaluate()
{
  // Not all external forces may have a corresponding energy term expression
  // if the deformable surface model is rather based on reaching an equilibrium
  // of internal and external forces. In this case, return infinity such that
  // the optimizer only checks the equilibrium condition
  return numeric_limits<double>::infinity();
}


} // namespace mirtk
