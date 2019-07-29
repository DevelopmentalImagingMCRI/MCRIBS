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

#ifndef MIRTK_SpringForce_H
#define MIRTK_SpringForce_H

#include "mirtk/SurfaceConstraint.h"
#include "mirtk/EdgeTable.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkDataArray.h"


namespace mirtk {


/**
 * Surface smoothness force based on discrete Laplace operator
 *
 * This internal deformable surface force is based on the umbrella
 * operator which approximates the discrete Laplace operator.
 *
 * It corresponds to the "spring forces" employed by FreeSurfer's surface tools:
 *
 *   Dale et al., Cortical Surface-Based Analysis I:
 *   Segmentation and Surface Reconstruction. NeuroImage, 9(2), 179–194 (1999).
 *
 *   Fischl et al., Cortical Surface-Based Analysis II: Inflation, Flattening,
 *   and a Surface-Based Coordinate System. NeuroImage, 9(2), 195–207 (1999).
 *
 * And the "tensile force" used by McInerney and Terzopoulos:
 *
 *   McInerney and Terzopoulos, Topology adaptive deformable surfaces for medical
 *   image volume segmentation. IEEE Transactions on Medical Imaging, 18(10),
 *   840–850, doi:10.1109/42.811261 (1999)
 *
 * \sa Other smoothness terms: CurvatureConstraint, QuadraticCurvatureConstraint
 */
class SpringForce : public SurfaceConstraint
{
  mirtkEnergyTermMacro(SpringForce, EM_SpringForce);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Weight of normal component
  mirtkPublicAttributeMacro(double, InwardNormalWeight);

  /// Weight of normal component
  mirtkPublicAttributeMacro(double, OutwardNormalWeight);

  /// Weight of tangential component
  mirtkPublicAttributeMacro(double, TangentialWeight);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const SpringForce &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  SpringForce(const char * = "", double = 1.0);

  /// Copy constructor
  SpringForce(const SpringForce &);

  /// Assignment operator
  SpringForce &operator =(const SpringForce &);

  /// Destructor
  virtual ~SpringForce();

  // ---------------------------------------------------------------------------
  // Evaluation

protected:

  /// Evaluate energy of internal force term
  virtual double Evaluate();

  /// Evaluate internal force w.r.t. transformation parameters or surface nodes
  virtual void EvaluateGradient(double *, double, double);

};


} // namespace mirtk

#endif // MIRTK_SpringForce_H
