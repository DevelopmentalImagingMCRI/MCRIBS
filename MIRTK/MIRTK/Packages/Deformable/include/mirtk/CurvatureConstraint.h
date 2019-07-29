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

#ifndef MIRTK_CurvatureConstraint_H
#define MIRTK_CurvatureConstraint_H

#include "mirtk/SurfaceConstraint.h"

#include "mirtk/List.h"
#include "mirtk/Array.h"

#include "vtkSmartPointer.h"
#include "vtkPoints.h"


namespace mirtk {


/**
 * Surface bending / smoothness constraint
 *
 * This internal deformable surface force term is similar to the discrete
 * Laplacian operator also used for iteratively smoothing surfaces (cf. the
 * "tensile force" used by McInerney & Terzopoulos), but not identical.
 *
 * McInerney & Terzopoulos, Topology adaptive deformable surfaces for medical
 * image volume segmentation. IEEE Transactions on Medical Imaging, 18(10),
 * 840–850, doi:10.1109/42.811261 (1999)
 *
 * Lachaud and Montanvert, Deformable meshes with automated topology changes
 * for coarse-to-fine three-dimensional surface extraction. Medical Image Analysis,
 * 3(2), 187–207, doi:10.1016/S1361-8415(99)80012-7 (1999)
 *
 * Park et al., A non-self-intersecting adaptive deformable surface for
 * complex boundary extraction from volumetric images, 25, 421–440 (2001).
 */
class CurvatureConstraint : public SurfaceConstraint
{
  mirtkEnergyTermMacro(CurvatureConstraint, EM_Curvature);

  // ---------------------------------------------------------------------------
  // Attributes

protected:

  /// Centroids of adjacent nodes
  mirtkAttributeMacro(vtkSmartPointer<vtkPoints>, Centroids);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const CurvatureConstraint &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  CurvatureConstraint(const char * = "", double = 1.0);

  /// Copy constructor
  CurvatureConstraint(const CurvatureConstraint &);

  /// Assignment operator
  CurvatureConstraint &operator =(const CurvatureConstraint &);

  /// Destructor
  virtual ~CurvatureConstraint();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Initialize internal force term
  virtual void Initialize();

  /// Reinitialize internal force term after change of input topology
  virtual void Reinitialize();

  /// Update internal force data structures
  virtual void Update(bool);

protected:

  /// Common (re-)initialization code of this class only (non-virtual function!)
  void Init();

  /// Compute penalty for current transformation estimate
  virtual double Evaluate();

  /// Compute internal force w.r.t. transformation parameters
  virtual void EvaluateGradient(double *, double, double);

};


} // namespace mirtk

#endif // MIRTK_CurvatureConstraint_H
