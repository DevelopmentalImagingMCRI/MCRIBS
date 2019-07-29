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

#ifndef MIRTK_NearOptimalIntrinsicSurfaceMapper_H
#define MIRTK_NearOptimalIntrinsicSurfaceMapper_H

#include "mirtk/IntrinsicSurfaceMapper.h"

#include "vtkSmartPointer.h"
#include "vtkDataArray.h"


namespace mirtk {


/**
 * Base class of filters that compute a near-optimal intrinsic surface map
 *
 * Surface map filters of this type find the optimal affine combination weights
 * of independently computed conformal and authalic surface maps with given fixed
 * boundary values which minimize a certain non-linear distortion criterion.
 *
 * - Desbrun, Meyer, and Alliez (2002). Intrinsic parameterizations of surface meshes.
 *   Computer Graphics Forum, 21(3), 209–218.
 */
class NearOptimalIntrinsicSurfaceMapper : public IntrinsicSurfaceMapper
{
  mirtkAbstractMacro(NearOptimalIntrinsicSurfaceMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Copy attributes of this class from another instance
  void CopyAttributes(const NearOptimalIntrinsicSurfaceMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  NearOptimalIntrinsicSurfaceMapper();

  /// Copy constructor
  NearOptimalIntrinsicSurfaceMapper(const NearOptimalIntrinsicSurfaceMapper &);

  /// Assignment operator
  NearOptimalIntrinsicSurfaceMapper &operator =(const NearOptimalIntrinsicSurfaceMapper &);

public:

  /// Destructor
  virtual ~NearOptimalIntrinsicSurfaceMapper();

  // ---------------------------------------------------------------------------
  // Execution

  /// Compute map values at interior points
  virtual void ComputeMap();

protected:

  /// Compute affine combination weight that minimizes a given distortion measure
  ///
  /// \param[in] u0 Discrete authalic  surface map values, i.e., lambda=0.
  /// \param[in] u1 Discrete conformal surface map values, i.e., lambda=1.
  ///
  /// \returns Affine combination weight \f$\lambda\f$, where final surface map values
  ///          are computed as \f$\lambda u1 + (1 - \lambda) * u0\f$.
  virtual double ComputeLambda(vtkDataArray *u0, vtkDataArray *u1) const = 0;

  /// Determine scale of surface mesh
  double Scale(vtkPolyData *) const;

  /// Determine scale of intrinsic parameterization
  double Scale(vtkDataArray *, vtkDataArray * = nullptr) const;

};


} // namespace mirtk

#endif // MIRTK_NearOptimalIntrinsicSurfaceMapper_H
