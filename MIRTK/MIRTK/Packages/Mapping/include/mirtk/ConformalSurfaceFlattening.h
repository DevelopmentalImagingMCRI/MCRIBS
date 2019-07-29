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

#ifndef MIRTK_ConformalSurfaceFlattening_H
#define MIRTK_ConformalSurfaceFlattening_H

#include "mirtk/SphericalSurfaceMapper.h"


namespace mirtk {


/**
 * Computes discrete conformal map of (closed) surface to the complex plane
 *
 * The conformal map is the solution of a sparse system of linear equations,
 * where the edge weights are identical to those used by Pinkall and Polthier
 * (1993) and later Eck et al. (1995). These cotangent weights are well known
 * in the finite element method (FEM) literature and correspond to a
 * discretization of the Laplace operator. The main difference of the
 * method by Angenent et al. (1999-2000) and Haker et al. (2000) compared to
 * the discrete conformal mapping of Eck et al. (1995) and Desbrun et al. (2002)
 * are the boundary constraints. While all these methods minimize the discrete
 * Dirichlet energy, the latter methods fix the boundary values of a non-closed
 * surface, whereas the method by Angenent and Haker is designed for closed
 * genus-0 surface meshes with a subsequent inverse stereographic projection
 * to the sphere. The intermediate parametric domain is a triangular subdomain
 * of the complex plane.
 *
 * \note This implementation is based on the corresponding ITK implementation
 *       originally contributed to by Gao et al. https://hdl.handle.net/1926/225
 *       which was later extended and integrated into ITK version 4 as
 *       https://itk.org/Doxygen/html/classitk_1_1ConformalFlatteningMeshFilter.html.
 *
 * - Pinkall and Polthier (1993). Computing Discrete Minimal Surfaces and Their Conjugates.
 *   Experiment. Math., 2(1), 15–36.
 * - Eck et al. (1995). Multiresolution analysis of arbitrary meshes. SIGGRAPH.
 * - Angenent et al. (1999). Conformal geometry and brain flattening. MICCAI, 271–278.
 * - Angenent et al. (n.d.). On the Laplace-Beltrami Operator and Brain Surface Flattening.
 * - Angenent et al. (2000). On Area Preserving Mappings of Minimal Distortion.
 *   In System Theory: Modeling, Analysis and Control, pp. 275–286.
 * - Haker et al. (2000). Conformal surface parameterization for texture mapping.
 *   IEEE Trans. Vis. Comput. Graphics, 6(2), 181–189.
 * - Desbrun, Meyer, and Alliez (2002). Intrinsic parameterizations of surface meshes.
 *   Computer Graphics Forum, 21(3), 209–218.
 */
class ConformalSurfaceFlattening : public SphericalSurfaceMapper
{
  mirtkObjectMacro(ConformalSurfaceFlattening);

  // ---------------------------------------------------------------------------
  // Attributes

  /// ID of cell containing polar point P used for inverse stereographic projection
  mirtkPublicAttributeMacro(int, PolarCellId);

  /// Whether to compose flattening map with inverse stereographic projection
  mirtkPublicAttributeMacro(bool, MapToSphere);

  /// Scaling factor applied to complex plane coordinates before projection to sphere
  mirtkPublicAttributeMacro(double, Scale);

  /// Radius of sphere to which flattened map is projected
  mirtkPublicAttributeMacro(double, Radius);

  /// Maximum number of iterations
  ///
  /// When the number of iterations is set to 1, a sparse direct solver is used.
  mirtkPublicAttributeMacro(int, NumberOfIterations);

  /// Tolerance for sparse linear solver
  mirtkPublicAttributeMacro(double, Tolerance);

  /// Computed map values at surface points
  mirtkAttributeMacro(vtkSmartPointer<vtkDataArray>, Values);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const ConformalSurfaceFlattening &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  ConformalSurfaceFlattening();

  /// Copy constructor
  ConformalSurfaceFlattening(const ConformalSurfaceFlattening &);

  /// Assignment operator
  ConformalSurfaceFlattening &operator =(const ConformalSurfaceFlattening &);

  /// Destructor
  virtual ~ConformalSurfaceFlattening();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Construct and solve symmetric system of linear equations
  virtual void ComputeMap();

  /// Finalize filter execution
  virtual void Finalize();

};


} // namespace mirtk

#endif // MIRTK_ConformalSurfaceFlattening_H
