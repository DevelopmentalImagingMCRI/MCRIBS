/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2015-2016 Imperial College London
 * Copyright 2015-2016 Andreas Schuh
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

#ifndef MIRTK_AsConformalAsPossibleMapper_H
#define MIRTK_AsConformalAsPossibleMapper_H

#include "mirtk/LinearTetrahedralMeshMapper.h"

#include "mirtk/Array.h"
#include "mirtk/Matrix3x3.h"


namespace mirtk {


/**
 * Compute as-conformal-as-possible (ACAP) piecewise linear volumetric map
 *
 * Paillé & Poulin (2012), As-conformal-as-possible discrete volumetric mapping,
 * Computers and Graphics, 36(5), 427–433.
 */
class AsConformalAsPossibleMapper : public LinearTetrahedralMeshMapper
{
  mirtkObjectMacro(AsConformalAsPossibleMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Uniform weight of scale and angle conformality
  mirtkPublicAttributeMacro(double, UniformWeight);

  /// Local orientation of tetrahedron (rotation matrix)
  mirtkAttributeMacro(Array<Matrix3x3>, Orientation);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const AsConformalAsPossibleMapper &);

public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  AsConformalAsPossibleMapper();

  /// Copy constructor
  AsConformalAsPossibleMapper(const AsConformalAsPossibleMapper &);

  /// Assignment operator
  AsConformalAsPossibleMapper &operator =(const AsConformalAsPossibleMapper &);

  /// Destructor
  virtual ~AsConformalAsPossibleMapper();

  // ---------------------------------------------------------------------------
  // Exection

protected:

  /// Initialize filter after input and parameters are set
  void Initialize();

  /// Finalize filter execution
  void Finalize();

  // ---------------------------------------------------------------------------
  // Auxiliary functions

  /// Calculate operator weight for given tetrahadron
  ///
  /// \param[in] cellId ID of tetrahedron.
  /// \param[in] v0     First  vertex/point of tetrahedron.
  /// \param[in] v1     Second vertex/point of tetrahedron.
  /// \param[in] v2     Third  vertex/point of tetrahedron.
  /// \param[in] v3     Fourth vertex/point of tetrahedron.
  /// \param[in] volume Volume of tetrahedron.
  ///
  /// \return Operator weight contribution of tetrahedron.
  virtual Matrix3x3 GetWeight(vtkIdType cellId,
                              const double v0[3],
                              const double v1[3],
                              const double v2[3],
                              const double v3[3],
                              double       volume) const;

};


} // namespace mirtk

#endif // MIRTK_AsConformalAsPossibleMapper_H
