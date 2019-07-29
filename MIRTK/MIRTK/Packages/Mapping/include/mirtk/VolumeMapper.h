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

#ifndef MIRTK_VolumeMapper_H
#define MIRTK_VolumeMapper_H

#include "mirtk/Object.h"

#include "mirtk/Mapping.h"

#include "vtkSmartPointer.h"
#include "vtkPointSet.h"
#include "vtkPolyData.h"
#include "vtkDataArray.h"


namespace mirtk {


/**
 * Base class of volumetric map solvers
 *
 * Subclasses implement solvers for the computation of a volumetric map given
 * map values on the boundary surface of the volume as boundary conditions.
 * Whether additional interior constraints are allowed depends on the respective
 * solver implementation.
 */
class VolumeMapper : public Object
{
  mirtkAbstractMacro(VolumeMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Input point set (e.g., surface mesh or tetrahedral mesh)
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPointSet>, InputSet);

  /// Input boundary map
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkDataArray>, InputMap);

  /// Boundary surface of input point set
  mirtkAttributeMacro(vtkSmartPointer<vtkPolyData>, Boundary);

  /// Boundary surface map
  mirtkAttributeMacro(vtkSmartPointer<vtkDataArray>, BoundaryMap);

  /// Volumetric map
  ///
  /// \note The output map is uninitialized! Mapping::Initialize must be
  ///       called before this map can be evaluated at map domain points.
  mirtkReadOnlyAttributeMacro(SharedPtr<Mapping>, Output);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const VolumeMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  VolumeMapper();

  /// Copy constructor
  VolumeMapper(const VolumeMapper &);

  /// Assignment operator
  VolumeMapper &operator =(const VolumeMapper &);

public:

  /// Destructor
  virtual ~VolumeMapper();

  // ---------------------------------------------------------------------------
  // Auxiliaries

  /// Dimension of codomain of volumetric map
  int NumberOfComponents() const;

  // ---------------------------------------------------------------------------
  // Execution

  /// Parameterize interior of input data set
  void Run();

protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Initialize boundary surface with corresponding boundary map as point data
  virtual void InitializeBoundary(vtkPointSet *, vtkDataArray *);

  /// Compute map value at free points
  virtual void Solve() = 0;

  /// Finalize filter execution
  virtual void Finalize();

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline int VolumeMapper::NumberOfComponents() const
{
  return _InputMap ? static_cast<int>(_InputMap->GetNumberOfComponents()) : 0;
}


} // namespace mirtk

#endif // MIRTK_VolumeMapper_H
