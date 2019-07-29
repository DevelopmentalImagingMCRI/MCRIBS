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

#ifndef MIRTK_MeshlessVolumeMapper_H
#define MIRTK_MeshlessVolumeMapper_H

#include "mirtk/VolumeMapper.h"

#include "mirtk/Array.h"
#include "mirtk/PointSet.h"
#include "mirtk/MeshlessMap.h"

#include "vtkAbstractCellLocator.h"


namespace mirtk {


/**
 * Base class of filters which compute a volumetric map of the interior of a
 * piecewise linear complex (PLC) using the method of fundamental solutions (MFS)
 *
 * Iteratively computes a volumetric map of the interior of the input point
 * set using the method of fundamental solutions (MFS). This implementation
 * is based on the (bi-)harmonic volumetric mapping methods presented in
 * (Li et al., 2010) and (Xu et al., 2013).
 *
 * Boundary (constraint) points and source (singularity) points are sampled
 * using the geometry adaptive sampling algorithm as outlined in Section 4.1
 * of (Li et al., 2010).
 *
 * The type of volumetric map (e.g., harmonic or biharmonic) depends on the
 * particular linear system which is defined by the subclass implementation of
 * the pure virtual base class functions GetCoefficientMatrix, AddRegularization,
 * and GetRightHandSide. This system is solved using the LU decomposition as
 * outlined in (Xu et al., 2013).
 *
 * - Li et al. (2010). Feature-aligned harmonic volumetric mapping using MFS.
 *   Computers and Graphics (Pergamon), 34(3), 242–251. doi:10.1016/j.cag.2010.03.004
 * - Xu et al. (2013). Biharmonic volumetric mapping using fundamental
 *   solutions. IEEE Transactions on Visualization and Computer Graphics,
 *   19(5), 787–798. doi:10.1109/TVCG.2012.173
 */
class MeshlessVolumeMapper : public VolumeMapper
{
  mirtkAbstractMacro(MeshlessVolumeMapper);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Ratio of number of contraints points divided by number of input surface points
  mirtkPublicAttributeMacro(double, BoundaryPointsRatio);

  /// Ratio of total number of source points divided by number of input surface points
  mirtkPublicAttributeMacro(double, SourcePointsRatio);

  /// Maximum number of source points in each subset
  mirtkPublicAttributeMacro(int, MaximumNumberOfSourcePoints);

  /// Number of iterations / approximation functions
  mirtkPublicAttributeMacro(int, NumberOfIterations);

  /// Number of lattice points for implicit surface representation
  mirtkPublicAttributeMacro(int, ImplicitSurfaceSize);

  /// Spacing of lattice points for implicit surface representation
  mirtkPublicAttributeMacro(double, ImplicitSurfaceSpacing);

  /// Distance of the offset surface from the input surface
  ///
  /// If negative, the absolute value is multiplied by the length of the
  /// diagonal of the input point set bounding box to get the absolute surface
  /// distance offset value.
  mirtkPublicAttributeMacro(double, DistanceOffset);

  /// Upper threshold of condition number of coefficient matrix
  mirtkPublicAttributeMacro(double, MaximumConditionNumber);

  /// Decimated offset surface from which to sample the source points
  mirtkReadOnlyAttributeMacro(vtkSmartPointer<vtkPolyData>, OffsetSurface);

  /// Locator used to find closest point on offset surface
  mirtkAttributeMacro(vtkSmartPointer<vtkAbstractCellLocator>, OffsetPointLocator);

  /// Partition of source points
  mirtkAttributeMacro(Array<Array<int> >, SourcePartition);

  /// Residual boundary map
  mirtkAttributeMacro(vtkSmartPointer<vtkDataArray>, ResidualMap);

protected:

  /// Get total number of boundary / constraints points
  int NumberOfBoundaryPoints() const;

  /// Get total number of source / singularity points
  int NumberOfSourcePoints() const;

  /// Get number of source points subsets
  int NumberOfSourcePointSets() const;

  /// Get number of source points in k-th subset
  int NumberOfSourcePoints(int) const;

  /// Get index of i-th point of k-th source points subset
  int SourcePointIndex(int k, int i) const;

  /// Copy attributes of this class from another instance
  void CopyAttributes(const MeshlessVolumeMapper &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  MeshlessVolumeMapper();

  /// Copy constructor
  MeshlessVolumeMapper(const MeshlessVolumeMapper &);

  /// Assignment operator
  MeshlessVolumeMapper &operator =(const MeshlessVolumeMapper &);

public:

  /// Destructor
  virtual ~MeshlessVolumeMapper();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Get closest point on offset surface
  ///
  /// \param[in] x Boundary point.
  /// \param[in] p Closest source point on offset surface.
  void GetClosestPointOnOffsetSurface(double x[3], double p[3]);

  /// Update boundary surface with corresponding boundary map as point data
  virtual void UpdateBoundary(vtkPolyData *);

  /// Sample boundary points from input surface
  virtual void PlaceBoundaryPoints();

  /// Compute and sample offset surface for placement of source points
  virtual void PlaceSourcePoints();

  /// Add new source point
  ///
  /// \param[in] q Source point coordinates.
  ///
  /// \returns Whether source point was added or too close to existing point.
  virtual bool AddSourcePoint(double q[3]);

  /// Evenly partition source points into smaller subsets
  virtual void PartitionSourcePoints();

  /// Initialize residual boundary map
  virtual void InitializeResidualMap();

  /// Update residual boundary map
  ///
  /// \returns Mean squared error of boundary map approximation.
  virtual double UpdateResidualMap(double * = nullptr,
                                   double * = nullptr,
                                   double * = nullptr);

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Compute meshless map coefficients
  virtual void Solve();

  // ---------------------------------------------------------------------------
  // Linear system

protected:

  /// Get coefficients matrix corresponding to the least squares fitting term(s)
  /// of the quadratic energy function at constraints points
  ///
  /// \param[in]  k     Index of source points subset.
  /// \param[out] coeff Coefficients matrix.
  virtual void GetCoefficients(int k, Matrix &coeff) const = 0;

  /// Get right-hand side of linear system
  ///
  /// \param[in]  k Index of source points subset.
  /// \param[out] b Right-hand side of linear system.
  virtual void GetConstraints(int k, Matrix &b) const = 0;

  /// Add solution of linear system to weights of volumetric map
  ///
  /// \param[in] k Index of source points subset.
  /// \param[in] w Solution of linear system.
  virtual void AddWeights(int k, const Matrix &w) = 0;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline int MeshlessVolumeMapper::NumberOfBoundaryPoints() const
{
  return static_cast<int>(_Boundary->GetNumberOfPoints());
}

// -----------------------------------------------------------------------------
inline int MeshlessVolumeMapper::NumberOfSourcePoints() const
{
  if (_Output) {
    MeshlessMap *map = dynamic_cast<MeshlessMap *>(_Output.get());
    return map->NumberOfSourcePoints();
  }
  return static_cast<int>(_OffsetSurface->GetNumberOfPoints());
}

// -----------------------------------------------------------------------------
inline int MeshlessVolumeMapper::NumberOfSourcePointSets() const
{
  return static_cast<int>(_SourcePartition.size());
}

// -----------------------------------------------------------------------------
inline int MeshlessVolumeMapper::NumberOfSourcePoints(int k) const
{
  return static_cast<int>(_SourcePartition[k].size());
}

// -----------------------------------------------------------------------------
inline int MeshlessVolumeMapper::SourcePointIndex(int k, int i) const
{
  return _SourcePartition[k][i];
}


} // namespace mirtk

#endif // MIRTK_MeshlessVolumeMapper_H
