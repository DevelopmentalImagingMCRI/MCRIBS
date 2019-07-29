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

#ifndef MIRTK_PiecewiseLinearMap_H
#define MIRTK_PiecewiseLinearMap_H

#include "mirtk/Mapping.h"
#include "mirtk/Point.h"

#include "vtkSmartPointer.h"
#include "vtkDataSet.h"
#include "vtkDataArray.h"
#include "vtkAbstractCellLocator.h"


namespace mirtk {


/**
 * Piecewise linear map defined at mesh nodes
 *
 * This map is defined by its values at discrete mesh nodes.
 * Intermediate values are interpolated using the weights of the respective
 * mesh cell which this point belongs to. If the point is not contained in
 * any mesh cell, the point is mapped to a constant outside value.
 *
 * A surface map is commonly represented by the texture coordinates of a
 * triangular surface mesh, while a tetrahedral mesh is commonly used to
 * parameterize a volumetric map. These maps are usually computed using a
 * finite element method (FEM).
 */
class PiecewiseLinearMap : public Mapping
{
  mirtkObjectMacro(PiecewiseLinearMap);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Mesh which discretizes the domain of this piecewise linear map
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkDataSet>, Domain);

  /// Point data array with map values at mesh points
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkDataArray>, Values);

  /// Locates cell within which a given point lies
  mirtkAttributeMacro(vtkSmartPointer<vtkAbstractCellLocator>, Locator);

  /// Maximum number cell points
  mirtkAttributeMacro(int, MaxCellSize);

  /// Squared distance tolerance used to locate cells
  /// \note Unused argument of vtkCellLocator::FindCell (as of VTK <= 7.0).
  static const double _Tolerance2;

  /// Copy attributes of this class from another instance
  void CopyAttributes(const PiecewiseLinearMap &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  PiecewiseLinearMap();

  /// Copy constructor
  PiecewiseLinearMap(const PiecewiseLinearMap &);

  /// Assignment operator
  PiecewiseLinearMap &operator =(const PiecewiseLinearMap &);

  /// Initialize map after inputs and parameters are set
  virtual void Initialize();

  /// Make deep copy of this volumetric map
  virtual Mapping *NewCopy() const;

  /// Destructor
  virtual ~PiecewiseLinearMap();

  // ---------------------------------------------------------------------------
  // Map domain

  // Import other overloads
  using Mapping::BoundingBox;

  /// Number of discrete points at which map values are given
  int NumberOfPoints() const;

  /// Get i-th domain mesh point
  ///
  /// \param[in] i Point index.
  ///
  /// \returns Point coordinates.
  class Point Point(int i) const;

  /// Get i-th domain mesh point
  ///
  /// \param[in]  i Point index.
  /// \param[out] x Point coordinate along x axis.
  /// \param[out] y Point coordinate along y axis.
  void GetPoint(int i, double &x, double &y) const;

  /// Get i-th domain mesh point
  ///
  /// \param[in]  i Point index.
  /// \param[out] x Point coordinate along x axis.
  /// \param[out] y Point coordinate along y axis.
  /// \param[out] z Point coordinate along z axis.
  void GetPoint(int i, double &x, double &y, double &z) const;

  /// Get i-th domain mesh point
  ///
  /// \param[in]  i Point index.
  /// \param[out] p Point coordinates.
  void GetPoint(int i, double p[3]) const;

  /// Get minimum axes-aligned bounding box of map domain
  ///
  /// \param[out] x1 Lower bound of map domain along x axis.
  /// \param[out] y1 Lower bound of map domain along y axis.
  /// \param[out] z1 Lower bound of map domain along z axis.
  /// \param[out] x2 Upper bound of map domain along x axis.
  /// \param[out] y2 Upper bound of map domain along y axis.
  /// \param[out] z2 Upper bound of map domain along z axis.
  virtual void BoundingBox(double &x1, double &y1, double &z1,
                           double &x2, double &y2, double &z2) const;

  // ---------------------------------------------------------------------------
  // Map codomain

  /// Dimension of codomain, i.e., number of output values
  virtual int NumberOfComponents() const;

  /// Mesh discretizing the map domain with mapped mesh points
  ///
  /// \note Use this function only when the dimension of the codomain is 2 or 3.
  ///
  /// \sa NumberOfComponents.
  ///
  /// \returns Mesh discretizing the codomain of the map, where the mesh points
  ///          correspond to the unmapped points of the map domain mesh or nullptr
  ///          when the dimension of the codomain is not 2 or 3.
  vtkSmartPointer<vtkDataSet> Codomain() const;

  // ---------------------------------------------------------------------------
  // Evaluation

  // Import other overloads
  using Mapping::Evaluate;

  /// Evaluate map at a given mesh point
  ///
  /// \param[in]  i Point index.
  /// \param[out] v Map value.
  void GetValue(int i, double *v) const;

  /// Evaluate map at a given mesh point
  ///
  /// \param[in] i Point index.
  /// \param[in] l Index of map value component.
  ///
  /// \returns The l-th component of the map value at the given mesh point.
  double Value(int i, int l = 0) const;

  /// Evaluate map at a given point
  ///
  /// \param[out] v Map value.
  /// \param[in]  x Coordinate of point along x axis at which to evaluate map.
  /// \param[in]  y Coordinate of point along y axis at which to evaluate map.
  /// \param[in]  z Coordinate of point along z axis at which to evaluate map.
  ///
  /// \returns Whether input point is inside map domain.
  virtual bool Evaluate(double *v, double x, double y, double z = 0) const;

  /// Evaluate map at a given point
  ///
  /// \param[in] x Coordinate of point along x axis at which to evaluate map.
  /// \param[in] y Coordinate of point along y axis at which to evaluate map.
  /// \param[in] z Coordinate of point along z axis at which to evaluate map.
  /// \param[in] l Index of map value component.
  ///
  /// \returns The l-th component of the map value evaluate at the given point
  ///          or the \c OutsideValue when input point is outside the map domain.
  virtual double Evaluate(double x, double y, double z = 0, int l = 0) const;

  // ---------------------------------------------------------------------------
  // I/O

  /// Read map from file
  virtual bool Read(const char *);

  /// Write map to file
  virtual bool Write(const char *) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline int PiecewiseLinearMap::NumberOfPoints() const
{
  return static_cast<int>(_Domain->GetNumberOfPoints());
}

// -----------------------------------------------------------------------------
inline void PiecewiseLinearMap::GetPoint(int i, double p[3]) const
{
  _Domain->GetPoint(static_cast<vtkIdType>(i), p);
}

// -----------------------------------------------------------------------------
inline void PiecewiseLinearMap::GetPoint(int i, double &x, double &y) const
{
  double p[3];
  _Domain->GetPoint(static_cast<vtkIdType>(i), p);
  x = p[0], y = p[1];
}

// -----------------------------------------------------------------------------
inline void PiecewiseLinearMap::GetPoint(int i, double &x, double &y, double &z) const
{
  double p[3];
  _Domain->GetPoint(static_cast<vtkIdType>(i), p);
  x = p[0], y = p[1], z = p[2];
}

// -----------------------------------------------------------------------------
inline Point PiecewiseLinearMap::Point(int i) const
{
  double p[3];
  this->GetPoint(i, p);
  return p;
}

// -----------------------------------------------------------------------------
inline void PiecewiseLinearMap::GetValue(int i, double *v) const
{
  _Values->GetTuple(static_cast<vtkIdType>(i), v);
}

// -----------------------------------------------------------------------------
inline double PiecewiseLinearMap::Value(int i, int l) const
{
  return _Values->GetComponent(static_cast<vtkIdType>(i), l);
}


} // namespace mirtk

#endif // MIRTK_PiecewiseLinearMap_H
