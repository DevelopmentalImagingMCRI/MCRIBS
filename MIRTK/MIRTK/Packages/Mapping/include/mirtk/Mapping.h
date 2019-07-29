/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
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

#ifndef MIRTK_Mapping_H
#define MIRTK_Mapping_H

#include "mirtk/Object.h"

#include "mirtk/Point.h"
#include "mirtk/Cfstream.h"
#include "mirtk/ImageAttributes.h"

#include "vtkSmartPointer.h"
#include "vtkPointSet.h"


namespace mirtk {


// Forward declaration of possible output map type
template <class TVoxel> class GenericImage;


/**
 * Surface map or volumetric map
 *
 * A mapping assigns either each point of a surface embedded in 3D Euclidean
 * space (a 2D manifold) or each point of a 3D volume a n-D target value,
 * where n is usually 1, 2, or 3.
 *
 * It may represent an embedding of a cortical surface mesh, a bijective map of
 * a surface to a disk, square, or sphere, as well as a bijective map from one
 * piecewise linear complex (PLC) to another, such as the mapping of the brain
 * volume to the unit ball (i.e., the interior of the unit sphere). An example
 * of a surface map from computer graphics is the assignment of 2D texture
 * coordinates to every vertex of the input mesh with linear interpolation
 * within each triangle to map a 2D texture onto the surface.
 *
 * A bijective mapping is a reparameterization of the input domain, where the
 * parameterization depends on the method used to compute the respective map.
 */
class Mapping : public Object
{
  mirtkAbstractMacro(Mapping);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Value assigned to points outside the map domain
  mirtkPublicAttributeMacro(double, OutsideValue);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const Mapping &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  Mapping();

  /// Copy constructor
  Mapping(const Mapping &);

  /// Assignment operator
  Mapping &operator =(const Mapping &);

public:

  /// Read mapping from file
  static Mapping *New(const char *);

  /// Destructor
  virtual ~Mapping();

  /// Initialize map after inputs and parameters are set
  virtual void Initialize();

  /// Make deep copy of this map
  virtual Mapping *NewCopy() const = 0;

  // ---------------------------------------------------------------------------
  // Map domain

  /// Dimension of map domain
  virtual int NumberOfArguments() const;

  /// Get minimum axes-aligned bounding box of map domain
  ///
  /// \param[out] x1 Lower bound of map domain along x axis.
  /// \param[out] y1 Lower bound of map domain along y axis.
  /// \param[out] x2 Upper bound of map domain along x axis.
  /// \param[out] y2 Upper bound of map domain along y axis.
  void BoundingBox(double &x1, double &y1,
                   double &x2, double &y2) const;

  /// Get minimum axes-aligned bounding box of map domain
  ///
  /// \param[out] x1 Lower bound of map domain along x axis.
  /// \param[out] y1 Lower bound of map domain along y axis.
  /// \param[out] z1 Lower bound of map domain along z axis.
  /// \param[out] x2 Upper bound of map domain along x axis.
  /// \param[out] y2 Upper bound of map domain along y axis.
  /// \param[out] z2 Upper bound of map domain along z axis.
  virtual void BoundingBox(double &x1, double &y1, double &z1,
                           double &x2, double &y2, double &z2) const = 0;

  /// Get minimum axes-aligned bounding box of input domain
  ///
  /// \param[out] bounds Bounds of input domain in VTK order, i.e.,
  ///                    [x1, x2, y1, y2, z1, z2].
  void BoundingBox(double bounds[6]) const;

  /// Get minimum axes-aligned bounding box of input domain
  ///
  /// \param[out] p1 Lower-left-front corner of input domain bounding box.
  /// \param[out] p2 Upper-right-back corner of input domain bounding box.
  void BoundingBox(Point &p1, Point &p2) const;

  /// Get regular lattice attributes of map domain
  ///
  /// \param[in] nx Number of lattice points along x axis.
  /// \param[in] ny Number of lattice points along y axis.
  /// \param[in] nz Number of lattice points along z axis.
  ImageAttributes Attributes(int nx, int ny = 0, int nz = 0) const;

  /// Get regular lattice attributes of map domain
  ///
  /// If no lattice spacing is specified, the length of the bounding box
  /// diagonal dividied by 256 is used.
  ///
  /// \param[in] dx Lattice spacing along x axis.
  /// \param[in] dy Lattice spacing along y axis.
  /// \param[in] dz Lattice spacing along z axis.
  ImageAttributes Attributes(double dx = .0, double dy = .0, double dz = .0) const;

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Dimension of codomain, i.e., number of output values
  virtual int NumberOfComponents() const = 0;

  /// Evaluate map at a given point
  ///
  /// \param[out] v Map value.
  /// \param[in]  x Coordinate of point along x axis at which to evaluate map.
  /// \param[in]  y Coordinate of point along y axis at which to evaluate map.
  /// \param[in]  z Coordinate of point along z axis at which to evaluate map.
  ///
  /// \returns Whether input point is inside map domain.
  virtual bool Evaluate(double *v, double x, double y, double z = 0) const = 0;

  /// Evaluate map at a given point
  ///
  /// \param[out] v Map value.
  /// \param[in]  p Point at which to evaluate map.
  ///
  /// \returns Whether input point is inside map domain.
  bool Evaluate(double *v, const double p[3]) const;

  /// Evaluate map at a given point
  ///
  /// \param[out] v Map value.
  /// \param[in]  p Point at which to evaluate map.
  ///
  /// \returns Whether input point is inside map domain.
  bool Evaluate(double *v, const Point &p) const;

  /// Evaluate map at a given point
  ///
  /// \param[in] x Coordinate of point along x axis at which to evaluate map.
  /// \param[in] y Coordinate of point along y axis at which to evaluate map.
  /// \param[in] z Coordinate of point along z axis at which to evaluate map.
  /// \param[in] l Index of map value component.
  ///
  /// \returns The l-th component of the map value evaluated at the given point
  ///          or the \c OutsideValue when input point is outside the map domain.
  virtual double Evaluate(double x, double y, double z = 0, int l = 0) const;

  /// Evaluate map at a given point
  ///
  /// \param[in] p Point at which to evaluate map.
  /// \param[in] l Index of map value component.
  ///
  /// \returns The l-th component of the map value evaluated at the given point
  ///          or the \c OutsideValue when input point is outside the map domain.
  double Evaluate(const double p[3], int l = 0) const;

  /// Evaluate map at a given point
  ///
  /// \param[in] p Point at which to evaluate map.
  /// \param[in] l Index of map value component.
  ///
  /// \returns The l-th component of the map value evaluated at the given point
  ///          or the \c OutsideValue when input point is outside the map domain.
  double Evaluate(const Point &, int l = 0) const;

  /// Evaluate map at each point of a regular lattice
  ///
  /// \param[out] f Defines lattice on which to evaluate the map. The map value
  ///               at each lattice point is stored at the respective voxel.
  ///               The number of map values stored in the output image is
  ///               determined by the temporal dimension of the image.
  /// \param[in]  l Index of first map value component to store in output image.
  /// \param[in]  m Piecewise linear complex (PLC) defining an arbitrary subset
  ///               of the lattice points at which to evaluate the map.
  virtual void Evaluate(GenericImage<float> &f, int l = 0, vtkSmartPointer<vtkPointSet> m = nullptr) const;

  /// Evaluate map at each point of a regular lattice
  ///
  /// \param[out] f Defines lattice on which to evaluate the map. The map value
  ///               at each lattice point is stored at the respective voxel.
  ///               The number of map values stored in the output image is
  ///               determined by the temporal dimension of the image.
  /// \param[in]  l Index of first map value component to store in output image.
  /// \param[in]  m Piecewise linear complex (PLC) defining an arbitrary subset
  ///               of the lattice points at which to evaluate the map.
  virtual void Evaluate(GenericImage<double> &f, int l = 0, vtkSmartPointer<vtkPointSet> m = nullptr) const;

  // ---------------------------------------------------------------------------
  // I/O

  /// Read map from file
  virtual bool Read(const char *);

  /// Write map to file
  virtual bool Write(const char *) const;

protected:

  /// Read map attributes and parameters from file stream
  virtual void ReadMap(Cifstream &);

  /// Write map attributes and parameters to file stream
  virtual void WriteMap(Cofstream &) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Map domain
// =============================================================================

// -----------------------------------------------------------------------------
inline int Mapping::NumberOfArguments() const
{
  double x1, y1, z1, x2, y2, z2;
  this->BoundingBox(x1, y1, z1, x2, y2, z2);
  return (z1 == z2 ? ((y1 == y2) ? 1 : 2) : 3);
}

// -----------------------------------------------------------------------------
inline void Mapping::BoundingBox(double &x1, double &y1, double &x2, double &y2) const
{
  double z1, z2;
  this->BoundingBox(x1, y1, z1, x2, y2, z2);
}

// -----------------------------------------------------------------------------
inline void Mapping::BoundingBox(double bounds[6]) const
{
  this->BoundingBox(bounds[0], bounds[2], bounds[4],
                    bounds[1], bounds[3], bounds[5]);
}

// -----------------------------------------------------------------------------
inline void Mapping::BoundingBox(Point &p1, Point &p2) const
{
  this->BoundingBox(p1._x, p1._y, p1._z, p2._x, p2._y, p2._z);
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
inline bool Mapping::Evaluate(double *v, const double p[3]) const
{
  return this->Evaluate(v, p[0], p[1], p[2]);
}

// -----------------------------------------------------------------------------
inline bool Mapping::Evaluate(double *v, const Point &p) const
{
  return this->Evaluate(v, p._x, p._y, p._z);
}

// -----------------------------------------------------------------------------
inline double Mapping::Evaluate(double x, double y, double z, int l) const
{
  double * const v = new double[this->NumberOfComponents()];
  this->Evaluate(v, x, y, z);
  const double s = v[l];
  delete[] v;
  return s;
}

// -----------------------------------------------------------------------------
inline double Mapping::Evaluate(const double p[3], int l) const
{
  return this->Evaluate(p[0], p[1], p[2], l);
}

// -----------------------------------------------------------------------------
inline double Mapping::Evaluate(const Point &p, int l) const
{
  return this->Evaluate(p._x, p._y, p._z, l);
}


} // namespace mirtk

#endif // MIRTK_Mapping_H
