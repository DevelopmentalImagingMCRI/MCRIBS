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

#ifndef MIRTK_MeshlessMap_H
#define MIRTK_MeshlessMap_H

#include "mirtk/Mapping.h"

#include "mirtk/Math.h"
#include "mirtk/Matrix.h"
#include "mirtk/PointSet.h"


namespace mirtk {


// Forward declarations
class Cifstream;
class Cofstream;


/**
 * Meshless map defined as superposition of kernel functions
 */
class MeshlessMap : public Mapping
{
  mirtkAbstractMacro(MeshlessMap);

  // ---------------------------------------------------------------------------
  // Attributes

private:

  /// Source points, i.e., centers of kernel functions
  mirtkPublicAttributeMacro(PointSet, SourcePoints);

  /// Coefficients of the map
  ///
  /// The columns contain the source point weights for each scalar map.
  mirtkPublicAttributeMacro(Matrix, Coefficients);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const MeshlessMap &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  MeshlessMap();

  /// Copy constructor
  MeshlessMap(const MeshlessMap &);

  /// Assignment operator
  MeshlessMap &operator =(const MeshlessMap &);

  /// Initialize map after inputs and parameters are set
  virtual void Initialize();

public:

  /// Destructor
  virtual ~MeshlessMap();

  // ---------------------------------------------------------------------------
  // Map domain

  // Import other overloads
  using Mapping::BoundingBox;

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
  // Source points

  /// Add source point with zero coefficient
  ///
  /// \returns Whether source point was added or too close to existing point.
  bool AddSourcePoint(double p[3], double tol = .0);

  /// Get number of source points
  int NumberOfSourcePoints() const;

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Dimension of codomain, i.e., number of output values
  virtual int NumberOfComponents() const;

  // ---------------------------------------------------------------------------
  // I/O

protected:

  /// Read map attributes and parameters from file stream
  virtual void ReadMap(Cifstream &);

  /// Write map attributes and parameters to file stream
  virtual void WriteMap(Cofstream &) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline bool MeshlessMap::AddSourcePoint(double p[3], double tol)
{
  if (tol > .0) {
    for (int i = 0; i < _SourcePoints.Size(); ++i) {
      const Point &q = _SourcePoints(i);
      if (fequal(p[0], q._x, tol) &&
          fequal(p[1], q._y, tol) &&
          fequal(p[2], q._z, tol)) {
        return false;
      }
    }
  }
  _SourcePoints.Add(p);
  _Coefficients.Resize(_SourcePoints.Size(), _Coefficients.Cols());
  return true;
}

// -----------------------------------------------------------------------------
inline int MeshlessMap::NumberOfSourcePoints() const
{
  return _SourcePoints.Size();
}

// -----------------------------------------------------------------------------
inline int MeshlessMap::NumberOfComponents() const
{
  return _Coefficients.Cols();
}


} // namespace mirtk

#endif // MIRTK_MeshlessMap_H
