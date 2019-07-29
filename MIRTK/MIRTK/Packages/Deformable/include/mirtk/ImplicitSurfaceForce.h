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

#ifndef MIRTK_ImplicitSurfaceForce_H
#define MIRTK_ImplicitSurfaceForce_H

#include "mirtk/SurfaceForce.h"

#include "mirtk/LinearInterpolateImageFunction.h"
#include "mirtk/FastLinearImageGradientFunction.h"

class vtkDataArray;


namespace mirtk {


/**
 * Base class of external forces based on a target implicit surface
 *
 * The input _Image of these force terms is the discrete distance function of
 * the implicit surface, e.g., a signed Euclidean distance transform of a binary
 * object mask.
 */
class ImplicitSurfaceForce : public SurfaceForce
{
  mirtkAbstractMacro(ImplicitSurfaceForce);

  // ---------------------------------------------------------------------------
  // Types
public:

  typedef GenericLinearInterpolateImageFunction<ImageType>  ImageFunction;
  typedef GenericFastLinearImageGradientFunction<ImageType> ImageGradient;

  /// Enumeration of implicit surface distance measures
  enum DistanceMeasureType
  {
    DM_Unknown,
    DM_Minimum,  ///< Minimum implicit surface distance
    DM_Normal    ///< Implicit surface distance in normal direction
  };

  // ---------------------------------------------------------------------------
  // Attributes
private:

  /// Type of implicit surface distance
  mirtkPublicAttributeMacro(DistanceMeasureType, DistanceMeasure);

  /// Signed distance offset
  mirtkPublicAttributeMacro(double, Offset);

  /// Minimum length of step when casting ray in normal direction
  mirtkPublicAttributeMacro(double, MinStepLength);

  /// Maximum implicit surface distance considered for ray casting
  mirtkPublicAttributeMacro(double, MaxDistance);

  /// Tolerance for implicit surface distance evaluation
  mirtkPublicAttributeMacro(double, Tolerance);

  /// Number of (normal) distance smoothing iterations
  ///
  /// \note Used only when DistanceMeasure is DM_Normal.
  mirtkPublicAttributeMacro(int, DistanceSmoothing);

  /// Whether to use heuristic to detect "holes" in implicit surface based on
  /// the normal distances and to fill in these distances by average distances
  /// of "hole" boundary points.
  ///
  /// \note Used only when DistanceMeasure is DM_Normal.
  mirtkPublicAttributeMacro(bool, FillInHoles);

  /// Continuous implicit surface distance function
  ImageFunction _Distance;

  /// Continuous implicit surface distance function gradient
  ImageGradient _DistanceGradient;

  /// Copy attributes of this class from another instance
  void CopyAttributes(const ImplicitSurfaceForce &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  ImplicitSurfaceForce(const char * = "", double = 1.0);

  /// Copy constructor
  ImplicitSurfaceForce(const ImplicitSurfaceForce &);

  /// Assignment operator
  ImplicitSurfaceForce &operator =(const ImplicitSurfaceForce &);

public:

  /// Destructor
  virtual ~ImplicitSurfaceForce();

  /// Initialize force term once input and parameters have been set
  virtual void Initialize();

  // ---------------------------------------------------------------------------
  // Configuration
protected:

  /// Set parameter value from string
  virtual bool SetWithPrefix(const char *, const char *);

  /// Set parameter value from string
  virtual bool SetWithoutPrefix(const char *, const char *);

public:

  // Import other overloads
  using SurfaceForce::Parameter;

  /// Get parameter name/value pairs
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Surface distance

  /// Get self-distance value at given world position along specified direction
  double SelfDistance(const double p[3], const double n[3]) const;

  /// Get distance value at given world position
  ///
  /// \param[in] p Point at which to evaluate implicit surface distance function.
  ///
  /// \returns Interpolated implicit surface distance function value.
  double Distance(const double p[3]) const;

  /// Get distance value at given world position along specified direction
  ///
  /// \param[in] p Starting point for ray casting.
  /// \param[in] n Direction of ray (incl. opposite direction).
  ///
  /// \returns Distance of closest intersection of ray cast from point \p p in
  ///          direction \p n (opposite directions) of length \c _MaxDistance.
  ///          If no intersection occurs, \c _MaxDistance is returned.
  double Distance(const double p[3], const double n[3]) const;

  /// Get (normalized) distance gradient at given world position
  ///
  /// \param[in]  p         Point at which to evaluate implicit surface distance gradient.
  /// \param[out] g         Gradient of implicit surface distance.
  /// \param[in]  normalize Whether to normalize the gradient vector.
  void DistanceGradient(const double p[3], double g[3], bool normalize = false) const;

protected:

  /// Get pointer to point data array of minimum implicit surface distances
  vtkDataArray *MinimumDistances() const;

  /// Initialize point data array used to store minimum implicit surface distances
  void InitializeMinimumDistances();

  /// Update minimum distances to implicit surface
  void UpdateMinimumDistances();

  /// Get pointer to point data array of implicit surface distances in normal directions
  vtkDataArray *NormalDistances() const;

  /// Initialize point data array used to store implicit surface distances in normal directions
  void InitializeNormalDistances();

  /// Update distances from implicit surface in normal direction
  void UpdateNormalDistances();

  /// Get pointer to point data array of implicit surface distances
  vtkDataArray *Distances() const;

  /// Initialize point data array used to store implicit surface distances
  void InitializeDistances();

  /// Update implicit surface distance measures
  void UpdateDistances();

};

////////////////////////////////////////////////////////////////////////////////
// ImplicitSurfaceForce::Measure from/to string conversion
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
/// Convert implicit surface distance measure enumeration value to string
template <>
inline string ToString(const enum ImplicitSurfaceForce::DistanceMeasureType &value,
                       int w, char c, bool left)
{
  const char *str;
  switch (value) {
    case ImplicitSurfaceForce::DM_Minimum: str = "Minimum"; break;
    case ImplicitSurfaceForce::DM_Normal:  str = "Normal";  break;
    default:                               str = "Unknown"; break;
  }
  return ToString(str, w, c, left);
}

// -----------------------------------------------------------------------------
/// Convert string to implicit surface distance measure enumeration value
template <>
inline bool FromString(const char *str, enum ImplicitSurfaceForce::DistanceMeasureType &value)
{
  string lstr = ToLower(str);
  if      (lstr == "minimum") value = ImplicitSurfaceForce::DM_Minimum;
  else if (lstr == "normal")  value = ImplicitSurfaceForce::DM_Normal;
  else                        value = ImplicitSurfaceForce::DM_Unknown;
  return (value != ImplicitSurfaceForce::DM_Unknown);
}


} // namespace mirtk

#endif // MIRTK_ImplicitSurfaceForce_H
