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

#ifndef MIRKT_ImageEdgeDistance_H
#define MIRKT_ImageEdgeDistance_H

#include "mirtk/SurfaceForce.h"

#include "mirtk/FastCubicBSplineInterpolateImageFunction.h"


namespace mirtk {


/**
 * External surface force which attracts the surface to nearby image edges
 */
class ImageEdgeDistance : public SurfaceForce
{
  mirtkEnergyTermMacro(ImageEdgeDistance, EM_ImageEdgeDistance);

  // ---------------------------------------------------------------------------
  // Types

public:

  /// Type of discrete intensity image
  typedef GenericImage<double> DiscreteImage;

  /// Type of image used to store local image statistics
  typedef GenericImage<float> LocalStatsImage;

  /// Type of interpolated image / image interpolation function
  typedef GenericFastCubicBSplineInterpolateImageFunction<DiscreteImage> ContinuousImage;

  /// Enumeration of edge force modes based on directional derivative of image intensities
  enum EdgeType
  {
    Extremum,             ///< Attract points to closest extrema of same sign
    ClosestMinimum,       ///< Attract points to closest minima
    ClosestMaximum,       ///< Attract points to closest maxima
    ClosestExtremum,      ///< Attract points to closest extrema
    StrongestMinimum,     ///< Attract points to strongest minima
    StrongestMaximum,     ///< Attract points to strongest maxima
    StrongestExtremum,    ///< Attract points to strongest extrema
    NeonatalWhiteSurface, ///< T2-weighted MRI  WM/cGM edge at neonatal age
    NeonatalPialSurface   ///< T2-weighted MRI cGM/CSF edge at neonatal age
  };

  // ---------------------------------------------------------------------------
  // Attributes

private:

  /// Type of edge which points are attracted to
  mirtkPublicAttributeMacro(enum EdgeType, EdgeType);

  /// Minimum foreground intensity value
  mirtkPublicAttributeMacro(double, Padding);

  /// Minimum object intensity value
  mirtkPublicAttributeMacro(double, MinIntensity);

  /// Maximum object intensity value
  mirtkPublicAttributeMacro(double, MaxIntensity);

  /// Minimum image edge gradient magnitude
  mirtkPublicAttributeMacro(double, MinGradient);

  /// Maximum image edge gradient magnitude to consider a point non-stationary
  mirtkPublicAttributeMacro(double, MaxGradient);

  /// Minimum T1-weighted image edge gradient magnitude
  ///
  /// An image edge in the T2-weighted image with negative gradient must have
  /// a corresponding positive gradient in the T1-weighted image. The minimum
  /// gradient magnitude for a T1-weighted edge to be considered is given by
  /// this parameter.
  mirtkPublicAttributeMacro(double, MinT1Gradient);

  /// Maximum T1-weighted image edge gradient magnitude
  ///
  /// An image edge in the T2-weighted image with negative gradient may not
  /// have a corresponding negative gradient in the T1-weighted image. The
  /// minimum magnitude of these negative T1-weighted intensity gradients
  /// to discard such T2-weighted edge is specified by this parameter.
  mirtkPublicAttributeMacro(double, MaxT1Gradient);

  /// Maximum edge point distance
  mirtkPublicAttributeMacro(double, MaxDistance);

  /// When positive, the force magnitude corresponds to a S-shaped map of
  /// estimated image edge distances, where points further than the specified
  /// threshold have constant force magnitude of one. Otherwise, the 95th
  /// percentile of point distances is used at each iteration.
  mirtkPublicAttributeMacro(double, DistanceThreshold);

  /// Radius of distance median filter
  mirtkPublicAttributeMacro(int, MedianFilterRadius);

  /// Number of edge distance smoothing iterations
  mirtkPublicAttributeMacro(int, DistanceSmoothing);

  /// Step length used for ray casting
  mirtkPublicAttributeMacro(double, StepLength);

  /// T1-weighted MR image
  mirtkPublicAggregateMacro(const RealImage, T1WeightedImage);

  /// White matter mask used by NeonatalWhiteSurface edge force
  mirtkPublicAggregateMacro(const BinaryImage, WhiteMatterMask);

  /// Grey matter mask used by NeonatalWhiteSurface edge force
  mirtkPublicAggregateMacro(const BinaryImage, GreyMatterMask);

  /// Approximate distance to boundary of cortex used to disable
  /// WM->dGM->cGM correction near the boundary to avoid accidental
  /// correction where a WM->cGM->BG transition is given instead
  ///
  /// \sa Application subdivide-brain-image.cc option -output-depth.
  mirtkPublicAggregateMacro(const RealImage, CorticalHullDistance);

  /// Distance to ventricles used to ignore edges caused by hyper-intense WM
  mirtkPublicAggregateMacro(const RealImage, VentriclesDistance);

  /// Distance to cerebellum used to avoid choosing second downhill part
  /// of WM->GM->BG transition over WM->GM part in case of following bright
  /// CSF near the cerebellum which is also surrounded by bright CSF
  mirtkPublicAggregateMacro(const RealImage, CerebellumDistance);

  /// Width of local white matter intensity statistics window
  mirtkPublicAttributeMacro(double, WhiteMatterWindowWidth);

  /// Width of local grey matter intensity statistics window
  mirtkPublicAttributeMacro(double, GreyMatterWindowWidth);

  /// Width of local grey matter intensity statistics window
  mirtkPublicAttributeMacro(double, T1GreyMatterWindowWidth);

  /// Global white matter intensity mean value
  mirtkAttributeMacro(double, GlobalWhiteMatterMean);

  /// Global white matter intensity variance value
  mirtkAttributeMacro(double, GlobalWhiteMatterVariance);

  /// Threshold that separates white and grey matter
  mirtkAttributeMacro(double, GlobalWhiteMatterThreshold);

  /// Global grey matter intensity mean value
  mirtkAttributeMacro(double, GlobalGreyMatterMean);

  /// Global grey matter intensity variance value
  mirtkAttributeMacro(double, GlobalGreyMatterVariance);

  /// Image with voxel-wise local white matter intensity mean values
  mirtkAttributeMacro(LocalStatsImage, LocalWhiteMatterMean);

  /// Image with voxel-wise local white matter intensity variance values
  mirtkAttributeMacro(LocalStatsImage, LocalWhiteMatterVariance);

  /// Image with voxel-wise local grey matter intensity mean values
  mirtkAttributeMacro(LocalStatsImage, LocalGreyMatterMean);

  /// Image with voxel-wise local grey matter intensity variance values
  mirtkAttributeMacro(LocalStatsImage, LocalGreyMatterVariance);

  /// Image with voxel-wise local grey matter T1 intensity mean values
  mirtkAttributeMacro(LocalStatsImage, LocalGreyMatterT1Mean);

  /// Image with voxel-wise local grey matter T1 intensity variance values
  mirtkAttributeMacro(LocalStatsImage, LocalGreyMatterT1Variance);

  /// Bounding box within which to allow WM->dGM->cGM correction
  mirtkAttributeMacro(Array<int>, CorticalDeepGreyMatterBoundingBox);

  /// Continuous T1-weighted image
  mirtkAttributeMacro(SharedPtr<ContinuousImage>, T1WeightedImageFunction);

  /// Continuous T2-weighted image
  mirtkAttributeMacro(SharedPtr<ContinuousImage>, T2WeightedImageFunction);

private:

  /// Copy attributes of this class from another instance
  void CopyAttributes(const ImageEdgeDistance &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  ImageEdgeDistance(const char * = "", double = 1.0);

  /// Copy constructor
  ImageEdgeDistance(const ImageEdgeDistance &);

  /// Assignment operator
  ImageEdgeDistance &operator =(const ImageEdgeDistance &);

  /// Destructor
  virtual ~ImageEdgeDistance();

  // ---------------------------------------------------------------------------
  // Configuration

protected:

  /// Set parameter value from string
  virtual bool SetWithoutPrefix(const char *, const char *);

public:

  // Import other overloads
  using SurfaceForce::Parameter;

  /// Get parameter name/value pairs
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Initialization

  /// Initialize external force once input and parameters have been set
  virtual void Initialize();

  // ---------------------------------------------------------------------------
  // Evaluation

public:

  /// Update moving input points and internal state of force term
  virtual void Update(bool = true);

protected:

  /// Evaluate external force term
  virtual double Evaluate();

  /// Evaluate external force
  virtual void EvaluateGradient(double *, double, double);

};

////////////////////////////////////////////////////////////////////////////////
// Enum <-> string conversion
////////////////////////////////////////////////////////////////////////////////

template <> bool FromString(const char *, enum ImageEdgeDistance::EdgeType &);
template <> string ToString(const enum ImageEdgeDistance::EdgeType &, int, char, bool);


} // namespace mirtk

#endif // MIRKT_ImageEdgeDistance_H
