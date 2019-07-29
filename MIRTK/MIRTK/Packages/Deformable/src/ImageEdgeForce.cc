/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Andreas Schuh
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

#include "mirtk/ImageEdgeForce.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/GaussianBlurring.h"
#include "mirtk/ConvolutionFunction.h"
#include "mirtk/FastLinearImageGradientFunction.h"
#include "mirtk/ObjectFactory.h"
#include "mirtk/VtkMath.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(ImageEdgeForce);


// =============================================================================
// Auxiliary functions
// =============================================================================

namespace ImageEdgeForceUtils {


// Type of edge map gradient interpolation/evaluation function
typedef GenericFastLinearImageGradientFunction<
           GenericImage<RealPixel>
        > EdgeGradient;

// -----------------------------------------------------------------------------
/// Detect image edges using the Sobel operator
void DetectEdges(RealImage &im, double sigma = .0)
{
  typedef ConvolutionFunction::ConvolveInX<RealPixel> ConvX;
  typedef ConvolutionFunction::ConvolveInY<RealPixel> ConvY;
  typedef ConvolutionFunction::ConvolveInZ<RealPixel> ConvZ;

  RealPixel h[3] = { .25, .5, .25};
  RealPixel g[3] = {-.5,  .0, .5};

  const ImageAttributes &attr = im.Attributes();
  RealImage gx(attr), gy(attr), gz(attr), gm(attr);

  if (sigma > .0) {
    GaussianBlurring<RealPixel> blur(sigma);
    blur.Input (&im);
    blur.Output(&gm);
    blur.Run();
    im.CopyFrom(gm);
  }

  ParallelForEachVoxel(ConvX(&im, g, 3, false), attr, im, gx);
  ParallelForEachVoxel(ConvY(&gx, h, 3, false), attr, gx, gm);
  ParallelForEachVoxel(ConvZ(&gm, h, 3, false), attr, gm, gx);

  ParallelForEachVoxel(ConvX(&im, h, 3, false), attr, im, gy);
  ParallelForEachVoxel(ConvY(&gy, g, 3, false), attr, gy, gm);
  ParallelForEachVoxel(ConvZ(&gm, h, 3, false), attr, gm, gy);

  ParallelForEachVoxel(ConvX(&im, h, 3, false), attr, im, gz);
  ParallelForEachVoxel(ConvY(&gz, h, 3, false), attr, gz, gm);
  ParallelForEachVoxel(ConvZ(&gm, g, 3, false), attr, gm, gz);

  for (int i = 0; i < gm.NumberOfVoxels(); ++i) {
    im(i) = sqrt(gx(i) * gx(i) + gy(i) * gy(i) + gz(i) * gz(i));
  }
}

// -----------------------------------------------------------------------------
/// Evaluate edge force
struct EvaluateGradient
{
  typedef ImageEdgeForce::GradientType Force;

  vtkPoints    *_Points;
  vtkDataArray *_Normals;
  EdgeGradient *_EdgeGradient;
  Force        *_Gradient; // Note: Opposite direction than force vector!
  double        _MaxNorm;

  EvaluateGradient() : _MaxNorm(.0) {}

  EvaluateGradient(const EvaluateGradient &other, split)
  :
    _Points(other._Points),
    _Normals(other._Normals),
    _EdgeGradient(other._EdgeGradient),
    _Gradient(other._Gradient),
    _MaxNorm(other._MaxNorm)
  {}

  void join(const EvaluateGradient &other)
  {
    if (other._MaxNorm > _MaxNorm) _MaxNorm = other._MaxNorm;
  }

  void operator ()(const blocked_range<vtkIdType> &re)
  {
    double p[3], n[3], g[3], norm;
    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      _Points->GetPoint(ptId, p);
      _EdgeGradient->WorldToImage(p[0], p[1], p[2]);
      _EdgeGradient->Evaluate(g,  p[0], p[1], p[2]);
      if (_Normals) {
        _Normals->GetTuple(ptId, n);
        norm = -vtkMath::Dot(g, n);
        vtkMath::MultiplyScalar(n, norm);
        _Gradient[ptId] = -Force(n[0], n[1], n[2]);
      } else {
//        vtkMath::Normalize(g);
        norm = vtkMath::Norm(g);
        _Gradient[ptId] = -Force(g[0], g[1], g[2]);
      }
      if (norm > _MaxNorm) _MaxNorm = norm;
    }
  }
};


} // namespace ImageEdgeForceUtils
using namespace ImageEdgeForceUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
ImageEdgeForce::ImageEdgeForce(const char *name, double weight)
:
  ExternalForce(name, weight),
  _Sigma(-1.0),
  _InNormalDirection(true)
{
  _ParameterPrefix.push_back("Image edge force ");
  _ParameterPrefix.push_back("Intensity edge force ");
  _ParameterPrefix.push_back("Edge force ");
}

// -----------------------------------------------------------------------------
void ImageEdgeForce::CopyAttributes(const ImageEdgeForce &other)
{
  _Sigma             = other._Sigma;
  _InNormalDirection = other._InNormalDirection;
  _EdgeField         = other._EdgeField;
}

// -----------------------------------------------------------------------------
ImageEdgeForce::ImageEdgeForce(const ImageEdgeForce &other)
:
  ExternalForce(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
ImageEdgeForce &ImageEdgeForce::operator =(const ImageEdgeForce &other)
{
  if (this != &other) {
    ExternalForce::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
ImageEdgeForce::~ImageEdgeForce()
{
}

// =============================================================================
// Configuration
// =============================================================================

// -----------------------------------------------------------------------------
bool ImageEdgeForce::SetWithoutPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Blurring") == 0) {
    return FromString(value, _Sigma);
  }
  if (strcmp(param, "In normal direction") == 0) {
    return FromString(value, _InNormalDirection);
  }
  return ExternalForce::SetWithoutPrefix(param, value);
}

// -----------------------------------------------------------------------------
ParameterList ImageEdgeForce::Parameter() const
{
  ParameterList params = ExternalForce::Parameter();
  InsertWithPrefix(params, "Blurring",            _Sigma);
  InsertWithPrefix(params, "In normal direction", _InNormalDirection);
  return params;
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
void ImageEdgeForce::Initialize()
{
  // Initialize base class
  ExternalForce::Initialize();
  if (_NumberOfPoints == 0) return;

  // Default smoothing sigma
  if (_Sigma < .0) {
    _Sigma = .7 * max(_Image->GetXSize(), _Image->GetYSize());
  }

  // Compute edge field
  _EdgeField = *_Image;
  DetectEdges(_EdgeField, _Sigma);
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void ImageEdgeForce::EvaluateGradient(double *gradient, double step, double weight)
{
  if (_NumberOfPoints == 0) return;

  memset(_Gradient, 0, _NumberOfPoints * sizeof(GradientType));

  EdgeGradient edge_gradient;
  edge_gradient.WrtWorld(true);
  edge_gradient.Input(&_EdgeField);
  edge_gradient.Initialize();

  ImageEdgeForceUtils::EvaluateGradient eval;
  eval._Points              = _PointSet->SurfacePoints();
  eval._Normals             = _InNormalDirection ? _PointSet->SurfaceNormals() : NULL;
  eval._EdgeGradient        = &edge_gradient;
  eval._Gradient            = _Gradient;
  parallel_reduce(blocked_range<vtkIdType>(0, _NumberOfPoints), eval);

  if (eval._MaxNorm > .0) {
    for (int i = 0; i < _NumberOfPoints; ++i) {
      _Gradient[i] /= eval._MaxNorm;
    }
  }

  ExternalForce::EvaluateGradient(gradient, step, weight / _NumberOfPoints);
}


} // namespace mirtk
