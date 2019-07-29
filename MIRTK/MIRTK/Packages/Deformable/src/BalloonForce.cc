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

#include "mirtk/BalloonForce.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/LinearInterpolateImageFunction.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/ObjectFactory.h"
#include "mirtk/VtkMath.h"

#include "vtkImageData.h"
#include "vtkImageStencilData.h"
#include "vtkImageStencilIterator.h"


namespace mirtk {


// See mirtk/Options.h
MIRTK_Common_EXPORT extern int debug;

// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(BalloonForce);


// =============================================================================
// Auxiliary functions
// =============================================================================

namespace BalloonForceUtils {

typedef GenericLinearInterpolateImageFunction<
           GenericImage<RegisteredImage::VoxelType>
         > ImageFunction;

// -----------------------------------------------------------------------------
/// Convert MIRTK image to VTK image
///
/// Note that vtkImageData has no implicit orientation. Therefore we just
/// convert the image in voxel coordinates (origin at 0 and voxel size 1x1x1).
vtkSmartPointer<vtkImageData> ConvertImage(const RegisteredImage *image)
{
  vtkSmartPointer<vtkImageData> imagedata = vtkSmartPointer<vtkImageData>::New();
  imagedata->SetOrigin(.0, .0, .0);
  imagedata->SetDimensions(image->X(), image->Y(), image->Z());
  imagedata->SetSpacing(1.0, 1.0, 1.0);
#if VTK_MAJOR_VERSION >= 6
  imagedata->AllocateScalars(VTK_FLOAT, 1);
#else
  imagedata->SetScalarType(VTK_FLOAT);
  imagedata->AllocateScalars();
#endif
  const int nvox = image->NumberOfSpatialVoxels();
  const RegisteredImage::VoxelType *ptr1 = image->Data();
  float *ptr2 = reinterpret_cast<float *>(imagedata->GetScalarPointer());
  for (int i = 0; i < nvox; ++i, ++ptr1, ++ptr2) {
    *ptr2 = static_cast<float>(*ptr1);
  }
  return imagedata;
}

// -----------------------------------------------------------------------------
/// Convert surface image stencil to binary mask, e.g., for debugging purposes
BinaryImage ImageStencilToMask(const ImageAttributes &attr,
                               vtkImageData          *data,
                               vtkImageStencilData   *stencil)
{
  BinaryImage mask(attr);
  const float * const start = reinterpret_cast<float *>(data->GetScalarPointer());
  vtkImageStencilIterator<float> it;
  it.Initialize(data, stencil, data->GetExtent());
  while (!it.IsAtEnd()) {
    if (it.IsInStencil()) {
      for (const float *cur = it.BeginSpan(); cur != it.EndSpan(); ++cur) {
        mask(cur - start) = true;
      }
    }
    it.NextSpan();
  }
  return mask;
}

// -----------------------------------------------------------------------------
/// Convert surface image stencil to binary mask, e.g., for debugging purposes
BinaryImage ImageStencilToMask(const ImageAttributes &attr,
                               vtkImageData          *data,
                               vtkImageStencilData   *stencil,
                               double                 p[3],
                               double                 radius)
{
  BinaryImage mask(attr);

  double rx = radius / attr._dx;
  double ry = radius / attr._dy;
  double rz = radius / attr._dz;

  double x = p[0], y = p[1], z = p[2];
  attr.WorldToLattice(x, y, z);

  int extent[6];
  extent[0] = ifloor(x - rx);
  extent[1] = iceil (x + rx);
  extent[2] = ifloor(y - ry);
  extent[3] = iceil (y + ry);
  extent[4] = ifloor(z - rz);
  extent[5] = iceil (z + rz);

  int i, i2, j, k, iter;
  for (k = extent[4]; k <= extent[5]; ++k)
  for (j = extent[2]; j <= extent[3]; ++j) {
    iter = 0;
    while (stencil->GetNextExtent(i, i2, extent[0], extent[1], j, k, iter)) {
      for (; i <= i2; ++i) {
        mask(i, j, k) = true;
      }
    }
  }

  return mask;
}

// -----------------------------------------------------------------------------
/// Compute point intensity thresholds based on local image statistics
struct ComputeLocalIntensityThresholds
{
  vtkPoints    *_Points;
  vtkDataArray *_Status;
  BaseImage    *_Image;
  BinaryImage  *_Mask;
  vtkDataArray *_LowerIntensity;
  vtkDataArray *_UpperIntensity;
  double        _LowerSigma;
  double        _UpperSigma;
  double        _RadiusX;
  double        _RadiusY;
  double        _RadiusZ;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    int    num, i, j, k, i1, i2, j1, j2, k1, k2;
    double p[3], value, mu, var, delta, sigma;
    vtkSmartPointer<vtkImageStencilData> roi = vtkSmartPointer<vtkImageStencilData>::New();
    vtkImageStencilIterator<float> it;

    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      mu = var = 0., num = 0;
      if (_Status == nullptr || _Status->GetComponent(ptId, 0) != 0.) {
        _Points->GetPoint(ptId, p);
        _Image->WorldToImage(p[0], p[1], p[2]);

        i1 = max(ifloor(p[0] - _RadiusX), 0);
        i2 = min(iceil (p[0] + _RadiusX), _Image->X()-1);
        j1 = max(ifloor(p[1] - _RadiusY), 0);
        j2 = min(iceil (p[1] + _RadiusY), _Image->Y()-1);
        k1 = max(ifloor(p[2] - _RadiusZ), 0);
        k2 = min(iceil (p[2] + _RadiusZ), _Image->Z()-1);

        for (k = k1; k <= k2; ++k)
        for (j = j1; j <= j2; ++j)
        for (i = i1; i <= i2; ++i) {
          if (_Mask->Get(i, j, k) != 0) {
            ++num;
            value = _Image->GetAsDouble(i, j, k);
            delta = value - mu;
            mu  += delta / num;
            var += delta * (value - mu);
          }
        }
        if (num > 2) var /= num - 1;
        else         var  = 0.;
        sigma = sqrt(var);
      }
      if (_LowerIntensity) {
        if (num == 0) {
          _LowerIntensity->SetComponent(ptId, 0, -inf);
        } else {
          _LowerIntensity->SetComponent(ptId, 0, mu - _LowerSigma * sigma);
        }
      }
      if (_UpperIntensity) {
        if (num == 0) {
          _LowerIntensity->SetComponent(ptId, 0, +inf);
        } else {
          _UpperIntensity->SetComponent(ptId, 0, mu + _UpperSigma * sigma);
        }
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Compute local statistics of intensities inside/outside surface mesh
struct ComputeLocalIntensityStatistics
{
  vtkPoints         *_Points;
  vtkDataArray      *_Status;
  const BaseImage   *_Image;
  const BinaryImage *_ForegroundMask;
  const BinaryImage *_BackgroundMask;
  vtkDataArray      *_ForegroundStatistics;
  vtkDataArray      *_BackgroundStatistics;
  double             _RadiusX;
  double             _RadiusY;
  double             _RadiusZ;

  enum Label { Outside = -1, BG = 0, FG = 1 };

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    int    num[2], i, j, k, i1, i2, j1, j2, k1, k2;
    double p[3], value, delta, mean[2], var[2];
    Label  c;

    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      mean[BG] = mean[FG] = var[BG] = var[FG] = .0;
      if (_Status == nullptr || _Status->GetComponent(ptId, 0) != 0.) {
        _Points->GetPoint(ptId, p);
        _Image->WorldToImage(p[0], p[1], p[2]);

        i1 = max(ifloor(p[0] - _RadiusX), 0);
        i2 = min(iceil (p[0] + _RadiusX), _Image->X()-1);
        j1 = max(ifloor(p[1] - _RadiusY), 0);
        j2 = min(iceil (p[1] + _RadiusY), _Image->Y()-1);
        k1 = max(ifloor(p[2] - _RadiusZ), 0);
        k2 = min(iceil (p[2] + _RadiusZ), _Image->Z()-1);

        num[BG] = num[FG] = 0;
        for (k = k1; k <= k2; ++k)
        for (j = j1; j <= j2; ++j)
        for (i = i1; i <= i2; ++i) {
          if (_ForegroundMask->Get(i, j, k) != 0) {
            c = FG;
          } else if (_BackgroundMask) {
            c = (_BackgroundMask->Get(i, j, k) != 0 ? BG : Outside);
          } else {
            c = (_Image->IsForeground(i, j, k) ? BG : Outside);
          }
          if (c != Outside) {
            value = _Image->GetAsDouble(i, j, k);
            delta = value - mean[c];
            num [c] += 1;
            mean[c] += delta / num[c];
            var [c] += delta * (value - mean[c]);
          }
        }
        if (num[BG] > 2) var[BG] /= num[BG] - 1;
        else             var[BG]  = 0.;
        if (num[BG] > 2) var[FG] /= num[FG] - 1;
        else             var[FG]  = 0.;
      }
      _BackgroundStatistics->SetComponent(ptId, 0, mean[BG]);
      _BackgroundStatistics->SetComponent(ptId, 1, sqrt(var[BG]));
      _ForegroundStatistics->SetComponent(ptId, 0, mean[FG]);
      _ForegroundStatistics->SetComponent(ptId, 1, sqrt(var[FG]));
    }
  }
};

// -----------------------------------------------------------------------------
/// Update balloon force magnitude and sign
struct UpdateMagnitude
{
  vtkPoints           *_Points;
  vtkDataArray        *_Status;
  const ImageFunction *_Image;
  bool                 _DeflateSurface;
  vtkDataArray        *_Intensity;
  double               _LowerIntensity;
  double               _UpperIntensity;
  vtkDataArray        *_LocalLowerIntensity;
  vtkDataArray        *_LocalUpperIntensity;
  vtkDataArray        *_BackgroundStatistics;
  double               _BackgroundSigmaFactor;
  vtkDataArray        *_ForegroundStatistics;
  double               _ForegroundSigmaFactor;
  vtkDataArray        *_Magnitude;
  double               _MagnitudeDamping;
  double               _MagnitudeThreshold;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    bool   inside;
    double m, p[3], v, mean, sigma, bgPb, fgPb;

    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      if (_Status && _Status->GetComponent(ptId, 0) == .0) {
        _Magnitude->SetComponent(ptId, 0, 0.);
        continue;
      }
      m = _Magnitude->GetComponent(ptId, 0);
      if (m == .0) continue;
      // Get intensity at current node position
      _Points->GetPoint(ptId, p);
      _Image->WorldToImage(p[0], p[1], p[2]);
      if (_Image->IsForeground(p[0], p[1], p[2])) {
        v = _Image->Evaluate(p[0], p[1], p[2]);
        if (_Intensity) _Intensity->SetComponent(ptId, 0, v);
        // Check global intensity thresholds
        inside = (_LowerIntensity <= v && v <= _UpperIntensity);
        // Check local intensity thresholds
        if (inside) {
          if (_LocalLowerIntensity && v < _LocalLowerIntensity->GetComponent(ptId, 0)) {
            inside = false;
          }
          if (_LocalUpperIntensity && v > _LocalUpperIntensity->GetComponent(ptId, 0)) {
            inside = false;
          }
        }
        // Check background/foreground statistics
        if (inside && _BackgroundStatistics && _ForegroundStatistics) {
          mean = _BackgroundStatistics->GetComponent(ptId, 0);
          if (mean == .0 || v < mean) {
            inside = false;
          } else {
            sigma = _BackgroundStatistics->GetComponent(ptId, 1) * _BackgroundSigmaFactor;
            if (sigma > .0) bgPb = exp(- pow(v - mean, 2) / (2.0 * sigma * sigma)) / sigma;
            else            bgPb = .0;
            mean  = _ForegroundStatistics->GetComponent(ptId, 0);
            sigma = _ForegroundStatistics->GetComponent(ptId, 1) * _ForegroundSigmaFactor;
            if (sigma > .0) {
              fgPb = exp(- pow(v - mean, 2) / (2.0 * sigma * sigma)) / sigma;
              if (bgPb > fgPb) {
                inside = false;
              }
            } else if (v != mean) {
              inside = false;
            }
          }
        }
        // Adjust sign of balloon force and damp magnitude if direction changes
        if (_DeflateSurface) {
          if (inside) {
            if (m < .0) m = - _MagnitudeDamping * m;
          } else {
            if (m > .0) m = - _MagnitudeDamping * m;
          }
        } else {
          if (inside) {
            if (m < .0) m = - _MagnitudeDamping * m;
          } else {
            if (m > .0) m = - _MagnitudeDamping * m;
          }
        }
      // Zero force outside image foreground
      } else {
        if (m > .0) m = - _MagnitudeDamping * m;
      }
      // Set new force magnitude
      if (abs(m) < _MagnitudeThreshold) m = .0;
      _Magnitude->SetComponent(ptId, 0, m);
    }
  }
};

// -----------------------------------------------------------------------------
/// Average balloon force magnitude across adjacent nodes
struct SmoothMagnitude
{
  const EdgeTable *_EdgeTable;
  vtkDataArray    *_Normals;
  vtkDataArray    *_Input;
  vtkDataArray    *_Output;
  double           _MinMagnitude;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    double     v, m, n0[3], n1[3], w, W;
    int        numAdjPts;
    const int *adjPts;

    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      m = abs(_Input->GetComponent(ptId, 0));
//      if (m != .0) {
        W = 1.0;
        _Normals->GetTuple(ptId, n0);
        _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPts);
        for (int i = 0; i < numAdjPts; ++i) {
          v = abs(_Input->GetComponent(adjPts[i], 0));
          _Normals->GetTuple(adjPts[i], n1);
          w = vtkMath::Dot(n0, n1);
          if (w > .5) {
            m += w * v;
            W += w;
          }
        }
        m /= W;
        m = copysign(m, _Input->GetComponent(ptId, 0));
//        m = copysign(max(m, _MinMagnitude), _Input->GetComponent(ptId, 0));
//      }
      _Output->SetComponent(ptId, 0, m);
    }
  }
};

// -----------------------------------------------------------------------------
/// Compute balloon force gradient (i.e., negative force)
struct ComputeGradient
{
  typedef BalloonForce::GradientType GradientType;

  vtkDataArray *_Normals;
  vtkDataArray *_Magnitude;
  GradientType *_Gradient;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    double n[3];
    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      _Normals->GetTuple(ptId, n);
      _Gradient[ptId] = -_Magnitude->GetComponent(ptId, 0) * GradientType(n);
    }
  }
};


} // namespace BalloonForceUtils
using namespace BalloonForceUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
BalloonForce::BalloonForce(const char *name, double weight)
:
  SurfaceForce(name, weight),
  _ForegroundMask(nullptr),
  _DeflateSurface(false),
  _LowerIntensity(-inf),
  _UpperIntensity(+inf),
  _LowerIntensitySigma(5.),
  _UpperIntensitySigma(5.),
  _ForegroundSigmaFactor(-1.),
  _BackgroundSigmaFactor(-1.),
  _Radius(-1.),
  _DampingFactor(.67),
  _MagnitudeThreshold(.1)
{
}

// -----------------------------------------------------------------------------
void BalloonForce::CopyAttributes(const BalloonForce &other)
{
  _ForegroundMask        = other._ForegroundMask;
  _DeflateSurface        = other._DeflateSurface;
  _LowerIntensity        = other._LowerIntensity;
  _UpperIntensity        = other._UpperIntensity;
  _LowerIntensitySigma   = other._LowerIntensitySigma;
  _UpperIntensitySigma   = other._UpperIntensitySigma;
  _ForegroundSigmaFactor = other._ForegroundSigmaFactor;
  _BackgroundSigmaFactor = other._BackgroundSigmaFactor;
  _Radius                = other._Radius;
  _DampingFactor         = other._DampingFactor;
  _MagnitudeThreshold    = other._MagnitudeThreshold;
}

// -----------------------------------------------------------------------------
BalloonForce::BalloonForce(const BalloonForce &other)
:
  SurfaceForce(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
BalloonForce &BalloonForce::operator =(const BalloonForce &other)
{
  if (this != &other) {
    SurfaceForce::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
BalloonForce::~BalloonForce()
{
}

// =============================================================================
// Configuration
// =============================================================================

// -----------------------------------------------------------------------------
bool BalloonForce::SetWithoutPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Deflate")         == 0 ||
      strcmp(param, "Deflate surface") == 0) {
    return FromString(value, _DeflateSurface);
  }
  if (strcmp(param, "Minimum intensity")         == 0 ||
      strcmp(param, "Lower threshold")           == 0 ||
      strcmp(param, "Lower intensity")           == 0 ||
      strcmp(param, "Lower intensity threshold") == 0) {
    return FromString(value, _LowerIntensity);
  }
  if (strcmp(param, "Maximum intensity")         == 0 ||
      strcmp(param, "Upper threshold")           == 0 ||
      strcmp(param, "Upper intensity")           == 0 ||
      strcmp(param, "Upper intensity threshold") == 0) {
    return FromString(value, _UpperIntensity);
  }
  if (strcmp(param, "Local intensity sigma")        == 0 ||
      strcmp(param, "Local intensity sigma factor") == 0) {
    double sigma;
    if (!FromString(value, sigma)) return false;
    _LowerIntensitySigma = _UpperIntensitySigma = sigma;
    return true;
  }
  if (strcmp(param, "Lower local intensity sigma")        == 0 ||
      strcmp(param, "Lower local intensity sigma factor") == 0) {
    return FromString(value, _LowerIntensitySigma);
  }
  if (strcmp(param, "Upper local intensity sigma")        == 0 ||
      strcmp(param, "Upper local intensity sigma factor") == 0) {
    return FromString(value, _UpperIntensitySigma);
  }
  if (strcmp(param, "Background intensity sigma")        == 0 ||
      strcmp(param, "Background intensity sigma factor") == 0) {
    return FromString(value, _BackgroundSigmaFactor);
  }
  if (strcmp(param, "Foreground intensity sigma")        == 0 ||
      strcmp(param, "Foreground intensity sigma factor") == 0) {
    return FromString(value, _ForegroundSigmaFactor);
  }
  if (strcmp(param, "Local window radius") == 0) {
    return FromString(value, _Radius);
  }
  if (strcmp(param, "Damping factor")           == 0 ||
      strcmp(param, "Magnitude damping factor") == 0 ||
      strcmp(param, "Magnitude damping")        == 0) {
    return FromString(value, _DampingFactor);
  }
  if (strcmp(param, "Magnitude threshold") == 0) {
    return FromString(value, _MagnitudeThreshold);
  }
  return SurfaceForce::SetWithoutPrefix(param, value);
}

// -----------------------------------------------------------------------------
ParameterList BalloonForce::Parameter() const
{
  ParameterList params = SurfaceForce::Parameter();
  InsertWithPrefix(params, "Deflate surface",             _DeflateSurface);
  InsertWithPrefix(params, "Lower intensity",             _LowerIntensity);
  InsertWithPrefix(params, "Upper intensity",             _UpperIntensity);
  InsertWithPrefix(params, "Lower local intensity sigma", _LowerIntensitySigma);
  InsertWithPrefix(params, "Upper local intensity sigma", _UpperIntensitySigma);
  InsertWithPrefix(params, "Background intensity sigma",  _BackgroundSigmaFactor);
  InsertWithPrefix(params, "Foreground intensity sigma",  _ForegroundSigmaFactor);
  InsertWithPrefix(params, "Local window radius",         _Radius);
  InsertWithPrefix(params, "Damping factor",              _DampingFactor);
  InsertWithPrefix(params, "Magnitude threshold",         _MagnitudeThreshold);
  return params;
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
void BalloonForce::ComputeLocalIntensityAttributes(bool thresholds, bool bg_fg_stats)
{
  if (!thresholds && !bg_fg_stats) return;

  BinaryImage fg_mask;
  if (_ForegroundMask == nullptr) {
    vtkSmartPointer<vtkImageData>        imagedata = ConvertImage(_Image);
    vtkSmartPointer<vtkPointSet>         surface   = WorldToImage(_PointSet->Surface(), _Image);
    vtkSmartPointer<vtkImageStencilData> stencil   = ImageStencil(imagedata, surface);
    fg_mask = ImageStencilToMask(_Image->Attributes(), imagedata, stencil);
  }

  if (thresholds) {
    const bool optional = true;
    ComputeLocalIntensityThresholds eval;
    eval._Points         = Points();
    eval._Status         = Status();
    eval._Image          = _Image;
    eval._Mask           = (_ForegroundMask ? _ForegroundMask : &fg_mask);
    eval._RadiusX        = _Radius / _Image->XSize();
    eval._RadiusY        = _Radius / _Image->YSize();
    eval._RadiusZ        = _Radius / _Image->ZSize();
    eval._LowerSigma     = _LowerIntensitySigma;
    eval._UpperSigma     = _UpperIntensitySigma;
    eval._LowerIntensity = PointData("Lower intensity", optional);
    eval._UpperIntensity = PointData("Upper intensity", optional);
    parallel_for(blocked_range<vtkIdType>(0, _NumberOfPoints), eval);
  }

  if (bg_fg_stats) {
    ComputeLocalIntensityStatistics eval;
    eval._Points               = _PointSet->Points();
    eval._Status               = _PointSet->Status();
    eval._Image                = _Image;
    eval._ForegroundMask       = (_ForegroundMask ? _ForegroundMask : &fg_mask);
    eval._BackgroundMask       = nullptr; // use foreground region of _Image
    eval._RadiusX              = _Radius / _Image->XSize();
    eval._RadiusY              = _Radius / _Image->YSize();
    eval._RadiusZ              = _Radius / _Image->ZSize();
    eval._BackgroundStatistics = PointData("Background statistics");
    eval._ForegroundStatistics = PointData("Foreground statistics");
    parallel_for(blocked_range<vtkIdType>(0, _NumberOfPoints), eval);
  }
}

// -----------------------------------------------------------------------------
void BalloonForce::Initialize()
{
  // Initialize base class
  SurfaceForce::Initialize();
  if (_NumberOfPoints == 0) return;

  // Check foreground mask attributes
  if (_ForegroundMask && !_ForegroundMask->HasSpatialAttributesOf(_Image)) {
    cerr << this->NameOfType() << "::Initialize: Mismatch of foreground mask and image attributes!" << endl;
    exit(1);
  }

  // Initial magnitude and direction of balloon force
  AddPointData("Signed magnitude")->FillComponent(0, _DeflateSurface ? -1.0 : 1.0);
  _DampingFactor      = max(.0, min(_DampingFactor,      1.0));
  _MagnitudeThreshold = max(.0, min(_MagnitudeThreshold, 1.0));

  // Radius of box window used for local intensity statistics
  // If zero, only the global intensity thresholds are considered
  if (_Radius < .0) {
    _Radius = 7.0 * max(max(_Image->XSize(),
                            _Image->YSize()),
                            _Image->ZSize());
  }

  // Add/remove point data arrays for local intensity thresholds
  // If sigma multiplier is negative, no local thresholds are used
  if (_Radius == .0 || _LowerIntensitySigma < .0) {
    RemovePointData("Lower intensity");
  } else {
    AddPointData("Lower intensity");
  }
  if (_Radius == .0 || _UpperIntensitySigma < .0) {
    RemovePointData("Upper intensity");
  } else {
    AddPointData("Upper intensity");
  }

  // Add/remove point data arrays for local inside/outside intensity statistics
  // If sigma multiplier is negative, no local inside/outside statistics are used
  if (_Radius == .0 || _ForegroundSigmaFactor < .0 || _BackgroundSigmaFactor < .0) {
    RemovePointData("Background statistics");
    RemovePointData("Foreground statistics");
  } else {
    AddPointData("Background statistics", 2);
    AddPointData("Foreground statistics", 2);
  }

  if (debug > 0) AddPointData("Intensity");
}

// -----------------------------------------------------------------------------
void BalloonForce::Update(bool gradient)
{
  // Note that _InitialUpdate is set to false by PointSetForce::Update
  const bool initial = _InitialUpdate;

  // Update base class
  SurfaceForce::Update(gradient);

  // Delayed initialization of local intensity thresholds
  // and update of local background/foreground statistics
  if (_Radius > .0) {
    const bool thresholds  = initial && (_LowerIntensitySigma >= 0. || _UpperIntensitySigma >= 0.);
    const bool bg_fg_stats = _BackgroundSigmaFactor > .0 && _ForegroundSigmaFactor > .0;
    ComputeLocalIntensityAttributes(thresholds, bg_fg_stats);
  }

  // Get (optional) point data (possibly interpolated during remeshing)
  const bool optional = true;
  vtkDataArray *magnitude       = PointData("Signed magnitude");
  vtkDataArray *intensity       = PointData("Intensity",             optional);
  vtkDataArray *lower_intensity = PointData("Lower intensity",       optional);
  vtkDataArray *upper_intensity = PointData("Upper intensity",       optional);
  vtkDataArray *bg_statistics   = PointData("Background statistics", optional);
  vtkDataArray *fg_statistics   = PointData("Foreground statistics", optional);

  // Initialize continuous intensity function
  BalloonForceUtils::ImageFunction image;
  image.Input(_Image);
  image.Initialize();

  // Update force magnitude and direction
  UpdateMagnitude update;
  update._Image                 = &image;
  update._Points                = Points();
  update._Status                = Status();
  update._DeflateSurface        = _DeflateSurface;
  update._Intensity             = intensity;
  update._LowerIntensity        = _LowerIntensity;
  update._UpperIntensity        = _UpperIntensity;
  update._LocalLowerIntensity   = lower_intensity;
  update._LocalUpperIntensity   = upper_intensity;
  update._BackgroundStatistics  = bg_statistics;
  update._BackgroundSigmaFactor = _BackgroundSigmaFactor;
  update._ForegroundStatistics  = fg_statistics;
  update._ForegroundSigmaFactor = _ForegroundSigmaFactor;
  update._Magnitude             = magnitude;
  update._MagnitudeDamping      = _DampingFactor;
  update._MagnitudeThreshold    = _MagnitudeThreshold;
  parallel_for(blocked_range<vtkIdType>(0, _NumberOfPoints), update);

  // Smooth magnitude such that adjacent nodes move coherently
  // (EXPERIMENTAL, see also PointSetForce::GradientAveraging)
  const int _MagnitudeSmoothing = 0;
  if (_MagnitudeSmoothing > 0) {
    vtkSmartPointer<vtkDataArray> smoothed_magnitude;
    smoothed_magnitude.TakeReference(magnitude->NewInstance());
    smoothed_magnitude->SetNumberOfComponents(magnitude->GetNumberOfComponents());
    smoothed_magnitude->SetNumberOfTuples(magnitude->GetNumberOfTuples());
    SmoothMagnitude smooth;
    smooth._Input        = magnitude;
    smooth._Output       = smoothed_magnitude;
    smooth._EdgeTable    = Edges();
    smooth._Normals      = Normals();
    smooth._MinMagnitude = _MagnitudeThreshold;
    blocked_range<vtkIdType> ptIds(0, magnitude->GetNumberOfTuples());
    for (int iter = 0; iter < _MagnitudeSmoothing; ++iter) {
      parallel_for(ptIds, smooth);
      swap(smooth._Input, smooth._Output);
    }
    if (smooth._Output != magnitude) {
      magnitude->CopyComponent(0, smoothed_magnitude, 0);
    }
  }
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void BalloonForce::EvaluateGradient(double *gradient, double step, double weight)
{
  if (_NumberOfPoints == 0) return;

  ComputeGradient eval;
  eval._Normals   = Normals();
  eval._Magnitude = PointData("Signed magnitude");
  eval._Gradient  = _Gradient;
  parallel_for(blocked_range<vtkIdType>(0, _NumberOfPoints), eval);

  SurfaceForce::EvaluateGradient(gradient, step, weight / _NumberOfPoints);
}


} // namespace mirtk
