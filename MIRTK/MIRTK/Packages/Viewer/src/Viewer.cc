/*=========================================================================

 Library   : Image Registration Toolkit (IRTK)
 Module    : $Id$
 Copyright : Imperial College, Department of Computing
 Visual Information Processing (VIP), 2008 onwards
 Date      : $Date$
 Version   : $Revision$
 Changes   : $Author$

 =========================================================================*/

#include <mirtk/RView.h>
#include <mirtk/Image.h>
#include <mirtk/Transformation.h>
#include <mirtk/Registration.h>

#ifdef __APPLE__
#include <OpenGl/gl.h>
#include <OpenGl/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif


// Define the maximum number of control points along each axis we can
// handle here.

#define MaxNumberOfCP 1000

// Define some global variables. These variables should NOT be used by
// any other classes. Perhaps they should be static members.

static int _NumberOfX;
static int _NumberOfY;
static double _BeforeX[MaxNumberOfCP][MaxNumberOfCP];
static double _BeforeY[MaxNumberOfCP][MaxNumberOfCP];
static double _BeforeZ[MaxNumberOfCP][MaxNumberOfCP];
static double _AfterX[MaxNumberOfCP][MaxNumberOfCP];
static double _AfterY[MaxNumberOfCP][MaxNumberOfCP];
static double _AfterZ[MaxNumberOfCP][MaxNumberOfCP];
static double _AfterGridX[MaxNumberOfCP][MaxNumberOfCP];
static double _AfterGridY[MaxNumberOfCP][MaxNumberOfCP];
static double _AfterGridZ[MaxNumberOfCP][MaxNumberOfCP];
static mirtk::Transformation::DOFStatus _CPStatus[MaxNumberOfCP][MaxNumberOfCP];

static int _NumberOfTagGridX;
static int _NumberOfTagGridY;
static double _AfterTagGridX [MaxNumberOfCP][MaxNumberOfCP];
static double _AfterTagGridY [MaxNumberOfCP][MaxNumberOfCP];
static double _AfterTagGridZ [MaxNumberOfCP][MaxNumberOfCP];
static double _BeforeTagGridX[MaxNumberOfCP][MaxNumberOfCP];
static double _BeforeTagGridY[MaxNumberOfCP][MaxNumberOfCP];
static double _BeforeTagGridZ[MaxNumberOfCP][MaxNumberOfCP];

// Define the default color scheme
#define COLOR_GRID                        glColor3f(1, 1, 0)
#define COLOR_ARROWS                      glColor3f(1, 1, 0)
#define COLOR_ISOLINES                    glColor3f(1, 1, 0)
#define COLOR_CONTOUR                     glColor4f(0, 1, 0, 0.5)
#define COLOR_CURSOR                      glColor3f(0, 1, 0)
#define COLOR_POINTS_ACTIVE               glColor3f(0, 1, 0)
#define COLOR_POINTS_PASSIVE              glColor3f(0, 0, 1)
#define COLOR_POINTS_UNKNOWN              glColor3f(1, 1, 0)
#define COLOR_CONTOUR_1                   glColor3f(1, 0, 0)
#define COLOR_CONTOUR_2                   glColor3f(0, 1, 0)
#define COLOR_CONTOUR_3                   glColor3f(0, 0, 1)
#define COLOR_CONTOUR_4                   glColor3f(1, 0, 1)
#define COLOR_CONTOUR_5                   glColor3f(0, 1, 1)
#define COLOR_TARGET_LANDMARKS            glColor3f(1, 0, 0)
#define COLOR_SOURCE_LANDMARKS            glColor3f(0, 0, 1)
#define COLOR_SELECTED_TARGET_LANDMARKS   glColor3f(1, 1, 0)
#define COLOR_SELECTED_SOURCE_LANDMARKS   glColor3f(0, 1, 1)

#ifdef HAVE_VTK

// vtk includes
#include <vtkPlane.h>
#include <vtkCutter.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkCell.h>
#include <vtkIdList.h>
#include <vtkGeometryFilter.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkUnstructuredGrid.h>

// object colour defines
#define COLOR_OBJECT           glColor3f(1, 0, 0)
#endif

GLubyte space[] =
{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

GLubyte letters[][13] =
{
		{ 0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xff, 0xc3, 0xc3, 0xc3, 0x66, 0x3c,
				0x18 },
		{ 0x00, 0x00, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7,
				0xfe },
		{ 0x00, 0x00, 0x7e, 0xe7, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xe7,
				0x7e },
		{ 0x00, 0x00, 0xfc, 0xce, 0xc7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc7, 0xce,
				0xfc },
		{ 0x00, 0x00, 0xff, 0xc0, 0xc0, 0xc0, 0xc0, 0xfc, 0xc0, 0xc0, 0xc0, 0xc0,
				0xff },
		{ 0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xfc, 0xc0, 0xc0, 0xc0,
				0xff },
		{ 0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xcf, 0xc0, 0xc0, 0xc0, 0xc0, 0xe7,
				0x7e },
		{ 0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xff, 0xc3, 0xc3, 0xc3, 0xc3,
				0xc3 },
		{ 0x00, 0x00, 0x7e, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18,
				0x7e },
		{ 0x00, 0x00, 0x7c, 0xee, 0xc6, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06,
				0x06 },
		{ 0x00, 0x00, 0xc3, 0xc6, 0xcc, 0xd8, 0xf0, 0xe0, 0xf0, 0xd8, 0xcc, 0xc6,
				0xc3 },
		{ 0x00, 0x00, 0xff, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0,
				0xc0 },
		{ 0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xdb, 0xff, 0xff, 0xe7,
				0xc3 },
		{ 0x00, 0x00, 0xc7, 0xc7, 0xcf, 0xcf, 0xdf, 0xdb, 0xfb, 0xf3, 0xf3, 0xe3,
				0xe3 },
		{ 0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xe7,
				0x7e },
		{ 0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7,
				0xfe },
		{ 0x00, 0x00, 0x3f, 0x6e, 0xdf, 0xdb, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0x66,
				0x3c },
		{ 0x00, 0x00, 0xc3, 0xc6, 0xcc, 0xd8, 0xf0, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7,
				0xfe },
		{ 0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x07, 0x7e, 0xe0, 0xc0, 0xc0, 0xe7,
				0x7e },
		{ 0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18,
				0xff },
		{ 0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3,
				0xc3 },
		{ 0x00, 0x00, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3,
				0xc3 },
		{ 0x00, 0x00, 0xc3, 0xe7, 0xff, 0xff, 0xdb, 0xdb, 0xc3, 0xc3, 0xc3, 0xc3,
				0xc3 },
		{ 0x00, 0x00, 0xc3, 0x66, 0x66, 0x3c, 0x3c, 0x18, 0x3c, 0x3c, 0x66, 0x66,
				0xc3 },
		{ 0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x3c, 0x3c, 0x66, 0x66,
				0xc3 },
		{ 0x00, 0x00, 0xff, 0xc0, 0xc0, 0x60, 0x30, 0x7e, 0x0c, 0x06, 0x03, 0x03,
				0xff } };

GLuint fontOffset;

// Little helper(s)
void status_glColor(int status)
{
	switch (status)
	{
	case mirtk::Status::Active:
		COLOR_POINTS_ACTIVE;
		break;
	case mirtk::Status::Passive:
		COLOR_POINTS_PASSIVE;
		break;
	default:
		COLOR_POINTS_UNKNOWN;
		break;
	}
}

Viewer::Viewer(RView *rview, ViewerMode viewerMode)
{
	_screenX1 = 0;
	_screenY1 = 0;
	_screenX2 = 0;
	_screenY2 = 0;
	_viewportX1 = 0;
	_viewportY1 = 0;
	_viewportX2 = 0;
	_viewportY2 = 0;

	// Pointer to registration viewer
	_rview = rview;

	// Mode of image viewer
	_viewerMode = viewerMode;
}

Viewer::~Viewer()
{
}

bool Viewer::UpdateTagGrid(mirtk::GreyImage *image, mirtk::Transformation *transformation, mirtk::PointSet landmark)
{
  if (landmark.Size() != 3) return false;

  mirtk::FreeFormTransformation   *affd = NULL;
  mirtk::MultiLevelTransformation *mffd = NULL;
  mirtk::Point p1, p2;
  double dx1, dy1, dx2, dy2, dz, x, y, tt, ts;
  int i, k, k1, k2, m, n;

  // Convert landmarks to image coordinates
  for (i = 0; i < landmark.Size(); i++) image->WorldToImage(landmark(i));

  // Check transformation
  mffd = dynamic_cast<mirtk::MultiLevelTransformation *>(transformation);

  if (mffd == NULL) {
    // Not a multi-level FFD, so let's try a single-level FFD
    affd = dynamic_cast<mirtk::FreeFormTransformation *>(transformation);
  } else {
    affd = mffd->GetLocalTransformation(mffd->NumberOfLevels() - 1);
  }

  // Find out time
  if (_rview->GetTarget() && 0 <= _rview->GetTargetFrame() && _rview->GetTargetFrame() < _rview->GetTarget()->GetT()) {
    tt = _rview->GetTarget()->ImageToTime(_rview->GetTargetFrame());
  } else {
    tt = affd->LatticeToTime(0);
  }
  if (_rview->GetSource() && 0 <= _rview->GetSourceFrame() && _rview->GetSourceFrame() < _rview->GetSource()->GetT()) {
    ts = _rview->GetSource()->ImageToTime(_rview->GetSourceFrame());
  } else {
    ts = affd->LatticeToTime(affd->GetT() - 1);
  }

  // Find out first corner of ROI
  landmark.BoundingBox(p1, p2);
  k1 = round(p1._z);
  k2 = round(p2._z);

  _NumberOfTagGridX = 13;
  _NumberOfTagGridY = 13;

  dx1 = (landmark(1)._x - landmark(0)._x) / 12.0;
  dy1 = (landmark(1)._y - landmark(0)._y) / 12.0;
  dx2 = (landmark(2)._x - landmark(1)._x) / 12.0;
  dy2 = (landmark(2)._y - landmark(1)._y) / 12.0;
  dz  = (k2 - k1) / 12;

  if (dz < 1) dz = 1;

  for (k = k1; k <= k2; k = k + dz) {
    for (m = 0; m < 13; m++) {
      for (n = 0; n < 13; n++) {
        x = landmark(0)._x + dx1 * n + dx2 * m;
        y = landmark(0)._y + dy1 * n + dy2 * m;

        _BeforeTagGridX[m][n] = x;
        _BeforeTagGridY[m][n] = y;
        _BeforeTagGridZ[m][n] = k;
        image->ImageToWorld(_BeforeTagGridX[m][n], _BeforeTagGridY[m][n], _BeforeTagGridZ[m][n]);

        _AfterTagGridX[m][n] = _BeforeTagGridX[m][n];
        _AfterTagGridY[m][n] = _BeforeTagGridY[m][n];
        _AfterTagGridZ[m][n] = _BeforeTagGridZ[m][n];
        if (_rview->GetSourceTransformApply()) {
          if (mffd != NULL) {
            mffd->LocalTransform(_AfterTagGridX[m][n], _AfterTagGridY[m][n], _AfterTagGridZ[m][n], ts, tt);
          } else if (affd != NULL) {
            affd->Transform     (_AfterTagGridX[m][n], _AfterTagGridY[m][n], _AfterTagGridZ[m][n], ts, tt);
          }
        }

        image->WorldToImage(_BeforeTagGridX[m][n], _BeforeTagGridY[m][n], _BeforeTagGridZ[m][n]);
        image->WorldToImage(_AfterTagGridX [m][n], _AfterTagGridY [m][n], _AfterTagGridZ [m][n]);
      }
    }
  }

  return true;
}

bool Viewer::Update1(mirtk::GreyImage *image, mirtk::MultiLevelTransformation *mffd, mirtk::FreeFormTransformation *affd, double ts, double tt)
{
	double x1, y1, z1, x2, y2, z2;
	int    i1, j1, k1, i2, j2, k2;
	int    index, i, j, k, m, n;

  // Check number of control points
  if ((affd->GetX() >= MaxNumberOfCP) || (affd->GetY() >= MaxNumberOfCP) || (affd->GetZ() >= MaxNumberOfCP)) {
    cerr << "Transformation has too many control points" << endl;
    exit(1);
  }

	// Find out first corner of ROI
	x1 = 0;
	y1 = 0;
	z1 = 0;
	image->ImageToWorld  (x1, y1, z1);
	affd ->WorldToLattice(x1, y1, z1);
	i1 = round(x1);
	j1 = round(y1);
	k1 = round(z1);

	if      (i1 < 0               ) i1 = 0;
	else if (i1 > affd->GetX() - 1) i1 = affd->GetX() - 1;
	if      (j1 < 0               ) j1 = 0;
	else if (j1 > affd->GetY() - 1) j1 = affd->GetY() - 1;
	if      (k1 < 0               ) k1 = 0;
	else if (k1 > affd->GetZ() - 1) k1 = affd->GetZ() - 1;

	// Find out second corner of ROI
	x2 = image->GetX() - 1;
	y2 = image->GetY() - 1;
	z2 = 0;
	image->ImageToWorld  (x2, y2, z2);
	affd ->WorldToLattice(x2, y2, z2);
	i2 = round(x2);
	j2 = round(y2);
	k2 = round(z2);

	if      (i2 < 0               ) i2 = 0;
	else if (i2 > affd->GetX() - 1) i2 = affd->GetX() - 1;
	if      (j2 < 0               ) j2 = 0;
	else if (j2 > affd->GetY() - 1) j2 = affd->GetY() - 1;
	if      (k2 < 0               ) k2 = 0;
	else if (k2 > affd->GetZ() - 1) k2 = affd->GetZ() - 1;

	// Swap if necessary
	if (i1 > i2) std::swap(i1, i2);
	if (j1 > j2) std::swap(j1, j2);
	if (k1 > k2) std::swap(k1, k2);

	switch (_viewerMode)
	{
	case Viewer_XY:
		_NumberOfX = i2 - i1 + 1;
		_NumberOfY = j2 - j1 + 1;
		break;
	case Viewer_XZ:
		_NumberOfX = i2 - i1 + 1;
		_NumberOfY = k2 - k1 + 1;
		break;
	case Viewer_YZ:
		_NumberOfX = j2 - j1 + 1;
		_NumberOfY = k2 - k1 + 1;
		break;
	default:
		break;
	}

	for (k = k1; k <= k2; k++) {
		for (j = j1; j <= j2; j++) {
			for (i = i1; i <= i2; i++) {
				// Calculate control points before and after deformation
				switch (_viewerMode)
				{
				case Viewer_XY:
					m = i - i1;
					n = j - j1;
					break;
				case Viewer_XZ:
					m = i - i1;
					n = k - k1;
					break;
				case Viewer_YZ:
					m = j - j1;
					n = k - k1;
					break;
				default:
					m = i;
					n = j;
					break;
				}

				_BeforeX[m][n] = i;
				_BeforeY[m][n] = j;
				_BeforeZ[m][n] = k;
				affd->LatticeToWorld(_BeforeX[m][n], _BeforeY[m][n], _BeforeZ[m][n]);

				_AfterX[m][n] = _BeforeX[m][n];
				_AfterY[m][n] = _BeforeY[m][n];
				_AfterZ[m][n] = _BeforeZ[m][n];
				if (mffd != NULL) {
          if (_rview->GetDisplayDeformationTotal()) {
            if (_rview->GetSourceTransformInvert()) {
              mffd->Inverse  (_AfterX[m][n], _AfterY[m][n], _AfterZ[m][n], ts, tt);
            } else {
              mffd->Transform(_AfterX[m][n], _AfterY[m][n], _AfterZ[m][n], ts, tt);
            }
          } else {
            if (_rview->GetSourceTransformInvert()) {
              mffd->LocalInverse  (_AfterX[m][n], _AfterY[m][n], _AfterZ[m][n], ts, tt);
            } else {
              mffd->LocalTransform(_AfterX[m][n], _AfterY[m][n], _AfterZ[m][n], ts, tt);
            }
          }
        } else {
          if (_rview->GetSourceTransformInvert()) {
            affd->Inverse  (_AfterX[m][n], _AfterY[m][n], _AfterZ[m][n], ts, tt);
          } else {
            affd->Transform(_AfterX[m][n], _AfterY[m][n], _AfterZ[m][n], ts, tt);
          }
        }

				image->WorldToImage(_BeforeX[m][n], _BeforeY[m][n], _BeforeZ[m][n]);
				image->WorldToImage(_AfterX[m][n], _AfterY[m][n], _AfterZ[m][n]);

				index = affd->LatticeToIndex(i, j, k);
        _CPStatus[m][n] = affd->IsActive(index) ? mirtk::Status::Active : mirtk::Status::Passive;

#ifdef HAVE__CPLABEL
				_CPLabel[m][n] = affd->GetLabel(index);
				if (_CPLabel[m][n] > _CPMaxLabel) _CPMaxLabel = _CPLabel[m][n];
				if (_CPLabel[m][n] < _CPMinLabel) _CPMinLabel = _CPLabel[m][n];
#endif

			}
		}
	}

	return true;
}

bool Viewer::Update2(mirtk::GreyImage *image, mirtk::MultiLevelTransformation *mffd, mirtk::FreeFormTransformation *affd, double ts, double tt)
{
	double dx, dy;
	int    i, j;

	dx = _rview->_DisplayDeformationGridResolution;
	dy = _rview->_DisplayDeformationGridResolution;
	_NumberOfX = round(static_cast<double>(this->GetWidth()  - 40) / dx);
	_NumberOfY = round(static_cast<double>(this->GetHeight() - 40) / dy);
	dx = (this->GetWidth()  - 40) / static_cast<double>(_NumberOfX);
	dy = (this->GetHeight() - 40) / static_cast<double>(_NumberOfY);

	for (j = 0; j < _NumberOfY; j++) {
		for (i = 0; i < _NumberOfX; i++) {
			_BeforeX[i][j] = i * dx + 20 + dx / 2.0;
			_BeforeY[i][j] = j * dy + 20 + dy / 2.0;
			_BeforeZ[i][j] = 0;
			_AfterX [i][j] = _BeforeX[i][j];
			_AfterY [i][j] = _BeforeY[i][j];
			_AfterZ [i][j] = 0;

			image->ImageToWorld(_AfterX[i][j], _AfterY[i][j], _AfterZ[i][j]);

			if (mffd != NULL) {
        if (_rview->GetDisplayDeformationTotal()) {
          if (_rview->GetSourceTransformInvert()) {
            mffd->Inverse  (_AfterX[i][j], _AfterY[i][j], _AfterZ[i][j], ts, tt);
          } else {
            mffd->Transform(_AfterX[i][j], _AfterY[i][j], _AfterZ[i][j], ts, tt);
          }
        } else {
          if (_rview->GetSourceTransformInvert()) {
            mffd->LocalInverse  (_AfterX[i][j], _AfterY[i][j], _AfterZ[i][j], ts, tt);
          } else {
            mffd->LocalTransform(_AfterX[i][j], _AfterY[i][j], _AfterZ[i][j], ts, tt);
          }
        }
			} else {
				if (_rview->GetSourceTransformInvert()) {
					affd->Inverse  (_AfterX[i][j], _AfterY[i][j], _AfterZ[i][j], ts, tt);
				} else {
					affd->Transform(_AfterX[i][j], _AfterY[i][j], _AfterZ[i][j], ts, tt);
				}
			}

			image->WorldToImage(_AfterX[i][j], _AfterY[i][j], _AfterZ[i][j]);

			//
			// FIXME was 
			// _CPStatus[i][j] = _Unknown;
			//
			_CPStatus[i][j] = mirtk::Status::Passive;
		}
	}

	return true;
}

bool Viewer::Update(mirtk::GreyImage *image, mirtk::Transformation *transformation)
{
  // Cast input transformation to single-/multi-level FFD
  mirtk::MultiLevelTransformation *mffd = dynamic_cast<mirtk::MultiLevelTransformation *>(transformation);
  mirtk::FreeFormTransformation   *affd = dynamic_cast<mirtk::FreeFormTransformation *>  (transformation);

  // Give up if no FFD given as input
  if (mffd == NULL && affd == NULL) {
    _NumberOfX = 0;
    _NumberOfY = 0;
    return false;
  }

  // If multi-level FFD, set affd to FFD of final level
  if (mffd != NULL) affd = mffd->GetLocalTransformation(mffd->NumberOfLevels() - 1);

  // Determine time parameters for transformation
  double tt, ts;
  if (_rview->GetTarget() && 0 <= _rview->GetTargetFrame() && _rview->GetTargetFrame() < _rview->GetTarget()->GetT()) {
    tt = _rview->GetTarget()->ImageToTime(_rview->GetTargetFrame());
  } else {
    tt = affd->LatticeToTime(0);
  }
  if (_rview->GetSource() && 0 <= _rview->GetSourceFrame() && _rview->GetSourceFrame() < _rview->GetSource()->GetT()) {
    ts = _rview->GetSource()->ImageToTime(_rview->GetSourceFrame());
  } else {
    ts = affd->LatticeToTime(affd->GetT() - 1);
  }

  // Compute (control) points before and after transformation
  // (application of mirtk::Transformation::Transform or mirtk::Transformation::Inverse)
  bool ok = false;
  if (_rview->_DisplayDeformationGridResolution == 0) {
    ok = this->Update1(image, mffd, affd, ts, tt);
  } else {
    ok = this->Update2(image, mffd, affd, ts, tt);
  }
  if (!ok) return false;

  // Deformation grid visualization
  if (_rview->GetDisplayDeformationGrid()) {
    // Copy points before transformation (space of target image)
    for (int j = 0; j < _NumberOfY; j++) {
      for (int i = 0; i < _NumberOfX; i++) {
        _AfterGridX[i][j] = _BeforeX[i][j];
        _AfterGridY[i][j] = _BeforeY[i][j];
        _AfterGridZ[i][j] = _BeforeZ[i][j];
      }
    }
    // If source image is being viewed, however, ...
    if (_rview->GetViewMode() == View_B) {
      // ... and if transformation is being applied, show inverse deformed grid
      if (_rview->GetSourceTransformApply()) {
        for (int j = 0; j < _NumberOfY; j++) {
          for (int i = 0; i < _NumberOfX; i++) {
            image->ImageToWorld(_AfterGridX[i][j], _AfterGridY[i][j], _AfterGridZ[i][j]);
            if (mffd != NULL) {
              if (_rview->GetDisplayDeformationTotal()) {
                if (_rview->GetSourceTransformInvert()) {
                  mffd->Transform(_AfterGridX[i][j], _AfterGridY[i][j], _AfterGridZ[i][j], tt, ts);
                } else {
                  mffd->Inverse  (_AfterGridX[i][j], _AfterGridY[i][j], _AfterGridZ[i][j], tt, ts);
                }
              } else {
                if (_rview->GetSourceTransformInvert()) {
                  mffd->LocalTransform(_AfterGridX[i][j], _AfterGridY[i][j], _AfterGridZ[i][j], tt, ts);
                } else {
                  mffd->LocalInverse  (_AfterGridX[i][j], _AfterGridY[i][j], _AfterGridZ[i][j], tt, ts);
                }
              }
            } else {
              if (_rview->GetSourceTransformInvert()) {
                affd->Transform(_AfterGridX[i][j], _AfterGridY[i][j], _AfterGridZ[i][j], tt, ts);
              } else {
                affd->Inverse  (_AfterGridX[i][j], _AfterGridY[i][j], _AfterGridZ[i][j], tt, ts);
              }
            }
            image->WorldToImage(_AfterGridX[i][j], _AfterGridY[i][j], _AfterGridZ[i][j]);
          }
        }
      }
    } else {
      // By default, just copy points after transformation (space of source image)
      for (int j = 0; j < _NumberOfY; j++) {
        for (int i = 0; i < _NumberOfX; i++) {
          _AfterGridX[i][j] = _AfterX[i][j];
          _AfterGridY[i][j] = _AfterY[i][j];
          _AfterGridZ[i][j] = _AfterZ[i][j];
        }
      }
    }
  }

  return true;
}

void Viewer::DrawCursor(CursorMode mode)
{
	int x, y;

	// Set color
	COLOR_CURSOR;

	// calculate width and height
	x = this->GetWidth();
	y = this->GetHeight();

	switch (mode)
	{
	case CrossHair:
		// Draw cross hair
		glBegin (GL_LINES);
		glVertex2f(_screenX1 + x / 2 - 10, _screenY1 + y / 2);
		glVertex2f(_screenX1 + x / 2 + 10, _screenY1 + y / 2);
		glVertex2f(_screenX1 + x / 2, _screenY1 + y / 2 - 10);
		glVertex2f(_screenX1 + x / 2, _screenY1 + y / 2 + 10);
		glEnd();
		break;
	case CursorX:
		// Draw cursor as broken 'X' (+)
		glBegin(GL_LINES);
		glVertex2f(_screenX1 + x / 2 - 10, _screenY1 + y / 2);
		glVertex2f(_screenX1 + x / 2 - 3, _screenY1 + y / 2);
		glVertex2f(_screenX1 + x / 2 + 3, _screenY1 + y / 2);
		glVertex2f(_screenX1 + x / 2 + 10, _screenY1 + y / 2);
		glVertex2f(_screenX1 + x / 2, _screenY1 + y / 2 - 10);
		glVertex2f(_screenX1 + x / 2, _screenY1 + y / 2 - 3);
		glVertex2f(_screenX1 + x / 2, _screenY1 + y / 2 + 3);
		glVertex2f(_screenX1 + x / 2, _screenY1 + y / 2 + 10);
		glEnd();
		break;
	case CursorV:
		// Draw cursor as 'V'
		glBegin(GL_LINES);
		glVertex2f(_screenX1 + x / 2 - 4, _screenY1 + y / 2 + 10);
		glVertex2f(_screenX1 + x / 2, _screenY1 + y / 2);
		glVertex2f(_screenX1 + x / 2, _screenY1 + y / 2);
		glVertex2f(_screenX1 + x / 2 + 4, _screenY1 + y / 2 + 10);
		glEnd();
		break;
	case CursorBar:
		// Draw cursor as bar with scales
		break;
	}
}

void Viewer::DrawIsolines(mirtk::GreyImage *image, int value)
{
	int i, j;

	// Set color
	COLOR_ISOLINES;
	glLineWidth(_rview->GetLineThickness());
	glBegin (GL_LINES);

	for (j = 0; j < this->GetHeight() - 1; j++) {
		for (i = 0; i < this->GetWidth() - 1; i++) {
			if (((image->Get(i, j, 0) <= value) && (image->Get(i + 1, j, 0) > value))
					|| ((image->Get(i, j, 0) > value)
							&& (image->Get(i + 1, j, 0) <= value))) {
				glVertex2f(_screenX1 + i + 0.5, _screenY1 + j - 0.5);
				glVertex2f(_screenX1 + i + 0.5, _screenY1 + j + 0.5);

			}
			if (((image->Get(i, j, 0) <= value) && (image->Get(i, j + 1, 0) > value))
					|| ((image->Get(i, j, 0) > value)
							&& (image->Get(i, j + 1, 0) <= value))) {
				glVertex2f(_screenX1 + i + 0.5, _screenY1 + j + 0.5);
				glVertex2f(_screenX1 + i - 0.5, _screenY1 + j + 0.5);
			}
		}
	}
	glEnd();
	glLineWidth(1);
}

void Viewer::DrawSegmentationContour(mirtk::GreyImage *image)
{
	int i, j;
	unsigned char r, g, b;
	glLineWidth(_rview->GetLineThickness());

	for (j = 1; j < this->GetHeight() - 1; j++) {
		for (i = 1; i < this->GetWidth() - 1; i++) {
			if ((image->Get(i, j, 0) > 0)
					&& (_rview->_segmentTable->_entry[image->Get(i, j, 0)]._visible
							== true)) {
				r = _rview->_segmentTable->_entry[image->Get(i, j, 0)]._color.r;
				g = _rview->_segmentTable->_entry[image->Get(i, j, 0)]._color.g;
				b = _rview->_segmentTable->_entry[image->Get(i, j, 0)]._color.b;
				if (image->Get(i, j, 0) != image->Get(i + 1, j, 0)) {
					glColor3ub(r, g, b);
					glBegin (GL_LINES);
					glVertex2f(_screenX1 + i, _screenY1 + j - 0.5);
					glVertex2f(_screenX1 + i, _screenY1 + j + 0.5);
					glEnd();
				}
				if (image->Get(i, j, 0) != image->Get(i - 1, j, 0)) {
					glColor3ub(r, g, b);
					glBegin (GL_LINES);
					glVertex2f(_screenX1 + i, _screenY1 + j - 0.5);
					glVertex2f(_screenX1 + i, _screenY1 + j + 0.5);
					glEnd();
				}
				if (image->Get(i, j, 0) != image->Get(i, j + 1, 0)) {
					glColor3ub(r, g, b);
					glBegin (GL_LINES);
					glVertex2f(_screenX1 + i + 0.5, _screenY1 + j);
					glVertex2f(_screenX1 + i - 0.5, _screenY1 + j);
					glEnd();
				}
				if (image->Get(i, j, 0) != image->Get(i, j - 1, 0)) {
					glColor3ub(r, g, b);
					glBegin (GL_LINES);
					glVertex2f(_screenX1 + i + 0.5, _screenY1 + j);
					glVertex2f(_screenX1 + i - 0.5, _screenY1 + j);
					glEnd();
				}
			}
		}
	}
	glLineWidth(1);
}

void Viewer::DrawTagGrid()
{
  int i, j;

  // Set color
  COLOR_GRID;
  glLineWidth(_rview->GetLineThickness());

  glBegin (GL_LINES);
  for (j = 0; j < _NumberOfTagGridY; j++) {
    for (i = 0; i < _NumberOfTagGridX - 1; i++) {
      glVertex2f(_screenX1 + _AfterTagGridX[i    ][j], _screenY1 + _AfterTagGridY[i    ][j]);
      glVertex2f(_screenX1 + _AfterTagGridX[i + 1][j], _screenY1 + _AfterTagGridY[i + 1][j]);
    }
  }
  for (j = 0; j < _NumberOfTagGridY - 1; j++) {
    for (i = 0; i < _NumberOfTagGridX; i++) {
      glVertex2f(_screenX1 + _AfterTagGridX[i][j    ], _screenY1 + _AfterTagGridY[i][j    ]);
      glVertex2f(_screenX1 + _AfterTagGridX[i][j + 1], _screenY1 + _AfterTagGridY[i][j + 1]);
    }
  }
  glEnd();
  glLineWidth(1);
}

void Viewer::DrawGrid()
{
	int i, j;

	// Set color
	COLOR_GRID;

	glLineWidth(_rview->GetLineThickness());

	glBegin (GL_LINES);
	for (j = 0; j < _NumberOfY; j++) {
		for (i = 0; i < _NumberOfX - 1; i++) {
			glVertex2f(_screenX1 + _AfterGridX[i    ][j], _screenY1 + _AfterGridY[i    ][j]);
			glVertex2f(_screenX1 + _AfterGridX[i + 1][j], _screenY1 + _AfterGridY[i + 1][j]);
		}
	}
	for (j = 0; j < _NumberOfY - 1; j++) {
		for (i = 0; i < _NumberOfX; i++) {
			glVertex2f(_screenX1 + _AfterGridX[i][j    ], _screenY1 + _AfterGridY[i][j    ]);
			glVertex2f(_screenX1 + _AfterGridX[i][j + 1], _screenY1 + _AfterGridY[i][j + 1]);
		}
	}
	glEnd();

	glLineWidth(1);
}

void Viewer::DrawArrows()
{
	int i, j;

	// Set color
	COLOR_ARROWS;

	for (j = 0; j < _NumberOfY; j++) {
		for (i = 0; i < _NumberOfX; i++) {
			glBegin (GL_LINES);
			glVertex2f(_screenX1 + _BeforeX[i][j], _screenY1 + _BeforeY[i][j]);
			glVertex2f(_screenX1 + _AfterX[i][j], _screenY1 + _AfterY[i][j]);
			glEnd();
			float dx = _AfterX[i][j] - _BeforeX[i][j];
			float dy = _AfterY[i][j] - _BeforeY[i][j];
			float fat_factor = 2.0;
			float line_len = sqrt(dx * dx + dy * dy);
			float archor_width = line_len / 6.0;
			if (line_len > 0.01) {
				float factor = fat_factor * archor_width / line_len;
				float add1dx = (dx * factor);
				float add1dy = (dy * factor);
				float add2dx = (dy * factor / (2 * fat_factor));
				float add2dy = (-dx * factor / (2 * fat_factor));
				mirtk::Point point[3];
				point[0]._x = _screenX1 + _AfterX[i][j];
				point[0]._y = _screenY1 + _AfterY[i][j];
				point[1]._x = point[0]._x - add1dx + add2dx;
				point[1]._y = point[0]._y - add1dy + add2dy;
				point[2]._x = point[0]._x - add1dx - add2dx;
				point[2]._y = point[0]._y - add1dy - add2dy;
				glBegin (GL_POLYGON);
				glVertex2f(point[0]._x, point[0]._y);
				glVertex2f(point[1]._x, point[1]._y);
				glVertex2f(point[2]._x, point[2]._y);
				glEnd();
			}
		}
	}
}

void Viewer::DrawPoints()
{
	int i, j;

	// Adjust pointsize
	glPointSize((GLfloat) 3);

	// Draw active and passive points
	glBegin (GL_POINTS);
	for (j = 0; j < _NumberOfY; j++) {
		for (i = 0; i < _NumberOfX; i++) {
			// Set color
			status_glColor(_CPStatus[i][j]);
			// Draw point
			glVertex2f(_screenX1 + _BeforeX[i][j], _screenY1 + _BeforeY[i][j]);
		}
	}
	glEnd();
}

void Viewer::DrawLandmarks(mirtk::PointSet &landmarks, std::set<int> &ids, mirtk::GreyImage *image, int bTarget, int bAll)
{
  glLineWidth(1.0);
  // Draw unselected landmarks first
  if (bAll) {
    for (int i = 0; i < landmarks.Size(); ++i) {
      // Skip selected landmarks
      if (ids.find(i) != ids.end()) continue;
      // Adjust colour
      if (bTarget) COLOR_TARGET_LANDMARKS;
      else         COLOR_SOURCE_LANDMARKS;
      // Get landmark point
      mirtk::Point p = landmarks(i);
      // Transform point
      if (!bTarget && _rview->GetSourceTransformApply()) {
        const bool inv = !_rview->GetSourceTransformInvert();
        if (inv) _rview->_sourceTransform->Inverse  (p);
        else     _rview->_sourceTransform->Transform(p);
      }
      // Draw point
      if (image->IsInFOV(p._x, p._y, p._z)) {
        image->WorldToImage(p);
        glBegin(GL_LINES);
        glVertex2f(_screenX1 + p._x - 8, _screenY1 + p._y);
        glVertex2f(_screenX1 + p._x + 8, _screenY1 + p._y);
        glVertex2f(_screenX1 + p._x, _screenY1 + p._y - 8);
        glVertex2f(_screenX1 + p._x, _screenY1 + p._y + 8);
        glEnd();
      }
    }
  }
  // Draw selected landmarks on top
  for (std::set<int>::const_iterator i = ids.begin(); i != ids.end(); ++i) {
    if (*i < 0 || *i >= landmarks.Size()) continue;
    // Adjust colour
    if (bTarget) COLOR_SELECTED_TARGET_LANDMARKS;
    else         COLOR_SELECTED_SOURCE_LANDMARKS;
    // Get landmark point
    mirtk::Point p = landmarks(*i);
    // Transform point
    if (!bTarget && _rview->GetSourceTransformApply()) {
      const bool inv = !_rview->GetSourceTransformInvert();
      if (inv) _rview->_sourceTransform->Inverse  (p);
      else     _rview->_sourceTransform->Transform(p);
    }
    // Draw point
    if (image->IsInFOV(p._x, p._y, p._z)) {
      image->WorldToImage(p);
      glBegin(GL_LINES);
      glVertex2f(_screenX1 + p._x - 8, _screenY1 + p._y);
      glVertex2f(_screenX1 + p._x + 8, _screenY1 + p._y);
      glVertex2f(_screenX1 + p._x, _screenY1 + p._y - 8);
      glVertex2f(_screenX1 + p._x, _screenY1 + p._y + 8);
      glEnd();
    }
  }
}

void Viewer::DrawCorrespondences(mirtk::PointSet &target, mirtk::PointSet &source, mirtk::GreyImage *image)
{
  // Adjust colour
  glLineWidth(1.0);
  glColor3f(0, 1, 0);

  // Draw lines connecting corresponding landmarks
  mirtk::Point p1, p2;
  for (int i = 0; i < target.Size(); ++i) {
    if (i >= source.Size()) continue;
    p1 = target(i);
    p2 = source(i);
    _rview->_sourceTransform->Inverse(p1._x, p1._y, p1._z);
    if (image->IsInFOV(p1._x, p1._y, p1._z) &&
        image->IsInFOV(p2._x, p2._y, p2._z)) {
      image->WorldToImage(p1);
      image->WorldToImage(p2);
      glBegin (GL_LINES);
      glVertex2f(_screenX1 + p1._x, _screenY1 + p1._y);
      glVertex2f(_screenX1 + p2._x, _screenY1 + p2._y);
      glEnd();
    }
  }
}

void Viewer::DrawCorrespondences(mirtk::PointSet &target, mirtk::PointSet &source, std::set<int> &ids, mirtk::GreyImage *image)
{
  // Adjust colour
  glLineWidth(1.0);
  glColor3f(0, 1, 0);

  // Draw lines connecting corresponding landmarks
  mirtk::Point p1, p2;
  for (std::set<int>::const_iterator id = ids.begin(); id != ids.end(); ++id) {
    if (*id < 0 || *id >= target.Size() || *id >= source.Size()) continue;
    p1 = target(*id);
    p2 = source(*id);
    _rview->_sourceTransform->Inverse(p1._x, p1._y, p1._z);
    if (image->IsInFOV(p1._x, p1._y, p1._z) &&
        image->IsInFOV(p2._x, p2._y, p2._z)) {
      image->WorldToImage(p1);
      image->WorldToImage(p2);
      glBegin (GL_LINES);
      glVertex2f(_screenX1 + p1._x, _screenY1 + p1._y);
      glVertex2f(_screenX1 + p2._x, _screenY1 + p2._y);
      glEnd();
    }
  }
}

void Viewer::DrawImage(Color *drawable)
{
	// Set raster position
	glRasterPos2f(_screenX1, _screenY1);

	// Draw pixelmap
	glDrawPixels(this->GetWidth(), this->GetHeight(), GL_RGB, GL_UNSIGNED_BYTE,
			drawable);
}

void Viewer::DrawROI(mirtk::GreyImage *image, double x1, double y1, double z1,
		double x2, double y2, double z2)
{
	image->WorldToImage(x1, y1, z1);
	image->WorldToImage(x2, y2, z2);
	glColor4f(1, 0, 0, 1);
	glBegin (GL_POLYGON);
	glVertex2f(_screenX1 + x1 - 2, _screenY1 + y1 - 2);
	glVertex2f(_screenX1 + x1 + 2, _screenY1 + y1 - 2);
	glVertex2f(_screenX1 + x1 + 2, _screenY1 + y1 + 2);
	glVertex2f(_screenX1 + x1 - 2, _screenY1 + y1 + 2);
	glEnd();
	glColor4f(0, 1, 0, 1);
	glBegin(GL_POLYGON);
	glVertex2f(_screenX1 + x2 - 2, _screenY1 + y2 - 2);
	glVertex2f(_screenX1 + x2 + 2, _screenY1 + y2 - 2);
	glVertex2f(_screenX1 + x2 + 2, _screenY1 + y2 + 2);
	glVertex2f(_screenX1 + x2 - 2, _screenY1 + y2 + 2);
	glEnd();
	glColor4f(1, 1, 0, 0.5);
	glEnable (GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glBegin(GL_POLYGON);
	glVertex2f(_screenX1 + x1, _screenY1 + y1);
	glVertex2f(_screenX1 + x1, _screenY1 + y2);
	glVertex2f(_screenX1 + x2, _screenY1 + y2);
	glVertex2f(_screenX1 + x2, _screenY1 + y1);
	glEnd();
	glDisable(GL_BLEND);
}

#ifdef HAVE_VTK

void Viewer::DrawObject(vtkPointSet **object, mirtk::GreyImage *image,
		int _DisplayObjectWarp,
		int _DisplayObjectGrid,
		mirtk::Transformation *transformation)
{
	int i;

	for (i = 0; i < _rview->_NoOfObjects; i++) {
		switch (i) {
			case 0:
			COLOR_CONTOUR_1;
			break;
			case 1:
			COLOR_CONTOUR_2;
			break;
			case 2:
			COLOR_CONTOUR_3;
			break;
			case 3:
			COLOR_CONTOUR_4;
			break;
			default:
			COLOR_CONTOUR_5;
			break;
		}

        glLineWidth(_rview->GetLineThickness());
		this->DrawObject(object[i], image, _DisplayObjectWarp, _DisplayObjectGrid, transformation);
	}
}

void Viewer::DrawObject(vtkPointSet *points, mirtk::GreyImage *image, int warp, int grid, mirtk::Transformation *transformation)
{
	int i, j;
	double p1[3], p2[3], p3[3], v1[3], v2[3], point[3], normal[3];
	mirtk::Point p;
	static vtkPlane *plane = vtkPlane::New();
	static vtkCutter *cutter = vtkCutter::New();

	if (points != NULL) {
		p1[0] = 0;
		p1[1] = 0;
		p1[2] = 0;
		p2[0] = image->GetX();
		p2[1] = 0;
		p2[2] = 0;
		p3[0] = 0;
		p3[1] = image->GetY();
		p3[2] = 0;
		image->ImageToWorld(p1[0], p1[1], p1[2]);
		image->ImageToWorld(p2[0], p2[1], p2[2]);
		image->ImageToWorld(p3[0], p3[1], p3[2]);
		v1[0] = p2[0] - p1[0];
		v1[1] = p2[1] - p1[1];
		v1[2] = p2[2] - p1[2];
		v2[0] = p3[0] - p1[0];
		v2[1] = p3[1] - p1[1];
		v2[2] = p3[2] - p1[2];
		normal[0] = v1[1] * v2[2] - v1[2] * v2[1];
		normal[1] = v1[2] * v2[0] - v1[0] * v2[2];
		normal[2] = v1[0] * v2[1] - v1[1] * v2[0];

		// Set up plane and cutter
		plane->SetOrigin(p1);
		plane->SetNormal(normal);
		cutter->SetCutFunction(plane);
#if VTK_MAJOR_VERSION >= 6
    cutter->SetInputData(points);
#else
		cutter->SetInput(points);
#endif

		// Reslice object
		cutter->Modified();
		cutter->Update();

		// Loop over cells
		for (i = 0; i < cutter->GetOutput()->GetNumberOfCells(); i++) {
			// Get pointIds from cell
			vtkIdList *ids = cutter->GetOutput()->GetCell(i)->GetPointIds();
			mirtk::PointSet pset;
			for (j = 0; j < ids->GetNumberOfIds(); j++) {

				// Get point from cell
				cutter->GetOutput()->GetPoints()->GetPoint(ids->GetId(j), point);

				// Transform the point if warp is on
				if(warp && transformation != NULL) {
					transformation->Transform(point[0], point[1], point[2]);
				}

				// Get point in image coordinates
				image->WorldToImage(point[0], point[1], point[2]);
				switch (_viewerMode) {
					case Viewer_XY:
					p = mirtk::Point(point[0], point[1], 0);
					break;
					case Viewer_YZ:
					p = mirtk::Point(point[0], point[1], 0);
					break;
					case Viewer_XZ:
					p = mirtk::Point(point[0], point[1], 0);
					break;
					default:
					break;
				}
				pset.Add(p);
			}

			// Now draw
			for (j = 0; j < pset.Size(); j++) {
				glBegin(GL_LINES);
				glVertex2f(_screenX1+pset(j)._x,
						_screenY1+pset(j)._y);
				glVertex2f(_screenX1+pset((j+1)%pset.Size())._x,
						_screenY1+pset((j+1)%pset.Size())._y);
				glEnd();
			}
		}
	}
}

#endif

void Viewer::DrawInfo(DisplayMode m)
{
	int x, y;
	static int first = true;

	if (first == true) {
		GLuint i, j;
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

		fontOffset = glGenLists(128);
		for (i = 0, j = 'A'; i < 26; i++, j++) {
			glNewList(fontOffset + j, GL_COMPILE);
			glBitmap(8, 13, 0.0, 2.0, 10.0, 0.0, letters[i]);
			glEndList();
		}
		glNewList(fontOffset + ' ', GL_COMPILE);
		glBitmap(8, 13, 0.0, 2.0, 10.0, 0.0, space);
		glEndList();
		first = false;
	}

	if (m == Native)
		return;

	// calculate width and height
	x = this->GetWidth();
	y = this->GetHeight();

	switch (_viewerMode)
	{
	case Viewer_XY:

		if (m == Neurological) {

			// Draw axis labels
			glRasterPos2i(_screenX1 + 5, _screenY1 + y / 2 - 5);
			glPushAttrib (GL_LIST_BIT);
			glListBase(fontOffset);
			glCallLists(strlen("L"), GL_UNSIGNED_BYTE, (GLubyte *) "R");
			glPopAttrib();

			glRasterPos2i(_screenX1 + x - 15, _screenY1 + y / 2 - 5);
			glPushAttrib(GL_LIST_BIT);
			glListBase(fontOffset);
			glCallLists(strlen("R"), GL_UNSIGNED_BYTE, (GLubyte *) "L");
			glPopAttrib();

		}
		else {

			// Draw axis labels
			glRasterPos2i(_screenX1 + 5, _screenY1 + y / 2 - 5);
			glPushAttrib (GL_LIST_BIT);
			glListBase(fontOffset);
			glCallLists(strlen("L"), GL_UNSIGNED_BYTE, (GLubyte *) "L");
			glPopAttrib();

			glRasterPos2i(_screenX1 + x - 15, _screenY1 + y / 2 - 5);
			glPushAttrib(GL_LIST_BIT);
			glListBase(fontOffset);
			glCallLists(strlen("R"), GL_UNSIGNED_BYTE, (GLubyte *) "R");
			glPopAttrib();

		}

		glRasterPos2i(_screenX1 + x / 2 - 5, _screenY1 + 5);
		glPushAttrib (GL_LIST_BIT);
		glListBase(fontOffset);
		glCallLists(strlen("P"), GL_UNSIGNED_BYTE, (GLubyte *) "P");
		glPopAttrib();

		glRasterPos2i(_screenX1 + x / 2 - 5, _screenY1 + y - 15);
		glPushAttrib(GL_LIST_BIT);
		glListBase(fontOffset);
		glCallLists(strlen("A"), GL_UNSIGNED_BYTE, (GLubyte *) "A");
		glPopAttrib();
		break;

	case Viewer_XZ:

		if (m == Neurological) {

			// Draw axis labels
			glRasterPos2i(_screenX1 + 5, _screenY1 + y / 2 - 5);
			glPushAttrib(GL_LIST_BIT);
			glListBase(fontOffset);
			glCallLists(strlen("R"), GL_UNSIGNED_BYTE, (GLubyte *) "R");
			glPopAttrib();

			glRasterPos2i(_screenX1 + x - 15, _screenY1 + y / 2 - 5);
			glPushAttrib(GL_LIST_BIT);
			glListBase(fontOffset);
			glCallLists(strlen("L"), GL_UNSIGNED_BYTE, (GLubyte *) "L");
			glPopAttrib();

		}
		else {

			// Draw axis labels
			glRasterPos2i(_screenX1 + 5, _screenY1 + y / 2 - 5);
			glPushAttrib(GL_LIST_BIT);
			glListBase(fontOffset);
			glCallLists(strlen("R"), GL_UNSIGNED_BYTE, (GLubyte *) "L");
			glPopAttrib();

			glRasterPos2i(_screenX1 + x - 15, _screenY1 + y / 2 - 5);
			glPushAttrib(GL_LIST_BIT);
			glListBase(fontOffset);
			glCallLists(strlen("L"), GL_UNSIGNED_BYTE, (GLubyte *) "R");
			glPopAttrib();

		}

		glRasterPos2i(_screenX1 + x / 2 - 5, _screenY1 + 5);
		glPushAttrib(GL_LIST_BIT);
		glListBase(fontOffset);
		glCallLists(strlen("I"), GL_UNSIGNED_BYTE, (GLubyte *) "I");
		glPopAttrib();

		glRasterPos2i(_screenX1 + x / 2 - 5, _screenY1 + y - 15);
		glPushAttrib(GL_LIST_BIT);
		glListBase(fontOffset);
		glCallLists(strlen("S"), GL_UNSIGNED_BYTE, (GLubyte *) "S");
		glPopAttrib();
		break;

	case Viewer_YZ:

		// Draw axis labels
		glRasterPos2i(_screenX1 + 5, _screenY1 + y / 2 - 5);
		glPushAttrib(GL_LIST_BIT);
		glListBase(fontOffset);
		glCallLists(strlen("P"), GL_UNSIGNED_BYTE, (GLubyte *) "P");
		glPopAttrib();

		glRasterPos2i(_screenX1 + x - 15, _screenY1 + y / 2 - 5);
		glPushAttrib(GL_LIST_BIT);
		glListBase(fontOffset);
		glCallLists(strlen("A"), GL_UNSIGNED_BYTE, (GLubyte *) "A");
		glPopAttrib();

		glRasterPos2i(_screenX1 + x / 2 - 5, _screenY1 + 5);
		glPushAttrib(GL_LIST_BIT);
		glListBase(fontOffset);
		glCallLists(strlen("I"), GL_UNSIGNED_BYTE, (GLubyte *) "I");
		glPopAttrib();

		glRasterPos2i(_screenX1 + x / 2 - 5, _screenY1 + y - 15);
		glPushAttrib(GL_LIST_BIT);
		glListBase(fontOffset);
		glCallLists(strlen("S"), GL_UNSIGNED_BYTE, (GLubyte *) "S");
		glPopAttrib();
		break;

	default:
		break;
	}
}

