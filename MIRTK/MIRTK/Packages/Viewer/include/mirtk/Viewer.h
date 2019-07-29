/*=========================================================================

  Library   : Image Registration Toolkit ()
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _VIEWER_H
#define _VIEWER_H

#ifdef __APPLE__
#include <OpenGl/gl.h>
#include <OpenGl/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <mirtk/IOConfig.h>

class RView;
class VoxelContour;
class MultiLevelTransformation;
class FreeFormTransformation;

class Viewer
{

  /// Screen corrdinates
  int _screenX1, _screenY1, _screenX2, _screenY2;

  /// Viewport coordinates
  double _viewportX1, _viewportY1, _viewportX2, _viewportY2;

  /// Pointer to registration viewer
  RView *_rview;

  /// Viewer mode
  ViewerMode _viewerMode;

public:

  /// Constructor
  Viewer(RView *, ViewerMode);

  /// Destructor
  virtual ~Viewer();

  /// Draw image viewer
  virtual void DrawImage(Color *);

  /// Draw isolines in image viewer
  virtual void DrawIsolines(mirtk::GreyImage *, int);

  /// Draw segmentation contours in image viewer
  virtual void DrawSegmentationContour(mirtk::GreyImage *);

  /// Draw cursor in image viewer
  virtual void DrawCursor(CursorMode mode);

  /// Draw control points
  virtual void DrawPoints();

  /// Draw tag grid derived from target landmarks
  virtual void DrawTagGrid();

  /// Draw control points as deformation grid
  virtual void DrawGrid();

  /// Draw control points as deformation arrows
  virtual void DrawArrows();

  /// Draw landmarks
  void DrawLandmarks(mirtk::PointSet &, std::set<int> &, mirtk::GreyImage *, int = true, int = true);

  /// Draw landmark correspondences
  void DrawCorrespondences(mirtk::PointSet &, mirtk::PointSet &, mirtk::GreyImage *);

  /// Draw landmark correspondences
  void DrawCorrespondences(mirtk::PointSet &, mirtk::PointSet &, std::set<int> &, mirtk::GreyImage *);

  /// Draw ROI
  void DrawROI(mirtk::GreyImage *image, double, double, double, double,
               double, double);

#if MIRTK_IO_WITH_VTK
  /// Draw multiple objects
  void DrawObject(vtkPointSet **, mirtk::GreyImage *, int = false, int = false, mirtk::Transformation* = NULL);

  /// Draw object
  void DrawObject(vtkPointSet *, mirtk::GreyImage *, int = false, int = false, mirtk::Transformation* = NULL);
#endif

  /// Draw information about L/R, A/P, S/I on the viewer
  void DrawInfo(DisplayMode);

  /// Update Grid Pattern
  bool UpdateTagGrid(mirtk::GreyImage *, mirtk::Transformation *, mirtk::PointSet);

  /// Update using control points
  bool Update1(mirtk::GreyImage *, mirtk::MultiLevelTransformation *, mirtk::FreeFormTransformation *, double, double);

  /// Update using display resolution
  bool Update2(mirtk::GreyImage *, mirtk::MultiLevelTransformation *, mirtk::FreeFormTransformation *, double, double);

  /// Update
  bool Update(mirtk::GreyImage *, mirtk::Transformation *);

  /// Get width of viewer
  int GetWidth();

  /// Get height of viewer
  int GetHeight();

  /// Set screen (in normalized coordinates)
  void SetViewport(double, double, double, double);

  /// Get view port (in normalized coordinates)
  void GetViewport(double &, double &, double &, double &);

  /// Set screen coordinates (in pixel coordinates)
  void SetScreen(int, int);

  /// Get view port (in pixel coordinates)
  void GetScreen(int &, int &, int &, int &);

  /// Set viewer mode
  void SetViewerMode(ViewerMode);

  /// Get viewer mode
  ViewerMode GetViewerMode();

  /// Clipping of a drawable to the viewport
  void  Clip();

};

inline int Viewer::GetWidth()
{
  return _screenX2 - _screenX1 + 1;
}

inline int Viewer::GetHeight()
{
  return _screenY2 - _screenY1 + 1;
}

inline void Viewer::SetViewport(double x1, double y1, double x2, double y2)
{
  _viewportX1 = x1;
  _viewportY1 = y1;
  _viewportX2 = x2;
  _viewportY2 = y2;
}

inline void Viewer::GetViewport(double &x1, double &y1, double &x2, double &y2)
{
  x1 = _viewportX1;
  y1 = _viewportY1;
  x2 = _viewportX2;
  y2 = _viewportY2;
}

inline void Viewer::SetScreen(int x, int y)
{
  _screenX1 = round(x * _viewportX1);
  _screenY1 = round(y * _viewportY1);
  _screenX2 = round(x * _viewportX2);
  _screenY2 = round(y * _viewportY2);
}

inline void Viewer::GetScreen(int &x1, int &y1, int &x2, int &y2)
{
  x1 = _screenX1;
  y1 = _screenY1;
  x2 = _screenX2;
  y2 = _screenY2;
}

inline void Viewer::SetViewerMode(ViewerMode viewerMode)
{
  _viewerMode = viewerMode;
}

inline ViewerMode Viewer::GetViewerMode()
{
  return _viewerMode;
}

inline void Viewer::Clip()
{
  glViewport(_screenX1, _screenY1, this->GetWidth(), this->GetHeight());
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(_screenX1, _screenX2, _screenY1, _screenY2);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

#endif



