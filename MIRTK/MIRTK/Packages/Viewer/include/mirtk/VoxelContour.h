/*=========================================================================

  Library   : Image Registration Toolkit ()
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _VOXELCONTOUR_H
#define _VOXELCONTOUR_H

#include <vector>

typedef enum ContourMode { FirstPoint, NewPoint, LastPoint} VoxelContourMode;

class VoxelContour
{

  /// Pointer to rview
  RView *_rview;

  /// Width of paint brush
  int _width;

  /// Total number of points
  int _totalSize;

  /// Number of points currently drawn
  int _currentSize;

  /// Current selection
  int _current;

  /// First point drawn
  int _firstx, _firsty, _firstz;

  /// Last point drawn
  int _lastx, _lasty, _lastz;

  /// Adds pointset
  void AddPointSet();

  /// Line drawing
  void lineBresenham(int x0, int y0, int z0, int x1, int y1, int z);

  /// Fill area
  void Fill(int seedX, int seedY, int seedZ);

  /// Region growing
  void RegionGrowing2D(int seedX, int seedY, int seedZ, double lowT, double highT);

  /// Region growing
  void RegionGrowing3D(int seedX, int seedY, int seedZ, double lowT, double highT);

  /// Region growing criteria
  bool RegionGrowingCriteria(int i, int j, int k, double lowT, double highT);

public:

  /// Pointer to segmentation
  mirtk::GreyImage *_raster;

  /// Constructor
  VoxelContour();

  /// Initialise contour
  void Initialise(RView *, mirtk::GreyImage *);

  /// Operator for access
  mirtk::Point &operator()(int);

  /// Add a single point (in world coordinates)
  void AddPoint(mirtk::Point p, int width);

  /// Add a single point (in pixel coordinates)
  void AddPoint(int x, int y, int z);

  /// Add a segment
  void AddPointSet(mirtk::Point p, int width);

  /// Add a single point and connect to first point
  void Close(mirtk::Point p, int width);

  /// Undo: Remove last segment which has been added
  void Undo();

  /// Return number of points in contour
  int Size();

  /// Clear contour
  void Clear();

  /// Region growing
  void RegionGrowing(mirtk::Point, int thresholdMin, int thresholdMax, RegionGrowingMode);

  /// Fill area
  void FillArea(mirtk::Point);

};

inline int VoxelContour::Size()
{
  return _totalSize + _currentSize;
}

#endif
