/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright (c) Imperial College London
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

#ifndef _VOXELCONTOUR_H
#define _VOXELCONTOUR_H

#include <vector>

#include <mirtk/ViewerExport.h>


enum ContourMode
{
  FirstPoint,
  NewPoint,
  LastPoint
};


class MIRTK_Viewer_EXPORT VoxelContour
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
  void LineBresenham(int x0, int y0, int z0, int x1, int y1, int z);

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
