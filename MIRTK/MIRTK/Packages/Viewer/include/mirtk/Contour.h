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

#ifndef _CONTOUR_H
#define _CONTOUR_H

#include <vector>

#include <mirtk/ViewerExport.h>
#include <mirtk/Image.h>


enum ContourMode { NewContour, NewPoint, CloseContour };

/// Class for storing the contour in viewer coordinates
class MIRTK_Viewer_EXPORT Contour
{

protected:

  /// Contour parts
  vector<PointSet> _pointSets;
  PointSet _allPoints;
  bool _updateAllPoints;

public:

  /// Constructor
  Contour();

  /// Destructor
  virtual ~Contour() {};

  virtual void Add(Point p);
  int IsEmpty();
  int Size();
  void Clear();
  int IsInside(double x, double y);
  virtual Point   &operator()(int);
  void AddNewSet(Point p);
  virtual void DeleteLastSet();
  void Print();

private:

  void AllPoints();

protected:

  void AddPointSet();

};

#endif
