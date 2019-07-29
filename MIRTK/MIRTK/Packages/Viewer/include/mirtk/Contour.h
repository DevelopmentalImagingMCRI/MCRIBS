/*=========================================================================

  Library   : Image Registration Toolkit ()
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _CONTOUR_H
#define _CONTOUR_H

#include <vector>

#include <mirtk/Image.h>


enum ContourMode { NewContour, NewPoint, CloseContour };

/// Class for storing the contour in viewer coordinates
class Contour
{

protected:

  /// Contour parts
  vector<PointSet> _pointSets;
  PointSet _allPoints;
  int _updateAllPoints;

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
