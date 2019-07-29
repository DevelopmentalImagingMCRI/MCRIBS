/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

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

#include <mirtk/RView.h>
#include <mirtk/HistogramWindow.h>

HistogramWindow::HistogramWindow(RView  *viewer)
{
  _v = viewer;
}

void HistogramWindow::CalculateHistograms()
{
  int i;

  if (_v->GetTarget()->IsEmpty()) {
    cerr << "No target image loaded." << endl;
    return;
  }

  // Calculate global histogram
  i = -1;
  CalculateHistogram(i);

  // Calculate histogram for each structure
  for (i = 0; i < _v->GetSegmentTable()->Size(); i++) {
    if ( _v->GetSegmentTable()->IsValid(i) == true) CalculateHistogram(i);
  }
}

void HistogramWindow::CalculateHistogram(int label_id)
{
  int i, j, k, l;
  double min, max;
  double value;

  if (_v->GetTarget()->IsEmpty()) {
    cerr<< "No target image loaded." << endl;
    return;
  }

  _v->GetTarget()->GetMinMaxAsDouble(&min, &max);

  if (label_id < 0) {
    _globalHistogram.PutMin(min);
    _globalHistogram.PutMax(max);
    _globalHistogram.PutNumberOfBins(HISTOGRAM_BINS);

    for (l = 0; l < _v->GetTarget()->GetT(); l++){
      for (k = 0; k < _v->GetTarget()->GetZ(); k++){
        for (j = 0; j < _v->GetTarget()->GetY(); j++){
          for (i = 0; i < _v->GetTarget()->GetX(); i++){
               value = _v->GetTarget()->GetAsDouble(i, j, k, l);
          	if (value!=0)_globalHistogram.AddSample(value);
          }
        }
      }
    }
  } else {
    if (_v->GetSegmentTable()->IsValid(label_id) == true) {
      _localHistogram[label_id].PutMin(min);
      _localHistogram[label_id].PutMax(max);
      _localHistogram[label_id].PutNumberOfBins(HISTOGRAM_BINS);
      for (l = 0; l < _v->GetTarget()->GetT(); l++){
        for (k = 0; k < _v->GetTarget()->GetZ(); k++){
          for (j = 0; j < _v->GetTarget()->GetY(); j++){
            for (i = 0; i < _v->GetTarget()->GetX(); i++){
                value = _v->GetTarget()->GetAsDouble(i, j, k, l);
            	if ((value!=0)&&(_v->GetSegmentation()->Get(i, j, k, l) == label_id)) _globalHistogram.AddSample(value);
            }
          }
        }
      }
    }
  }
}

