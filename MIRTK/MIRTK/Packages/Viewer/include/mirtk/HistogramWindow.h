/*=========================================================================

  Library   : Image Registration Toolkit ()
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _HISTOGRAMWINDOW_H
#define _HISTOGRAMWINDOW_H

#include <mirtk/Histogram1D.h>
#include <mirtk/RView.h>

#define HISTOGRAM_BINS 256

class HistogramWindow
{

  friend class Fl_HistogramWindow;

protected:

  /// Pointer to registration viewer
  RView *_v;

  /// Global histogram for entire image
  mirtk::Histogram1D<int> _globalHistogram;

  /// Global histogram for single segmentation
  mirtk::Histogram1D<int> _localHistogram[SHRT_MAX+1];

public:

  /// Constructor
  HistogramWindow(RView *);

  /// Destructor
  virtual ~HistogramWindow() {};

  /// Compute histograms for everything
  void CalculateHistograms();

protected:

  /// Compute histogram for single segmentation
  void CalculateHistogram(int);
};

#endif
