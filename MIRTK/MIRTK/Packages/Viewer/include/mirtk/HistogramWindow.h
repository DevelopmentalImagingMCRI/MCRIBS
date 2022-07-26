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

#ifndef _HISTOGRAMWINDOW_H
#define _HISTOGRAMWINDOW_H

#include <mirtk/ViewerExport.h>
#include <mirtk/Histogram1D.h>
#include <mirtk/RView.h>

#define HISTOGRAM_BINS 256

class MIRTK_Viewer_EXPORT HistogramWindow
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
