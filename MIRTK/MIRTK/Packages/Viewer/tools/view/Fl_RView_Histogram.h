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

#ifndef _FL_HISTOGRAMWINDOW_H
#define _FL_HISTOGRAMWINDOW_H

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/fl_draw.H>

#include <mirtk/RView.h>

class Fl_HistogramWindow : public Fl_Window
{

protected:

  /// Pointer to histogram window
  HistogramWindow _histogramWindow;

  /// Maximum in histogram
  int _maxHistogram;

public:

  /// Pointer to the registration viewer
  RView *_v;

  /// Constructor
  Fl_HistogramWindow(int, int, int, int, const char *, RView *);

  /// Destructor
  ~Fl_HistogramWindow();

  /// Default draw function
  void draw();

  /// Compute position
  void position(int, int, double& , double&);

  /// Recalculate histogram
  void recalculate();

};

#endif
