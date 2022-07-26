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

#ifndef _FL_RVIEW_H
#define _FL_RVIEW_H

#include <mirtk/OpenGl.h>

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Gl_Window.H>

#include <mirtk/RView.h>


class Fl_RView : public Fl_Gl_Window
{

public:

  /// Pointer to the registration viewer
  RView *v;

  /// Constructor
  Fl_RView(int, int, int, int, const char *);

  /// Default draw function
  void draw();

  /// Default function to handle events
  int handle(int);

  #if FL_API_VERSION < 10304
  int pixel_w()
  {
    return w();
  }

  int pixel_h()
  {
    return h();
  }

  float pixels_per_unit()
  {
    return 1;
  }
  #endif

  int pixel_event_x()
  {
    return int(pixels_per_unit() * Fl::event_x() + 0.5);
  }

  int pixel_event_y()
  {
    return int(pixels_per_unit() * Fl::event_y() + 0.5);
  }

  int pixel_event_dx()
  {
    return int(pixels_per_unit() * Fl::event_dx() + 0.5);
  }

  int pixel_event_dy()
  {
    return int(pixels_per_unit() * Fl::event_dy() + 0.5);
  }

};

#endif
