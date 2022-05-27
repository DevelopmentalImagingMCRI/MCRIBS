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

#ifndef FL_SIMPLE_BROWSER_H
#define FL_SIMPLE_BROWSER_H

#include <FL/Fl.H>
#include <FL/Fl_Browser.H>
#include <FL/Fl_Hold_Browser.H>

/// This class provides a simple alternative to Fl_Browser
class Fl_Simple_Browser : public Fl_Hold_Browser
{

public:

  /// Normal FLTK constructor
  Fl_Simple_Browser(int x, int y, int w, int h, const char *label = 0);

  /// Override of Fl_Hold_Browser::handle()
  int handle(int);

};

inline Fl_Simple_Browser::Fl_Simple_Browser(int x, int y, int w, int h, const char *label) : Fl_Hold_Browser(x, y, w, h, label)
{
}

#endif
