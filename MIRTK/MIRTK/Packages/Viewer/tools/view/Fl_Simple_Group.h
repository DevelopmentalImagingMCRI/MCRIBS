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

#ifndef _FL_SIMPLE_GROUP_H
#define _FL_SIMPLE_GROUP_H

/* fltk includes */
#include <FL/Fl.H>
#include <FL/fl_draw.H>
#include <FL/Fl_Group.H>

//! This class provides a simple aesthetic alternative to Fl_Group
class Fl_Simple_Group : public Fl_Group
{

public:

  //! Normal FLTK constructor
  Fl_Simple_Group(int, int, int, int, const char *);

  //! Override of Fl_Group::draw()
  void draw();

};

#endif
