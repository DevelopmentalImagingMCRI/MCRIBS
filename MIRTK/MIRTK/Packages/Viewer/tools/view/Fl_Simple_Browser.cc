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

#include "Fl_Simple_Browser.h"

#include <string.h>
#include <stdio.h>

int Fl_Simple_Browser::handle(int event)
{
  switch (event) {
  case FL_KEYBOARD:
  case FL_SHORTCUT:
    if (Fl::event_key() == FL_BackSpace) {
      do_callback();
      return 1;
    }
    return this->Fl_Browser_::handle(event);
  default:
    return this->Fl_Browser_::handle(event);
  }
  return 1;
}

