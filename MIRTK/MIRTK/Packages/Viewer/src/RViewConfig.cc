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

#include <mirtk/OpenGl.h>
#include <mirtk/RView.h>

#include <mirtk/Image.h>
#include <mirtk/Transformation.h>
#include <mirtk/Registration.h>


RViewConfig View_XY[] = {
  {  0.0, 0.0, 1.0, 1.0, Viewer_XY },
  { -1.0, 0.0, 0.0, 0.0, Viewer_None }
};

RViewConfig View_XZ[] = {
  {  0.0, 0.0, 1.0, 1.0, Viewer_XZ },
  { -1.0, 0.0, 0.0, 0.0, Viewer_None }
};

RViewConfig View_YZ[] = {
  {  0.0, 0.0, 1.0, 1.0, Viewer_YZ },
  { -1.0, 0.0, 0.0, 0.0, Viewer_None }
};

RViewConfig View_XY_XZ_v[] = {
  {  0.0, 0.5, 1.0, 1.0, Viewer_XY },
  {  0.0, 0.0, 1.0, 0.5, Viewer_XZ },
  { -1.0, 0.0, 0.0, 0.0, Viewer_None }
};

RViewConfig View_XY_XZ_h[] = {
  {  0.0, 0.0, 0.5, 1.0, Viewer_XY },
  {  0.5, 0.0, 1.0, 1.0, Viewer_XZ },
  { -1.0, 0.0, 0.0, 0.0, Viewer_None }
};

RViewConfig View_XY_YZ_v[] = {
  {  0.0, 0.5, 1.0, 1.0, Viewer_XY },
  {  0.0, 0.0, 1.0, 0.5, Viewer_YZ },
  { -1.0, 0.0, 0.0, 0.0, Viewer_None }
};

RViewConfig View_XY_YZ_h[] = {
  {  0.0, 0.0, 0.5, 1.0, Viewer_XY },
  {  0.5, 0.0, 1.0, 1.0, Viewer_YZ },
  { -1.0, 0.0, 0.0, 0.0, Viewer_None }
};

RViewConfig View_XZ_YZ_v[] = {
  {  0.0, 0.5, 1.0, 1.0, Viewer_XZ },
  {  0.0, 0.0, 1.0, 0.5, Viewer_YZ },
  { -1.0, 0.0, 0.0, 0.0, Viewer_None }

};

RViewConfig View_XZ_YZ_h[] = {
  {  0.0, 0.0, 0.5, 1.0, Viewer_XZ },
  {  0.5, 0.0, 1.0, 1.0, Viewer_YZ },
  { -1.0, 0.0, 0.0, 0.0, Viewer_None }
};

RViewConfig View_XY_XZ_YZ[] = {
  {  0.0, 0.5, 0.5, 1.0, Viewer_XY },
  {  0.5, 0.5, 1.0, 1.0, Viewer_XZ },
  {  0.0, 0.0, 0.5, 0.5, Viewer_YZ },
  { -1.0, 0.0, 0.0, 0.0, Viewer_None }
};

RViewConfig View_AB_XY_v[] = {
  {  0.0, 0.0, 0.5, 1.0, Viewer_XY },
  {  0.5, 0.0, 1.0, 1.0, Viewer_XY },
  { -1.0, 0.0, 0.0, 0.0, Viewer_None }
};

RViewConfig View_AB_XZ_v[] = {
  {  0.0, 0.0, 0.5, 1.0, Viewer_XZ },
  {  0.5, 0.0, 1.0, 1.0, Viewer_XZ },
  { -1.0, 0.0, 0.0, 0.0, Viewer_None }
};

RViewConfig View_AB_YZ_v[] = {
  {  0.0, 0.0, 0.5, 1.0, Viewer_YZ },
  {  0.5, 0.0, 1.0, 1.0, Viewer_YZ },
  { -1.0, 0.0, 0.0, 0.0, Viewer_None }
};

RViewConfig View_AB_XY_XZ_v[] = {
  {  0.0, 0.5, 0.5, 1.0, Viewer_XY },
  {  0.5, 0.5, 1.0, 1.0, Viewer_XY },
  {  0.0, 0.0, 0.5, 0.5, Viewer_XZ },
  {  0.5, 0.0, 1.0, 0.5, Viewer_XZ },
  { -1.0, 0.0, 0.0, 0.0, Viewer_None }
};

RViewConfig View_AB_XY_h[] = {
  {  0.0, 0.0, 1.0, 0.5, Viewer_XY },
  {  0.0, 0.5, 1.0, 1.0, Viewer_XY },
  { -1.0, 0.0, 0.0, 0.0, Viewer_None }
};

RViewConfig View_AB_XZ_h[] = {
  {  0.0, 0.0, 1.0, 0.5, Viewer_XZ },
  {  0.0, 0.5, 1.0, 1.0, Viewer_XZ },
  { -1.0, 0.0, 0.0, 0.0, Viewer_None }
};

RViewConfig View_AB_YZ_h[] = {
  {  0.0, 0.0, 1.0, 0.5, Viewer_YZ },
  {  0.0, 0.5, 1.0, 1.0, Viewer_YZ },
  { -1.0, 0.0, 0.0, 0.0, Viewer_None }
};

RViewConfig View_AB_XY_XZ_h[] = {
  {  0.0, 0.5, 0.5, 1.0, Viewer_XY },
  {  0.5, 0.5, 1.0, 1.0, Viewer_XZ },
  {  0.0, 0.0, 0.5, 0.5, Viewer_XY },
  {  0.5, 0.0, 1.0, 0.5, Viewer_XZ },
  { -1.0, 0.0, 0.0, 0.0, Viewer_None }
};

