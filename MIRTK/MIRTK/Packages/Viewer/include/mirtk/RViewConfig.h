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

#ifndef _RVIEWCONFIG_H
#define _RVIEWCONFIG_H

#include <mirtk/ViewerExport.h>


struct MIRTK_Viewer_EXPORT RViewConfig
{
  double xmin;
  double ymin;
  double xmax;
  double ymax;
  ViewerMode mode;
};

enum ConfigViewerMode
{
  _View_XY, _View_XZ, _View_YZ,
  _View_XY_XZ_v, _View_XY_YZ_v, _View_XZ_YZ_v,
  _View_XY_XZ_h, _View_XY_YZ_h, _View_XZ_YZ_h,
  _View_XY_XZ_YZ,
  _View_AB_XY_v, _View_AB_XZ_v, _View_AB_YZ_v,
  _View_AB_XY_h, _View_AB_XZ_h, _View_AB_YZ_h,
  _View_AB_XY_XZ_v, _View_AB_XY_XZ_h
};


/// RView configuration with single reslice plane in the X-Y plane
MIRTK_Viewer_EXPORT extern RViewConfig View_XY[];

/// RView configuration with single reslice plane in the X-Z plane
MIRTK_Viewer_EXPORT extern RViewConfig View_XZ[];

/// RView configuration with single reslice plane in the Y-Z plane
MIRTK_Viewer_EXPORT extern RViewConfig View_YZ[];

/// RView configuration with two vertical reslice planes (X-Y and X-Z plane).
MIRTK_Viewer_EXPORT extern RViewConfig View_XY_XZ_v[];

/// RView configuration with two vertical reslice planes (X-Y and Y-Z plane).
MIRTK_Viewer_EXPORT extern RViewConfig View_XY_YZ_v[];

/// RView configuration with two vertical reslice planes (X-Z and Y-Z plane).
MIRTK_Viewer_EXPORT extern RViewConfig View_XZ_YZ_v[];

/// RView configuration with two horizontal reslice planes (X-Y and Z-X plane).
MIRTK_Viewer_EXPORT extern RViewConfig View_XY_XZ_h[];

/// RView configuration with two horizontal reslice planes (X-Y and Y-Z plane).
MIRTK_Viewer_EXPORT extern RViewConfig View_XY_YZ_h[];

/// RView configuration with two horizontal reslice planes (X-Z and Y-Z plane).
MIRTK_Viewer_EXPORT extern RViewConfig View_XZ_YZ_h[];

/// RView configuration with three reslice planes
MIRTK_Viewer_EXPORT extern RViewConfig View_XY_XZ_YZ[];

/// RView configuration with single reslice plane for each image,
/// where the reslice planes of the target are left of those of the source
MIRTK_Viewer_EXPORT extern RViewConfig View_AB_XY_v[];

/// RView configuration with single reslice plane for each image,
/// where the reslice planes of the target are left of those of the source
MIRTK_Viewer_EXPORT extern RViewConfig View_AB_XZ_v[];

/// RView configuration with single reslice plane for each image,
/// where the reslice planes of the target are left of those of the source
MIRTK_Viewer_EXPORT extern RViewConfig View_AB_YZ_v[];

/// RView configuration with two horizontal reslice planes for both images
/// where the reslice planes of the target are left of those of the source
MIRTK_Viewer_EXPORT extern RViewConfig View_AB_XY_XZ_v[];

/// RView configuration with single reslice plane for each image,
/// where the reslice planes of the target are above those of the source
MIRTK_Viewer_EXPORT extern RViewConfig View_AB_XY_h[];

/// RView configuration with single reslice plane for each image,
/// where the reslice planes of the target are above those of the source
MIRTK_Viewer_EXPORT extern RViewConfig View_AB_XZ_h[];

/// RView configuration with single reslice plane for each image,
/// where the reslice planes of the target are above those of the source
MIRTK_Viewer_EXPORT extern RViewConfig View_AB_YZ_h[];

/// RView configuration with two horizontal reslice planes for both images
/// where the reslice planes of the target are above those of the source
MIRTK_Viewer_EXPORT extern RViewConfig View_AB_XY_XZ_h[];


#endif
