/*=========================================================================

  Library   : Image Registration Toolkit ()
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _RVIEWCONFIG_H
#define _RVIEWCONFIG_H

typedef struct {

  double xmin;
  double ymin;
  double xmax;
  double ymax;
  ViewerMode mode;

} RViewConfig;

typedef enum { _View_XY, _View_XZ, _View_YZ,
               _View_XY_XZ_v, _View_XY_YZ_v, _View_XZ_YZ_v,
               _View_XY_XZ_h, _View_XY_YZ_h, _View_XZ_YZ_h,
               _View_XY_XZ_YZ,
               _View_AB_XY_v, _View_AB_XZ_v, _View_AB_YZ_v,
               _View_AB_XY_h, _View_AB_XZ_h, _View_AB_YZ_h,
               _View_AB_XY_XZ_v, _View_AB_XY_XZ_h
             } ConfigViewerMode;


/// RView configuration with single reslice plane in the X-Y plane
extern RViewConfig View_XY[];

/// RView configuration with single reslice plane in the X-Z plane
extern RViewConfig View_XZ[];

/// RView configuration with single reslice plane in the Y-Z plane
extern RViewConfig View_YZ[];

/// RView configuration with two vertical reslice planes (X-Y and X-Z plane).
extern RViewConfig View_XY_XZ_v[];

/// RView configuration with two vertical reslice planes (X-Y and Y-Z plane).
extern RViewConfig View_XY_YZ_v[];

/// RView configuration with two vertical reslice planes (X-Z and Y-Z plane).
extern RViewConfig View_XZ_YZ_v[];

/// RView configuration with two horizontal reslice planes (X-Y and Z-X plane).
extern RViewConfig View_XY_XZ_h[];

/// RView configuration with two horizontal reslice planes (X-Y and Y-Z plane).
extern RViewConfig View_XY_YZ_h[];

/// RView configuration with two horizontal reslice planes (X-Z and Y-Z plane).
extern RViewConfig View_XZ_YZ_h[];

/// RView configuration with three reslice planes
extern RViewConfig View_XY_XZ_YZ[];

/// RView configuration with single reslice plane for each image,
/// where the reslice planes of the target are left of those of the source
extern RViewConfig View_AB_XY_v[];

/// RView configuration with single reslice plane for each image,
/// where the reslice planes of the target are left of those of the source
extern RViewConfig View_AB_XZ_v[];

/// RView configuration with single reslice plane for each image,
/// where the reslice planes of the target are left of those of the source
extern RViewConfig View_AB_YZ_v[];

/// RView configuration with two horizontal reslice planes for both images
/// where the reslice planes of the target are left of those of the source
extern RViewConfig View_AB_XY_XZ_v[];

/// RView configuration with single reslice plane for each image,
/// where the reslice planes of the target are above those of the source
extern RViewConfig View_AB_XY_h[];

/// RView configuration with single reslice plane for each image,
/// where the reslice planes of the target are above those of the source
extern RViewConfig View_AB_XZ_h[];

/// RView configuration with single reslice plane for each image,
/// where the reslice planes of the target are above those of the source
extern RViewConfig View_AB_YZ_h[];

/// RView configuration with two horizontal reslice planes for both images
/// where the reslice planes of the target are above those of the source
extern RViewConfig View_AB_XY_XZ_h[];


#endif
