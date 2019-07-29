/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <mirtk/RView.h>

#include <mirtk/Image.h>
#include <mirtk/Transformation.h>
#include <mirtk/Registration.h>

#ifdef __APPLE__
#include <OpenGl/gl.h>
#include <OpenGl/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif



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

