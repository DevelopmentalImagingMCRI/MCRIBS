# ============================================================================
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2008-2017 Imperial College London
# Copyright 2016-2018 Andreas Schuh
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ============================================================================

##############################################################################
# @file  CMakeLists.txt
# @brief Build configuration of MIRTK Viewer GUI.
##############################################################################

mirtk_add_executable(
  view
  SOURCES
    Fl_RView.cc
    Fl_RView.h
    Fl_RViewUI.cc
    Fl_RViewUI.h
    Fl_RView_Histogram.cc
    Fl_RView_Histogram.h
    Fl_RView_Landmark.cc
    Fl_RView_Landmark.h
    Fl_RView_Registration.cc
    Fl_RView_Registration.h
    Fl_RView_Segmentation.cc
    Fl_RView_Segmentation.h
    Fl_RView_Transformation.cc
    Fl_RView_Transformation.h
    Fl_RView_Visualisation.cc
    Fl_RView_Visualisation.h
    Fl_Simple_Browser.cc
    Fl_Simple_Browser.h
    Fl_Simple_Group.cc
    Fl_Simple_Group.h
    Fl_Value_Slider2.cc
    Fl_Value_Slider2.h
  DEPENDS
    LibCommon
    LibNumerics
    LibImage
    LibIO
    LibPointSet
    LibTransformation
    LibViewer
    ${FLTK_LIBRARIES}
    ${GLUT_LIBRARIES}
    ${OPENGL_LIBRARIES}
    ${VTK_LIBRARIES}
)
