# ============================================================================
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright Imperial College London
# Copyright 2018 Andreas Schuh
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
# @brief Build configuration of MIRTK Viewer library.
##############################################################################

set(HEADERS
  ${BINARY_INCLUDE_DIR}/mirtk/ViewerExport.h
  Color.h
  ColorRGBA.h
  Contour.h
  LookupTable.h
  RView.h
  RViewConfig.h
  Viewer.h
  HistogramWindow.h
  Segment.h
  SegmentTable.h
  VoxelContour.h
)

set(SOURCES
  Color.cc
  ColorRGBA.cc
  LookupTable.cc
  RView.cc
  RViewConfig.cc
  Viewer.cc
  HistogramWindow.cc
  Segment.cc
  SegmentTable.cc
  VoxelContour.cc
)

set(DEPENDS
  LibCommon
  LibNumerics
  LibImage
  LibPointSet
  LibIO
  LibTransformation
  ${OPENGL_LIBRARIES}
  ${GLUT_LIBRARIES}
  ${VTK_LIBRARIES}
)

mirtk_add_library()
