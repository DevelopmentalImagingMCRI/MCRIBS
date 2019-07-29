# ============================================================================
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2013-2015 Imperial College London
# Copyright 2013-2015 Andreas Schuh
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
# @brief Build configuration of MIRTK Deformable library.
##############################################################################

set(HEADERS
  ${BINARY_INCLUDE_DIR}/mirtk/DeformableExport.h
  BalloonForce.h
  CurvatureConstraint.h
  DeformableConfig.h
  DeformableSurfaceDebugger.h
  DeformableSurfaceLogger.h
  DeformableSurfaceModel.h
  EulerMethod.h
  EulerMethodWithDamping.h
  EulerMethodWithMomentum.h
  ExternalForce.h
  ExternalForceTerm.h
  GaussCurvatureConstraint.h
  ImageEdgeDistance.h
  ImageEdgeForce.h
  ImplicitSurfaceForce.h
  ImplicitSurfaceDistance.h
  InflationForce.h
  InflationStoppingCriterion.h
  InternalForce.h
  InternalForceTerm.h
  MaximumCurvatureConstraint.h
  MeanCurvatureConstraint.h
  MetricDistortion.h
  MinActiveStoppingCriterion.h
  NonSelfIntersectionConstraint.h
  NormalForce.h
  PointSetForce.h
  QuadraticCurvatureConstraint.h
  RepulsiveForce.h
  SpringForce.h
  StretchingForce.h
  SurfaceConstraint.h
  SurfaceForce.h
)

set(SOURCES
  BalloonForce.cc
  CurvatureConstraint.cc
  DeformableConfig.cc
  DeformableSurfaceDebugger.cc
  DeformableSurfaceLogger.cc
  DeformableSurfaceModel.cc
  EulerMethod.cc
  EulerMethodWithDamping.cc
  EulerMethodWithMomentum.cc
  ExternalForce.cc
  GaussCurvatureConstraint.cc
  ImageEdgeDistance.cc
  ImageEdgeForce.cc
  ImplicitSurfaceForce.cc
  ImplicitSurfaceDistance.cc
  InflationForce.cc
  InflationStoppingCriterion.cc
  InternalForce.cc
  MaximumCurvatureConstraint.cc
  MeanCurvatureConstraint.cc
  MetricDistortion.cc
  MinActiveStoppingCriterion.cc
  NonSelfIntersectionConstraint.cc
  NormalForce.cc
  PointSetForce.cc
  QuadraticCurvatureConstraint.cc
  RepulsiveForce.cc
  SpringForce.cc
  StretchingForce.cc
  SurfaceConstraint.cc
  SurfaceForce.cc
)

set(DEPENDS
  LibCommon
  LibNumerics
  LibImage
  LibPointSet
  LibIO
  ${VTK_LIBRARIES}
)

if (TARGET TBB::tbb)
  list(APPEND DEPENDS TBB::tbb)
  if (TARGET TBB::malloc)
    list(APPEND DEPENDS TBB::malloc)
  endif ()
endif ()

mirtk_add_library(AUTO_REGISTER)
