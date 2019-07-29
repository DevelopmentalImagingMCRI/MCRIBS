# ============================================================================
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2013-2016 Imperial College London
# Copyright 2013-2016 Andreas Schuh
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
# @brief Build configuration of MIRTK Mapping library.
##############################################################################

set(HEADERS ${BINARY_INCLUDE_DIR}/mirtk/${PROJECT_NAME}Export.h)
set(SOURCES)

# List of library classes where indentation reflects inheritance
# Note: Header only classes must be included with the .h file name extension.
set(CLASSES
  # Map types
  Mapping
    MeshlessMap
      MeshlessHarmonicMap
        MeshlessBiharmonicMap
    PiecewiseLinearMap
  # Surface boundary parameterization
  BoundarySegmentParameterizer
    UniformBoundarySegmentParameterizer
    ChordLengthBoundarySegmentParameterizer
    SubdividedBoundarySegmentParameterizer
  # Surface boundary mapping
  BoundaryMapper
    BoundarySegmentMapper
      BoundaryToDiskMapper
      BoundaryToPolygonMapper
      BoundaryToSquareMapper
  # Surface mapping
  SurfaceMapper
    FixedBoundarySurfaceMapper
      LinearFixedBoundarySurfaceMapper
        SymmetricWeightsSurfaceMapper
          UniformSurfaceMapper
          ChordLengthSurfaceMapper
          HarmonicSurfaceMapper
          AuthalicSurfaceMapper
          IntrinsicSurfaceMapper
            NearOptimalIntrinsicSurfaceMapper
              IntrinsicLeastAreaDistortionSurfaceMapper
              IntrinsicLeastEdgeLengthDistortionSurfaceMapper
        NonSymmetricWeightsSurfaceMapper
          ShapePreservingSurfaceMapper
          MeanValueSurfaceMapper
      #FlatSurfaceMapper
      #  RegularGridSurfaceMapper
      #    HarmonicRegularGridSurfaceMapper
    FreeBoundarySurfaceMapper
      LeastSquaresConformalSurfaceMapper
    SphericalSurfaceMapper
      ConformalSurfaceFlattening
    #  SphericalMultiDimensionalScaling
  # Volumetric mapping
  VolumeMapper
    TetrahedralMeshMapper
      LinearTetrahedralMeshMapper
        HarmonicTetrahedralMeshMapper
        AsConformalAsPossibleMapper
    MeshlessVolumeMapper
      MeshlessHarmonicVolumeMapper
)

if (MIRTK_Numerics_WITH_eigs)
  list(APPEND CLASSES SpectralConformalSurfaceMapper)
endif ()

# Add source files implementing each class to HEADERS and SOURCES lists
foreach (class IN LISTS CLASSES)
  if (class MATCHES "\\.h$")
    list(APPEND HEADERS ${class})
  else ()
    list(APPEND HEADERS ${class}.h)
    list(APPEND SOURCES ${class}.cc)
  endif ()
endforeach ()

# List of link dependencies
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
endif ()

# Add library target
mirtk_add_library()
