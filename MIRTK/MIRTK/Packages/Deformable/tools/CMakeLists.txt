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
# @brief Build configuration of MIRTK Deformable tools.
##############################################################################

mirtk_add_executable(
  deform-mesh
  DEPENDS
    LibCommon
    LibNumerics
    LibImage
    LibIO
    LibPointSet
    LibTransformation
    LibDeformable
    ${VTK_LIBRARIES}
)

mirtk_add_executable(recon-neonatal-cortex DEPENDS ${BASIS_PYTHON_LIBRARY_TARGET})
