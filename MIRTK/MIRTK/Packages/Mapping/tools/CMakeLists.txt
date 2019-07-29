# ============================================================================
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2015-2016 Imperial College London
# Copyright 2015-2016 Andreas Schuh
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
# @brief Build configuration of MIRTK Mapping tools.
##############################################################################

macro(add_mapping_tool cmd)
  mirtk_add_executable(${cmd} DEPENDS Lib${PROJECT_NAME} ${VTK_LIBRARIES})
endmacro ()

# Compute map
add_mapping_tool(calculate-boundary-map)
add_mapping_tool(calculate-surface-map)
add_mapping_tool(calculate-volume-map)
add_mapping_tool(compose-maps)

# Evaluate map
add_mapping_tool(evaluate-surface-map)
add_mapping_tool(evaluate-volume-map)

# Apply map
#add_mapping_tool(map-surface)
#add_mapping_tool(map-volume)
