# ==============================================================================
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
# ==============================================================================

################################################################################
# @file  BasisProject.cmake
# @brief Sets basic information about the MIRTK module and calls basis_project().
#
# This file defines basic information about a project by calling 
# the basis_project() function. This basic information, also known as metadata, 
# is used by CMake BASIS to setup the project. The dependencies to other modules
# have to be specified here such that the top-level IRTK project can analyze the
# inter-module dependencies, as well as dependencies on third-party libraries.
#
# @sa https://cmake-basis.github.io/standard/modules.html
#
# @ingroup BasisSettings
################################################################################

# Note: The #<*> dependency patterns are required by the basisproject tool and
#       should be kept on a separate line as last commented argument of the
#       corresponding options of the basis_project() command. The TEMPLATE
#       option and set argument are also required by this tool and should not
#       be changed manually. The argument is updated by basisproject --update.

basis_project (

  # ----------------------------------------------------------------------------
  # meta-data
  NAME        "Mapping"
  VERSION     "0.0.0" # version of this module
  SOVERSION   "0"     # API yet unstable
  PACKAGE     "MIRTK"
  AUTHORS     "Andreas Schuh"
  DESCRIPTION "Brain mapping library of the Medical Image Registration ToolKit."
  COPYRIGHT   "2013-2016 Imperial College London, Andreas Schuh"
  LICENSE     "Apache License Version 2.0"
  CONTACT     "Andreas Schuh <andreas.schuh.84@gmail.com>"
  TEMPLATE    "mirtk-module/1.0"

  # ----------------------------------------------------------------------------
  # dependencies
  DEPENDS
    MIRTK{Common,Numerics,Image,PointSet,IO}
    Boost-1.48 # {math_c99} used, but headers only
    Eigen3
    "VTK-8|7|6{
      vtkCommonCore,
      vtkCommonDataModel,
      vtkFiltersCore,
      vtkFiltersHybrid,
      vtkFiltersModeling,
      vtkImagingCore,
      vtkIOXML
    }"
    #<dependency>
  OPTIONAL_DEPENDS
    TBB{tbb}
    #<optional-dependency>
  TEST_DEPENDS
    #<test-dependency>
  OPTIONAL_TEST_DEPENDS
    #<optional-test-dependency>

)
