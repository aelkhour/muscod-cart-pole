# File     CMakeLists.txt
# Author   Christian Hoffmann, Antonio El Khoury
# Date     2009, 2012
#
# This file is part of the MUSCOD_TEST package. MUSCOD_TEST is proprietary software of
#   Simulation and Optimization Workgroup
#   Interdisciplinary Center for Scientific Computing (IWR)
#   University of Heidelberg, Germany.
#
####################################################################################################
#
# Main CMake script for building the MUSCOD_SAMPLE optimal control problem
#
####################################################################################################

####################################################################################################
#### BASIC CONFIGURATION
####################################################################################################
CMAKE_MINIMUM_REQUIRED ( VERSION 2.6 )

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/cpack.cmake)

PROJECT( "muscod_sample_application" CXX )

LIST( APPEND CMAKE_MODULE_PATH
	${PROJECT_SOURCE_DIR}/CMake
	${PROJECT_SOURCE_DIR}/CMake/Modules
)

INCLUDE( SimOptDefaults )

####################################################################################################
#### GET USED PACKAGES
####################################################################################################
FIND_PACKAGE( MUSCOD REQUIRED )  # find and include "MUSCODConfig.cmake"
IF( MUSCOD_FOUND )
	INCLUDE( ${MUSCOD_USE_FILE} )
ENDIF( MUSCOD_FOUND )

####################################################################################################
#### BUILD LIBRARY OF THE OPTIMAL CONTROL PROBLEM
####################################################################################################

# NOTE: Here you need to add all source files that are required for your
#       optimal control problem.
ADD_LIBRARY (cart-pole
  SRC/cart-pole.cc
  )

####################################################################################################
#### SELF-TESTING
####################################################################################################
INCLUDE( TESTDartConfig ) # NOTE: line must be called BEFORE INCLUDE( Dart )!
INCLUDE( Dart )       # Makes CTest and Dart available, including call of ENABLE_TESTING()
INCLUDE( TestSuite )   # Defines the actual tests.
