###########################################
## For building with the rest of OpenSim ##
###########################################

cmake_minimum_required(VERSION 3.2)

# Define project
PROJECT (Millard12EqWithAff)

# INCLUDE_DIRECTORIES(${OpenSim_SOURCE_DIR} ${OpenSim_SOURCE_DIR}/Vendors)

SET(OPENSIM_INSTALL_DIR $ENV{OPENSIM_HOME} 
		CACHE PATH "Top-level directory of OpenSim install")

# Change name of build target
SET(TARGET Millard12EqWithAff CACHE TYPE STRING)
# Identify the cpp file(s) that were to be built
FILE(GLOB SOURCE_FILES *.h *.cpp)
SET(SOURCE ${SOURCE_FILES})

# To add Debug feature add ";Debug" after Release on the line below
SET(CMAKE_CONFIGURATION_TYPES "RelWithDebInfo;Release" 
      CACHE STRING "Semicolon separated list of supported configuration types, only supports Debug, Release, MinSizeRel, and RelWithDebInfo, anything else will be ignored." FORCE )

# Location of headers
# SET(SIMTK_HEADERS_DIR ${OPENSIM_INSTALL_DIR}/sdk/include/SimTK/include)
SET(SIMTKcommon_HEADERS_DIR ${OPENSIM_INSTALL_DIR}/sdk/Simbody/include/)
SET(SIMBODY_HEADERS_DIR ${OPENSIM_INSTALL_DIR}/sdk/Simbody/include)
SET(SPDLOG_HEADERS_DIR ${OPENSIM_INSTALL_DIR}/sdk/spdlog/include)
SET(OPENSIM_HEADERS_DIR ${OPENSIM_INSTALL_DIR}/sdk/include)
# INCLUDE_DIRECTORIES(${SIMTK_HEADERS_DIR} ${SIMTKcommon_HEADERS_DIR} ${OPENSIM_HEADERS_DIR})
INCLUDE_DIRECTORIES(${SIMBODY_HEADERS_DIR} ${SPDLOG_HEADERS_DIR} ${OPENSIM_HEADERS_DIR})
# Libraries and dlls
SET(OPENSIM_LIBS_DIR ${OPENSIM_INSTALL_DIR}/sdk/lib ${OPENSIM_INSTALL_DIR}/lib ${OPENSIM_INSTALL_DIR}/sdk/Simbody/lib ${OPENSIM_INSTALL_DIR}/sdk/spdlog/lib)
SET(OPENSIM_DLLS_DIR ${OPENSIM_INSTALL_DIR}/bin ${OPENSIM_INSTALL_DIR}/sdk/Simbody/bin)
LINK_DIRECTORIES(${OPENSIM_LIBS_DIR} ${OPENSIM_DLLS_DIR})

ADD_EXECUTABLE(${TARGET} ${SOURCE})

TARGET_LINK_LIBRARIES(${TARGET}
	debug osimCommon_d optimized osimCommon
	debug osimSimulation_d optimized osimSimulation
	debug osimAnalyses_d optimized osimAnalyses
	debug osimActuators_d optimized osimActuators
	debug osimTools_d optimized osimTools
	debug SimTKcommon_d optimized SimTKcommon
	debug SimTKmath_d optimized SimTKmath
	debug SimTKsimbody_d optimized SimTKsimbody
)

# This block copies the additional files into the running directory
# For example vtp, obj files. Add to the end for more extentions
FILE(GLOB DATA_FILES *.vtp *.obj *.osim)
FOREACH (dataFile ${DATA_FILES})
 ADD_CUSTOM_COMMAND(
    TARGET ${TARGET}
    COMMAND ${CMAKE_COMMAND}
    ARGS -E copy
    ${dataFile}
	${Millard12EqWithAff_BINARY_DIR})
 ENDFOREACH (dataFile) 

MARK_AS_ADVANCED(CMAKE_INSTALL_PREFIX)
MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH)
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)
