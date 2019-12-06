cmake_minimum_required(VERSION 3.6)

set(EIGEN3_INC_DIR "../eigen")
set(LINKSTATIC true)

project(PhySim)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include_directories(${EIGEN3_INC_DIR})

include_directories(${SS_SPQR_INC_DIR})
include_directories(${CHOLMOD_INC_DIR})
include_directories(${SSPARSE_INC_DIR})

# Add project files

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB SRCFILES 
	source/PhySim/*.cpp 
	source/PhySim/Energies/*.cpp
	source/PhySim/Geometry/*.cpp
	source/PhySim/Models/*.cpp
	source/PhySim/Solvers/*.cpp
	source/PhySim/Utils/*.cpp
	)
	
file(GLOB INCFILES 
	include/PhySim/*.h 
	include/PhySim/Energies/*.h
	include/PhySim/Geometry/*.h
	include/PhySim/Models/*.h
	include/PhySim/Solvers/*.h
	include/PhySim/Utils/*.h
	)

source_group(TREE ${PROJECT_SOURCE_DIR}/source/PhySim PREFIX Source FILES ${SRCFILES})
source_group(TREE ${PROJECT_SOURCE_DIR}/include/PhySim PREFIX Header FILES ${INCFILES})

# Create library

set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -openmp" )

if (LINKSTATIC)
set( CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -MTd" )
set( CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -MT" )
else (LINKSTATIC)
set( CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -MDd" )
set( CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -MD" )
endif (LINKSTATIC)

add_library(${PROJECT_NAME} STATIC ${SRCFILES} ${INCFILES})