cmake_minimum_required( VERSION 3.6 FATAL_ERROR )
set( CMAKE_VERBOSE_MAKEFILE ON )

set( CMAKE_CXX_STANDARD 17 )

if ( WIN32 )
	add_definitions(/bigobj)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /bigobj")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /bigobj")
endif ( WIN32 )

set( CMAKE_BUILD_TYPE Release )

project( SmpApp )

option(WITH_MATLAB "Use matlab engine" OFF)
option(WITH_BOOST "Use boost" OFF)
option(WITH_GEOGRAM "Use Geogram remesher" OFF)

include (${CMAKE_CURRENT_SOURCE_DIR}/CMakeLocal.txt)

set(CEREAL_INCLUDE_DIR "external/cereal/include")
set(EIGEN3_INCLUDE_DIR "external/eigen")
set(PHYSIM_INCLUDE_DIR "external/physim/include")
set(CETM_DIR "external/cetm")

if ( WIN32 )
	set ( LIB_PREFIX "" )
	set ( LIB_SUFFFIX "lib" )
    
    if ( WITH_MATLAB )
        if (NOT MATLAB_ENGINE)
            set(MATLAB_ENGINE "/extern/lib/win64/microsoft/libMatlabEngine.lib")
        endif ()
        if (NOT MATLAB_DATA)
            set(MATLAB_DATA "/extern/lib/win64/microsoft/libMatlabDataArray.lib")
        endif ()
    endif ( WITH_MATLAB)
endif ( WIN32 )

if ( UNIX )
	set ( LIB_PREFIX "lib" )
	set ( LIB_SUFFFIX "a" )
    
    if ( WITH_MATLAB )
        if (NOT MATLAB_ENGINE)
            set(MATLAB_ENGINE "/extern/bin/maci64/libMatlabEngine.dylib")
        endif ()
        if (NOT MATLAB_DATA)        
            set(MATLAB_DATA "/extern/bin/maci64/libMatlabDataArray.dylib")
        endif ()
    endif ( WITH_MATLAB)
endif ( UNIX )

get_filename_component( CINDER_PATH ${CINDER_DIR} ABSOLUTE )
get_filename_component( APP_PATH ${CMAKE_CURRENT_SOURCE_DIR} ABSOLUTE )

include( "${CINDER_PATH}/proj/cmake/modules/cinderMakeApp.cmake" )

if (WITH_MATLAB)
	add_definitions(-DUSE_MATLAB)
endif()

if (WITH_GEOGRAM)
	add_definitions(-DUSE_GEOGRAM)
endif()

file ( GLOB HEADERS
	${APP_PATH}/smp_addons.h
	${APP_PATH}/specimen.h
	${APP_PATH}/interface.h
	${APP_PATH}/smpup.h
	${APP_PATH}/smp_mesh.h
	${APP_PATH}/smp_physim.h
	${APP_PATH}/smp_simulator.h
	${APP_PATH}/bracket_model.h
	${APP_PATH}/bracket_configurer.h
	${APP_PATH}/utils/utils.h
	${APP_PATH}/utils/umath.h
	${APP_PATH}/utils/single_task_thread.h
	${APP_PATH}/utils/self_serializer.h
	${APP_PATH}/smpapp/smpapp.h
	${APP_PATH}/smpapp/smpapp_visual.h
	${APP_PATH}/smpapp/smpup_visual.h
	${APP_PATH}/smpapp/smp_imgui.h
	${APP_PATH}/smpapp/smp_imgui_handler.h
	${APP_PATH}/smpapp/smp_simulation_handler.h
)

file( GLOB SOURCES
	${APP_PATH}/specimen.cpp
	${APP_PATH}/smpup.cpp
	${APP_PATH}/smp_mesh.cpp
	${APP_PATH}/smp_physim.cpp
	${APP_PATH}/smp_simulator.cpp
	${APP_PATH}/bracket_model.cpp
	${APP_PATH}/bracket_configurer.cpp
	${APP_PATH}/utils/utils.cpp
	${APP_PATH}/utils/umath.cpp
	${APP_PATH}/smpapp/smpapp.cpp
	${APP_PATH}/smpapp/smpapp_visual.cpp
	${APP_PATH}/smpapp/smpup_visual.cpp
	${APP_PATH}/smpapp/smp_imgui.cpp
	${APP_PATH}/smpapp/smp_imgui_popups.cpp
	${APP_PATH}/smpapp/smp_imgui_main_menu.cpp
	${APP_PATH}/smpapp/smp_imgui_windows.cpp
	${APP_PATH}/smpapp/smp_imgui_handler.cpp
	${APP_PATH}/smpapp/smp_simulation_handler.cpp
)

set ( INCLUDES
	${EIGEN3_INCLUDE_DIR}
	${LIBIGL_DIR}/include
	${LIBIGL_DIR}/include/igl/copyleft
	${LIBIGL_DIR}/external/triangle
    ${CETM_DIR}/include
    ${PHYSIM_INCLUDE_DIR}
	${CEREAL_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
)

if (WITH_MATLAB)
	file( GLOB HEADERS
		${HEADERS}
		utils/smp_matlab.h
	)
	file( GLOB SOURCES
		${SOURCES}
		utils/smp_matlab.cpp
	)
	set ( INCLUDES
		${INCLUDES}
		${MATLAB_DIR}/extern/include
	)
	set ( LIBRARIES
		${LIBRARIES}
		${MATLAB_DIR}${MATLAB_ENGINE}
		${MATLAB_DIR}${MATLAB_DATA}
	)
endif()

if (WITH_GEOGRAM)
	file( GLOB HEADERS
		${HEADERS}
		utils/smp_geogram.h
		${APP_PATH}/smpapp/smp_remesh_handler.h
	)
	file( GLOB SOURCES
		${SOURCES}
		utils/smp_geogram.cpp
		${APP_PATH}/smpapp/smp_remesh_handler.cpp
	)
	set ( INCLUDES
		${INCLUDES}
		${GEOGRAM_INCLUDE_DIR}
	)
	set ( LIBRARIES
		${LIBRARIES}
	)
endif()

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} PREFIX Header FILES ${HEADERS})
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} PREFIX Source FILES ${SOURCES})

set( SOURCES ${HEADERS} ${SOURCES})

add_subdirectory(external/PhySim "${CMAKE_CURRENT_BINARY_DIR}/PhySim_build")
add_subdirectory(external/cetm "${CMAKE_CURRENT_BINARY_DIR}/cetm_build")

ci_make_app(
    APP_NAME    "SmpApp"
    CINDER_PATH ${CINDER_PATH}
    SOURCES     ${SOURCES}
    INCLUDES    ${APP_PATH}/include ${INCLUDES}
	LIBRARIES   ${LIBRARIES}
)

target_compile_definitions(SmpApp PRIVATE
	EIGEN_DENSEBASE_PLUGIN="eigen_dense_base_addons.h"
)

add_dependencies(${PROJECT_NAME} PhySim)
add_dependencies(${PROJECT_NAME} cetm)

target_link_libraries( SmpApp PRIVATE PhySim )
target_link_libraries( SmpApp PRIVATE cetm )

target_link_libraries( SmpApp PRIVATE debug ${LIBIGL_BUILD_DIR}/lib/Debug/${LIB_PREFIX}triangle.${LIB_SUFFFIX} )
target_link_libraries( SmpApp PRIVATE optimized ${LIBIGL_BUILD_DIR}/lib/Release/${LIB_PREFIX}triangle.${LIB_SUFFFIX} )

if ( WITH_GEOGRAM )
	target_link_libraries( SmpApp PRIVATE debug ${GEOGRAM_LIB_DEBUG}/${LIB_PREFIX}geogram.${LIB_SUFFFIX})
	target_link_libraries( SmpApp PRIVATE optimized ${GEOGRAM_LIB_RELEASE}/${LIB_PREFIX}geogram.${LIB_SUFFFIX})
endif()

# Boost
if ( WITH_BOOST )
	add_definitions(-DUSE_BOOST)
    
    set(boost_min_ver 1.50.0)
    set(boost_libs system filesystem)
    find_package(Boost ${boost_min_ver})

    if(Boost_FOUND)
        find_package(Boost ${boost_min_ver} COMPONENTS ${boost_libs})
    endif()

    target_link_libraries(SmpApp PRIVATE ${Boost_LIBRARIES})
endif ( WITH_BOOST )
