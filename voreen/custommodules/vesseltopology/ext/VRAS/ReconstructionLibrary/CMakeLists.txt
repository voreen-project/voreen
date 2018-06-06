
cmake_minimum_required(VERSION 3.5)
project(VascResc CXX)

set ( VasRec_VERSION_MAJOR 0)
set ( VasRec_VERSION_MINOR 1)

set(CORE_DIR "${PROJECT_SOURCE_DIR}/Core") 
set(MATHCORE_DIR "${PROJECT_SOURCE_DIR}/MathCore")
set(IO_DIR "${PROJECT_SOURCE_DIR}/IO")
set(RECONSTRUCT2D_DIR "${PROJECT_SOURCE_DIR}/Reconstruct2D")
set(RECONSTRUCT3D_DIR "${PROJECT_SOURCE_DIR}/Reconstruct3D")

find_package(GSL REQUIRED)

file( GLOB_RECURSE Core_Headers "${CORE_DIR}/*.h")
set(CORE_SOURCES ${Core_Headers})


file( GLOB_RECURSE Math_Headers "${MATHCORE_DIR}/*.h")
set(MATH_SOURCES ${Math_Headers})

file( GLOB_RECURSE IO_Headers "${IO_DIR}/*.h")
file( GLOB_RECURSE IO_Impl "${IO_DIR}/*.cxx")
set(IO_SOURCES ${IO_Impl} ${IO_Headers})


file( GLOB_RECURSE RECONSTRUCT2D_Headers  "${RECONSTRUCT2D_DIR}/*.h")
set(RECONSTRUCT2D_SOURCES ${RECONSTRUCT2D_Headers} )


message(STATUS "2D SOURCES >> ${RECONSTRUCT2D_SOURCES}")
file( GLOB_RECURSE RECONSTRUCT3D_Impl  "${RECONSTRUCT3D_DIR}/*.cxx")

set(RECONSTRUCT3D_SOURCES ${RECONSTRUCT3D_Headers}  ${RECONSTRUCT3D_Impl})


find_package( ITK REQUIRED)
include(${ITK_USE_FILE})

find_package( VTK REQUIRED)
include(${VTK_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR}/Core)
include_directories(MathCore)
include_directories(IO)
include_directories(RECONSTRUCT2D_DIR)
include_directories(RECONSTRUCT3D_DIR)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11  -g  -O3 -Wall  ")

add_library( VascResc STATIC ${CORE_SOURCES} ${MATH_SOURCES} ${IO_SOURCES} ${RECONSTRUCT2D_SOURCES} ${RECONSTRUCT3D_SOURCES} )

target_link_libraries( VascResc ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${GSL_LIBRARIES} )

set_target_properties(VascResc PROPERTIES LINKER_LANGUAGE CXX)
