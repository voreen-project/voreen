PROJECT(volumeorigintest)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0 FATAL_ERROR)
INCLUDE(../../../cmake/commonconf.cmake)

MESSAGE(STATUS "Configuring VolumeOriginTest Application")

ADD_EXECUTABLE(volumeorigintest volumeorigintest.cpp)
ADD_DEFINITIONS(${VRN_DEFINITIONS} ${VRN_MODULE_DEFINITIONS})
INCLUDE_DIRECTORIES(${VRN_INCLUDE_DIRECTORIES})
TARGET_LINK_LIBRARIES(volumeorigintest tgt voreen_core ${VRN_EXTERNAL_LIBRARIES} )

