CMAKE_MINIMUM_REQUIRED(VERSION 3.16.3 FATAL_ERROR)
PROJECT(processorcreatetest)
INCLUDE(../../../cmake/commonconf.cmake)

MESSAGE(STATUS "Configuring ProcessorCreateTest Application")

ADD_EXECUTABLE(processorcreatetest processorcreatetest.cpp)
ADD_DEFINITIONS(${VRN_DEFINITIONS} ${VRN_MODULE_DEFINITIONS})
INCLUDE_DIRECTORIES(${VRN_INCLUDE_DIRECTORIES})
TARGET_LINK_LIBRARIES(processorcreatetest tgt voreen_core ${VRN_EXTERNAL_LIBRARIES} )

