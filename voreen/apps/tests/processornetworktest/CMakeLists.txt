PROJECT(processornetworktest)
CMAKE_MINIMUM_REQUIRED(VERSION 3.5.1 FATAL_ERROR)
INCLUDE(../../../cmake/commonconf.cmake)

MESSAGE(STATUS "Configuring ProcessorNetworkTest Application")

ADD_EXECUTABLE(processornetworktest processornetworktest.cpp)
ADD_DEFINITIONS(${VRN_DEFINITIONS} ${VRN_MODULE_DEFINITIONS})
INCLUDE_DIRECTORIES(${VRN_INCLUDE_DIRECTORIES})
TARGET_LINK_LIBRARIES(processornetworktest tgt voreen_core ${VRN_EXTERNAL_LIBRARIES} )

