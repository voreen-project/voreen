CMAKE_MINIMUM_REQUIRED(VERSION 3.16.3 FATAL_ERROR)
PROJECT(networkevaluatortest)
INCLUDE(../../../cmake/commonconf.cmake)

MESSAGE(STATUS "Configuring NetworkEvaluatorTest Application")

ADD_EXECUTABLE(networkevaluatortest networkevaluatortest.cpp)
ADD_DEFINITIONS(${VRN_DEFINITIONS} ${VRN_MODULE_DEFINITIONS})
INCLUDE_DIRECTORIES(${VRN_INCLUDE_DIRECTORIES})
TARGET_LINK_LIBRARIES(networkevaluatortest tgt voreen_core ${VRN_EXTERNAL_LIBRARIES} )

