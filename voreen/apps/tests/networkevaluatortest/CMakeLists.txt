PROJECT(networkevaluatortest)
CMAKE_MINIMUM_REQUIRED(VERSION 3.5.1 FATAL_ERROR)
INCLUDE(../../../cmake/commonconf.cmake)

MESSAGE(STATUS "Configuring NetworkEvaluatorTest Application")

ADD_EXECUTABLE(networkevaluatortest networkevaluatortest.cpp)
ADD_DEFINITIONS(${VRN_DEFINITIONS} ${VRN_MODULE_DEFINITIONS})
INCLUDE_DIRECTORIES(${VRN_INCLUDE_DIRECTORIES})
TARGET_LINK_LIBRARIES(networkevaluatortest tgt voreen_core ${VRN_EXTERNAL_LIBRARIES} )

