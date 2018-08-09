PROJECT(descriptiontest)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0 FATAL_ERROR)
INCLUDE(../../../cmake/commonconf.cmake)

MESSAGE(STATUS "Configuring DescriptionTest Application")

ADD_EXECUTABLE(descriptiontest descriptiontest.cpp)
ADD_DEFINITIONS(${VRN_DEFINITIONS} ${VRN_MODULE_DEFINITIONS})
INCLUDE_DIRECTORIES(${VRN_INCLUDE_DIRECTORIES})
TARGET_LINK_LIBRARIES(descriptiontest tgt voreen_core ${VRN_EXTERNAL_LIBRARIES} )

