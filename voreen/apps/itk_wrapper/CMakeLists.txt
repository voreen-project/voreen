PROJECT(itk_wrapper)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0 FATAL_ERROR)
INCLUDE(../../cmake/commonconf.cmake)

MESSAGE(STATUS "Configuring ITK_Wrapper Application")

SET(VRN_ITKWRAPPER_SOURCES
    itk_wrapper.cpp
    baseclasses.cpp
    template.cpp
)

SET(VRN_ITKWRAPPER_HEADERS
    itk_wrapper.h
    baseclasses.h
    template.h
)

ADD_EXECUTABLE(itk_wrapper 
    ${VRN_ITKWRAPPER_SOURCES} ${VRN_ITKWRAPPER_HEADERS} 
) 
ADD_DEFINITIONS(${VRN_DEFINITIONS} ${VRN_MODULE_DEFINITIONS})
INCLUDE_DIRECTORIES(${VRN_INCLUDE_DIRECTORIES})
TARGET_LINK_LIBRARIES(itk_wrapper tgt voreen_core ${VRN_EXTERNAL_LIBRARIES})
