CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

##
## Uncomment the following lines according to
## the external libs that are used by the voreen_core lib.
##
#SET(VRN_WITH_DEVIL 1)
#SET(VRN_WITH_ZLIB 1)
#SET(VRN_WITH_TIFF 1)
#SET(VRN_WITH_FFMPEG 1)
#SET(VRN_WITH_PYTHON 1)
#SET(VRN_WITH_DCMTK 1)

FIND_PATH(VRN_DIR config.txt)
OPTION(VRN_WITH_DEVIL               "VRN_WITH_DEVIL"                ON)
OPTION(VRN_WITH_ZLIB                "VRN_WITH_ZLIB"                 ON)
OPTION(VRN_WITH_TIFF                "VRN_WITH_TIFF"                 OFF)
OPTION(VRN_WITH_FFMPEG              "VRN_WITH_FFMPEG"               OFF)
OPTION(VRN_WITH_PYTHON              "VRN_WITH_PYTHON"               OFF)
OPTION(VRN_WITH_DCMTK               "VRN_WITH_DCMTK"                OFF)

SET(CMAKE_EXE_LINKER_FLAGS /NODEFAULTLIB:libc.lib)

# Set the project/library name here.  Classes should be declared as
# "class vtkVoreenActor_EXPORT vtkFoo", where vtkVoreenActor is the name of the
# project set here.
PROJECT(vtkVoreenActor)

# Set your list of sources here.  Do not change the name of the
# vtkVoreenActor_SRCS variable.
SET(vtkVoreenActor_SRCS
  vtkVoreenActor.cpp
  vtkVoreenActor.h
  vtkdummycanvas.cpp
  vtkdummyCanvas.h
)

# List the kits from VTK that are needed by this project.
SET(vtkVoreenActor_LIBS
  vtkCommon
  vtkRendering 
  vtkGraphics 
  vtkHybrid 
  vtkImaging 
  vtkIO  
  vtkFiltering 
  vtkCommon
  vtkVolumeRendering

  glu32 	
  netapi32
  Version
  WbemUuid

  ${VRN_DIR}/$(ConfigurationName)/voreen_core.lib
  ${VRN_DIR}/ext/glew/lib/$(ConfigurationName)/glew32s.lib
)
  
IF(VRN_WITH_DEVIL)
  SET(vtkVoreenActor_LIBS ${vtkVoreenActor_LIBS} 
    ${VRN_DIR}/ext/il/lib/DevIL.lib
  )
ENDIF(VRN_WITH_DEVIL)

IF(VRN_WITH_ZLIB)
  SET(vtkVoreenActor_LIBS ${vtkVoreenActor_LIBS} 
    ${VRN_DIR}/ext/zlib/lib/zdll.lib
  )
ENDIF(VRN_WITH_ZLIB)

IF(VRN_WITH_TIFF)
  SET(vtkVoreenActor_LIBS ${vtkVoreenActor_LIBS} 
    ${VRN_DIR}/ext/tiff/lib/libtiff.lib
  )
ENDIF(VRN_WITH_TIFF)

IF(VRN_WITH_FONTRENDERING)
  SET(vtkVoreenActor_LIBS ${vtkVoreenActor_LIBS} 
    ${VRN_DIR}/ext/freetype/lib/freetype.lib
    ${VRN_DIR}/ext/ftgl/lib/ftgl.lib
  )
  SET(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} 
    /FORCE:MULTIPLE
  ) 
ENDIF(VRN_WITH_FONTRENDERING)

IF(VRN_WITH_FFMPEG)
  SET(vtkVoreenActor_LIBS ${vtkVoreenActor_LIBS} 
    ${VRN_DIR}/ext/ffmpeg/win32/avcodec.lib
    ${VRN_DIR}/ext/ffmpeg/win32/avdevice.lib
    ${VRN_DIR}/ext/ffmpeg/win32/avformat.lib
    ${VRN_DIR}/ext/ffmpeg/win32/avutil.lib
    ${VRN_DIR}/ext/ffmpeg/win32/swscale.lib
  )
ENDIF(VRN_WITH_FFMPEG)

IF(VRN_WITH_PYTHON)
  SET(vtkVoreenActor_LIBS ${vtkVoreenActor_LIBS} 
    #  ../../ext/python/lib/$(ConfigurationName)/VS2005/python26
    ${VRN_DIR}/ext/python/lib/$(ConfigurationName)/python26.lib
  )
ENDIF(VRN_WITH_PYTHON)

IF(VRN_WITH_DCMTK)
  SET(vtkVoreenActor_LIBS ${vtkVoreenActor_LIBS} 
    ${VRN_DIR}/src/modules/dicom/ext/dcmtk/win32/$(ConfigurationName)/dcmimage.lib
    ${VRN_DIR}/src/modules/dicom/ext/dcmtk/win32/$(ConfigurationName)/dcmimgle.lib
    ${VRN_DIR}/src/modules/dicom/ext/dcmtk/win32/$(ConfigurationName)/dcmnet.lib
    ${VRN_DIR}/src/modules/dicom/ext/dcmtk/win32/$(ConfigurationName)/dcmdata.lib
    ${VRN_DIR}/src/modules/dicom/ext/dcmtk/win32/$(ConfigurationName)/dcmjpeg.lib
    ${VRN_DIR}/src/modules/dicom/ext/dcmtk/win32/$(ConfigurationName)/ijg8.lib
    ${VRN_DIR}/src/modules/dicom/ext/dcmtk/win32/$(ConfigurationName)/ijg12.lib
    ${VRN_DIR}/src/modules/dicom/ext/dcmtk/win32/$(ConfigurationName)/ijg16.lib
    ${VRN_DIR}/src/modules/dicom/ext/dcmtk/win32/$(ConfigurationName)/ofstd.lib
    ${VRN_DIR}/src/modules/dicom/ext/dcmtk/win32/$(ConfigurationName)/dcmtls.lib
  )
ENDIF(VRN_WITH_DCMTK)

SET( CMAKE_CXX_FLAGS "-D TGT_WITHOUT_DEFINES -D VRN_DIR=\"${VRN_DIR}\"")

#-----------------------------------------------------------------------------
# Most users should not need to change anything below this line.

# Need to include class headers and the configuration header.
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR})

# Voreen headers
INCLUDE_DIRECTORIES(${VRN_DIR}/include/)
INCLUDE_DIRECTORIES(${VRN_DIR}/ext/)
INCLUDE_DIRECTORIES(${VRN_DIR}/ext/glew/include/)

# Setup output paths.
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH
    "Single output directory for building all libraries.")
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH
    "Single output directory for building all executables.")
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

# Find and load VTK settings.
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

# Give user option of building shared or static libraries.  Default to
# the choice made for VTK.
OPTION(BUILD_SHARED_LIBS "Build with shared libraries."
       ${VTK_BUILD_SHARED_LIBS})

# Set vtkVoreenActor_BUILD_SHARED_LIBS to 1 or 0 for use in creating the
# configuration header.
IF(BUILD_SHARED_LIBS)
  SET(vtkVoreenActor_BUILD_SHARED_LIBS_CONFIG 1)
ELSE(BUILD_SHARED_LIBS)
  SET(vtkVoreenActor_BUILD_SHARED_LIBS_CONFIG 0)
ENDIF(BUILD_SHARED_LIBS)

IF(VTK_FOUND)
  # If this is a build tree, provide an option for putting this
  # project's executables and libraries in with VTK's.
  IF (EXISTS ${VTK_DIR}/bin)
    OPTION(USE_VTK_OUTPUT_PATHS
     "Use VTK's output directory for this project's executables and libraries."
     OFF)
    MARK_AS_ADVANCED (USE_VTK_OUTPUT_PATHS)
    IF (USE_VTK_OUTPUT_PATHS)
      SET (LIBRARY_OUTPUT_PATH ${VTK_DIR}/bin)
      SET (EXECUTABLE_OUTPUT_PATH ${VTK_DIR}/bin)
    ENDIF (USE_VTK_OUTPUT_PATHS)
  ENDIF (EXISTS ${VTK_DIR}/bin)

  # Create the header to configure the classes.
  CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/${PROJECT_NAME}Configure.h.in
                 ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Configure.h
                 @ONLY IMMEDIATE)

  # Use wrapping hints for this project.
  SET(VTK_WRAP_HINTS "${PROJECT_SOURCE_DIR}/hints")

  # Create the instantiator for these classes.
  SET(VTK_USE_INSTANTIATOR_NEW 1)
  VTK_MAKE_INSTANTIATOR3("${PROJECT_NAME}Instantiator"
    vtkVoreenActorInstantiator_SRCS
    "${vtkVoreenActor_SRCS}"
    "VTK_${PROJECT_NAME}_EXPORT"
    ${PROJECT_BINARY_DIR}
    "${PROJECT_NAME}Configure.h")

  # Create the library.
  ADD_LIBRARY(${PROJECT_NAME} ${vtkVoreenActor_SRCS} ${vtkVoreenActorInstantiator_SRCS})
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${vtkVoreenActor_LIBS})

  # Create Tcl wrappers if VTK is wrapped in Tcl.
  IF(VTK_WRAP_TCL)
    INCLUDE(${VTK_CMAKE_DIR}/vtkWrapTcl.cmake)
    VTK_WRAP_TCL3(${PROJECT_NAME}TCL vtkVoreenActorTCL_SRCS "${vtkVoreenActor_SRCS}" "")
    ADD_LIBRARY(${PROJECT_NAME}TCL ${vtkVoreenActorTCL_SRCS})
    TARGET_LINK_LIBRARIES(${PROJECT_NAME}TCL ${PROJECT_NAME})
    FOREACH(c ${vtkVoreenActor_LIBS})
      TARGET_LINK_LIBRARIES(${PROJECT_NAME}TCL ${c}TCL)
    ENDFOREACH(c)

    # Create the list of kits that are required.
    SET(vtkVoreenActor_KITS "")
    FOREACH(c ${vtkVoreenActor_LIBS})
      SET(vtkVoreenActor_KITS "${vtkVoreenActor_KITS} ${c}")
    ENDFOREACH(c)
    FOREACH(c ${vtkVoreenActor_SRCS})
      GET_FILENAME_COMPONENT(vtkVoreenActor_LAST_CLASS ${c} NAME_WE)
    ENDFOREACH(c)

    IF(CMAKE_CONFIGURATION_TYPES)
      FOREACH(config ${CMAKE_CONFIGURATION_TYPES})
        SET(vtkVoreenActor_TCL_LIBRARY_DIR "${LIBRARY_OUTPUT_PATH}/${config}")
        CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/pkgIndex.tcl.in
                       ${PROJECT_BINARY_DIR}/${config}/pkgIndex.tcl
                       @ONLY IMMEDIATE)
      ENDFOREACH(config)
    ELSE(CMAKE_CONFIGURATION_TYPES)
      SET(vtkVoreenActor_TCL_LIBRARY_DIR "${LIBRARY_OUTPUT_PATH}")
      CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/pkgIndex.tcl.in
                     ${PROJECT_BINARY_DIR}/pkgIndex.tcl
                     @ONLY IMMEDIATE)
    ENDIF(CMAKE_CONFIGURATION_TYPES)

  ENDIF(VTK_WRAP_TCL)

  # Create Python wrappers if VTK is wrapped in Python.
  IF(VTK_WRAP_PYTHON)
    INCLUDE(${VTK_CMAKE_DIR}/vtkWrapPython.cmake)
    VTK_WRAP_PYTHON3(${PROJECT_NAME}Python vtkVoreenActorPYTHON_SRCS "${vtkVoreenActor_SRCS}")
    ADD_LIBRARY(${PROJECT_NAME}PythonD ${vtkVoreenActorPYTHON_SRCS})
    ADD_LIBRARY(${PROJECT_NAME}Python MODULE ${PROJECT_NAME}PythonInit.cxx)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME}PythonD ${PROJECT_NAME})
    FOREACH(c ${vtkVoreenActor_LIBS})
      TARGET_LINK_LIBRARIES(${PROJECT_NAME}PythonD ${c}PythonD)
    ENDFOREACH(c)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME}Python ${PROJECT_NAME}PythonD)
  ENDIF(VTK_WRAP_PYTHON)

  # Create Java wrappers if VTK is wrapped in Java.
  IF(VTK_WRAP_JAVA)
    INCLUDE(${VTK_CMAKE_DIR}/vtkWrapJava.cmake)
    SET(VTK_WRAP_JAVA3_INIT_DIR "${PROJECT_SOURCE_DIR}")
    SET(VTK_JAVA_HOME ${PROJECT_BINARY_DIR}/java/${PROJECT_NAME})
    MAKE_DIRECTORY(${VTK_JAVA_HOME})
    VTK_WRAP_JAVA3(${PROJECT_NAME}Java vtkVoreenActorJAVA_SRCS "${vtkVoreenActor_SRCS}")
    ADD_LIBRARY(${PROJECT_NAME}Java ${vtkVoreenActorJAVA_SRCS})
    TARGET_LINK_LIBRARIES(${PROJECT_NAME}Java ${PROJECT_NAME})
    FOREACH(c ${vtkVoreenActor_LIBS})
      TARGET_LINK_LIBRARIES(${PROJECT_NAME}Java ${c}Java)
    ENDFOREACH(c)
  ENDIF(VTK_WRAP_JAVA)

  ADD_EXECUTABLE(TestProject main.cpp)
  TARGET_LINK_LIBRARIES(TestProject vtkVoreenActor)

ENDIF(VTK_FOUND)
