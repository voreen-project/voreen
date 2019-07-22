
################################################################################
# External dependency: ITK library
################################################################################

IF(WIN32)
    # Defines the ITK version on windows. Must be manually updated.
    SET(ITK_VERSION 4.12)

    # Define used libs:
    LIST(APPEND ITK_LIB_NAMES 
        "ITKBiasCorrection" "ITKBioCell" "ITKCommon" "ITKDeprecated" "ITKDICOMParser" "itkdouble-conversion"
        "ITKEXPAT" "ITKFEM" "itkgdcmcharls" "itkgdcmCommon" "itkgdcmDICT" "itkgdcmDSED" "itkgdcmIOD" "itkgdcmjpeg8"
        "itkgdcmjpeg12" "itkgdcmjpeg16" "itkgdcmMEXD" "itkgdcmMSFF" "itkgdcmopenjpeg" "itkgdcmsocketxx"
        "ITKgiftiio" "ITKIOBioRad" "ITKIOBMP" "ITKIOCSV" "ITKIOGDCM" "ITKIOGE"
        "ITKIOGIPL" "ITKIOHDF5" "ITKIOImageBase" "ITKIOIPL" "ITKIOJPEG" "ITKIOLSM" "ITKIOMesh" "ITKIOMeta"
        "ITKIOMRC" "ITKIONIFTI" "ITKIONRRD" "ITKIOPNG" "ITKIOSiemens" "ITKIOSpatialObjects" "ITKIOStimulate"
        "ITKIOTIFF" "ITKIOTransformBase" "ITKIOTransformHDF5" "ITKIOTransformInsightLegacy" "ITKIOTransformMatlab"
        "ITKIOVTK" "ITKIOXML" "itkjpeg" "ITKKLMRegionGrowing" "ITKLabelMap" "ITKMesh" "ITKMetaIO" "itknetlib"
        "itkNetlibSlatec" "ITKniftiio" "ITKNrrdIO" "ITKOptimizers" "ITKOptimizersv4" "ITKPath" "itkpng"
        "ITKPolynomials" "ITKQuadEdgeMesh" "ITKReview" "ITKSpatialObjects" "ITKStatistics" "itksys" "itktestlib" 
        "itktiff" "ITKTransform" "itkv3p_netlib" "itkvcl" "ITKVideoCore" "ITKVideoIO" "itkvnl_algo" "itkvnl"
        "ITKVNLInstantiation" "ITKVTK" "ITKWatersheds" "itkzlib" "ITKznz" 
    )
    LIST(APPEND ITK_LIB_NAMES_VERSIONLESS
        "libitkhdf5_cpp" "libitkhdf5"
    )
    
    # Define used dlls:
    LIST(APPEND ITK_DLL_NAMES 
        "ITKCommon" "ITKIOBioRad" "ITKIOBMP" "ITKIOCSV" "ITKIOGDCM" "ITKIOGE" "ITKStatistics" "ITKWatersheds"
        "ITKIOGIPL" "ITKIOHDF5" "ITKIOImageBase" "ITKIOIPL" "ITKIOJPEG" "ITKIOLSM" "ITKIOMesh" "ITKIOMeta"
        "ITKIOMRC" "ITKIONIFTI" "ITKIONRRD" "ITKIOPNG" "ITKIOSiemens" "ITKIOStimulate" "ITKOptimizers"
        "ITKIOTIFF" "ITKIOTransformBase" "ITKIOTransformHDF5" "ITKIOTransformInsightLegacy" "ITKIOTransformMatlab"
        "ITKIOVTK" "ITKIOXML" "ITKTransform"
    )
    
    # check, if itk is present
    IF(NOT EXISTS ${MOD_DIR}/ext/InsightToolkit-${ITK_VERSION}/include/itkConfigure.h)
        MESSAGE(FATAL_ERROR "ITK headers not found (modules/itk/ext/InsightToolkit-${ITK_VERSION}/include/itkConfigure.h). "
            "Copy ITK ${ITK_VERSION} library to modules/itk/ext/ (see http://voreen.uni-muenster.de)")
    ENDIF()

    # set debug and release libraries
    IF(VRN_MSVC2015)            
        SET(ITK_LIB_DIR "${MOD_DIR}/ext/InsightToolkit-${ITK_VERSION}/lib/msvc2015")
    ELSEIF(VRN_MSVC2017)            
        SET(ITK_LIB_DIR "${MOD_DIR}/ext/InsightToolkit-${ITK_VERSION}/lib/msvc2017")
    ENDIF()

    # set includes
    SET(MOD_INCLUDE_DIRECTORIES 
        ${MOD_DIR}/ext/InsightToolkit-${ITK_VERSION}/include
    )

    # set libraries
    FOREACH(elem ${ITK_LIB_NAMES})
        LIST(APPEND MOD_DEBUG_LIBRARIES    "${ITK_LIB_DIR}/debug/${elem}-${ITK_VERSION}.lib")
        LIST(APPEND MOD_RELEASE_LIBRARIES  "${ITK_LIB_DIR}/release/${elem}-${ITK_VERSION}.lib")
    ENDFOREACH()
    
    # set libraries without version number
    FOREACH(elem ${ITK_LIB_NAMES_VERSIONLESS})
        LIST(APPEND MOD_DEBUG_LIBRARIES    "${ITK_LIB_DIR}/debug/${elem}_D.lib")
        LIST(APPEND MOD_RELEASE_LIBRARIES  "${ITK_LIB_DIR}/release/${elem}.lib")
    ENDFOREACH()

    # set dlls
    FOREACH(elem ${ITK_DLL_NAMES})
        LIST(APPEND MOD_DEBUG_DLLS         "${ITK_LIB_DIR}/debug/${elem}-${ITK_VERSION}.dll")
        LIST(APPEND MOD_RELEASE_DLLS       "${ITK_LIB_DIR}/release/${elem}-${ITK_VERSION}.dll")
    ENDFOREACH()

ELSEIF(UNIX)
    FIND_PACKAGE(ITK REQUIRED)
    INCLUDE(${ITK_USE_FILE})
    SET(MOD_LIBRARIES ${ITK_LIBRARIES})
ENDIF()


################################################################################
# Core module resources 
################################################################################
SET(MOD_CORE_MODULECLASS ITKModule)
SET(MOD_REQUIRE_OPENGL_COMPATIBILITY_PROFILE ON)

SET(MOD_CORE_SOURCES 
    ${MOD_DIR}/utils/itkwrapper.cpp
    #${MOD_DIR}/io/itkvolumereader.cpp
    ${MOD_DIR}/processors/anisotropicdiffusion.cpp
    ${MOD_DIR}/processors/doublethreshold.cpp
    ${MOD_DIR}/processors/gradientvectorflow.cpp
    ${MOD_DIR}/processors/itkprocessor.cpp
    ${MOD_DIR}/processors/volumefilter_itk.cpp
    ${MOD_DIR}/processors/valuedregionalmaximaimagefilter.cpp
    ${MOD_DIR}/processors/vesselness.cpp
    ${MOD_DIR}/processors/mutualinformationregistration.cpp
    ${MOD_DIR}/processors/fastmarchingimagefilter.cpp
    ${MOD_DIR}/processors/curveslevelsetimagefilterworkflow.cpp
    ${MOD_DIR}/processors/geodesicactivecontourlevelsetimagefilterworkflow.cpp
    ${MOD_DIR}/processors/geodesicactivecontourshapepriorlevelsetimagefilterworkflow.cpp
    ${MOD_DIR}/processors/narrowbandcurveslevelsetimagefilterworkflow.cpp
    ${MOD_DIR}/processors/narrowbandthresholdsegmentationlevelsetimagefilterworkflow.cpp
    ${MOD_DIR}/processors/shapedetectionlevelsetimagefilterworkflow.cpp
    ${MOD_DIR}/processors/watershedimagefilter.cpp
    ${MOD_DIR}/processors/thresholdlabelerimagefilter.cpp
    ${MOD_DIR}/processors/labelmapoverlayimagefilter.cpp
    ${MOD_DIR}/processors/labelmapcontouroverlayimagefilter.cpp
    ${MOD_DIR}/processors/bayesianclassifierimagefilter.cpp
)

SET(MOD_CORE_HEADERS 
    ${MOD_DIR}/utils/itkwrapper.h
    #${MOD_DIR}/io/itkvolumereader.h
    ${MOD_DIR}/processors/anisotropicdiffusion.h
    ${MOD_DIR}/processors/doublethreshold.h
    ${MOD_DIR}/processors/gradientvectorflow.h
    ${MOD_DIR}/processors/itkprocessor.h
    ${MOD_DIR}/processors/volumefilter_itk.h
    ${MOD_DIR}/processors/valuedregionalmaximaimagefilter.h
    ${MOD_DIR}/processors/vesselness.h
    ${MOD_DIR}/processors/mutualinformationregistration.h
    ${MOD_DIR}/processors/fastmarchingimagefilter.h
    ${MOD_DIR}/processors/curveslevelsetimagefilterworkflow.h
    ${MOD_DIR}/processors/geodesicactivecontourlevelsetimagefilterworkflow.h
    ${MOD_DIR}/processors/geodesicactivecontourshapepriorlevelsetimagefilterworkflow.h
    ${MOD_DIR}/processors/narrowbandcurveslevelsetimagefilterworkflow.h
    ${MOD_DIR}/processors/narrowbandthresholdsegmentationlevelsetimagefilterworkflow.h
    ${MOD_DIR}/processors/shapedetectionlevelsetimagefilterworkflow.h
    ${MOD_DIR}/processors/watershedimagefilter.h
    ${MOD_DIR}/processors/thresholdlabelerimagefilter.h
    ${MOD_DIR}/processors/labelmapoverlayimagefilter.h
    ${MOD_DIR}/processors/labelmapcontouroverlayimagefilter.h
    ${MOD_DIR}/processors/bayesianclassifierimagefilter.h
)
