################################################################################
# Core module resources
################################################################################
SET(MOD_CORE_MODULECLASS FlowreenModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/flowreenmodule.cpp

    # datastructures
    ${MOD_DIR}/datastructures/streamline.cpp
    ${MOD_DIR}/datastructures/streamlinebundle.cpp
    ${MOD_DIR}/datastructures/streamlinelist.cpp
    ${MOD_DIR}/datastructures/streamlinelistbase.cpp
    ${MOD_DIR}/datastructures/streamlinelistdecorator.cpp
    ${MOD_DIR}/datastructures/streamlinelistobserver.cpp

    # ports
    ${MOD_DIR}/ports/streamlinelistport.cpp

    # processors
    ${MOD_DIR}/processors/flowdirectionoverlay.cpp
    ${MOD_DIR}/processors/streamlinerenderer3d.cpp
    ${MOD_DIR}/processors/streamline/streamlinecombine.cpp
    ${MOD_DIR}/processors/streamline/streamlinecreator.cpp
    ${MOD_DIR}/processors/streamline/streamlinerotation.cpp
    ${MOD_DIR}/processors/streamline/streamlinesave.cpp
    ${MOD_DIR}/processors/streamline/streamlineselector.cpp
    ${MOD_DIR}/processors/streamline/streamlinesource.cpp
    ${MOD_DIR}/processors/streamline/streamlinetoboundingbox.cpp
    ${MOD_DIR}/processors/streamline/pathlinecreator.cpp

    # utils
    ${MOD_DIR}/utils/streamlinebundledetectorbackgroundthread.cpp
    ${MOD_DIR}/utils/streamlinecreatorbackgroundthread.cpp
    ${MOD_DIR}/utils/pathlinecreatorbackgroundthread.cpp
)

IF(VRN_OPENGL_COMPATIBILITY_PROFILE)
    #check if core profile compatible or port
    LIST(APPEND MOD_CORE_SOURCES
        # datastructures
        ${MOD_DIR}/datastructures/deprecated/flow2d.cpp
        ${MOD_DIR}/datastructures/deprecated/flow3d.cpp
        ${MOD_DIR}/datastructures/deprecated/simpletexture.cpp
        ${MOD_DIR}/datastructures/deprecated/streamlinetexture.cpp
        ${MOD_DIR}/datastructures/deprecated/volumeflow3d.cpp
        ${MOD_DIR}/datastructures/deprecated/volumeoperatorflowmagnitude.cpp

        # io
        ${MOD_DIR}/io/flowreader.cpp

        # processors
        ${MOD_DIR}/processors/flowarrowrenderer2d.cpp
        ${MOD_DIR}/processors/flowarrowrenderer3d.cpp
        ${MOD_DIR}/processors/flowmagnitudes3d.cpp
        ${MOD_DIR}/processors/floworthogonalslicerenderer.cpp
        ${MOD_DIR}/processors/flowreenadapter.cpp
        ${MOD_DIR}/processors/flowreenprocessor.cpp
        ${MOD_DIR}/processors/flowslicerenderer.cpp
        ${MOD_DIR}/processors/flowslicerenderer2d.cpp
        ${MOD_DIR}/processors/flowslicerenderer3d.cpp
        ${MOD_DIR}/processors/flowstreamlinestexture3d.cpp
        ${MOD_DIR}/processors/pathlinerenderer3d.cpp
        ${MOD_DIR}/processors/pathlinerenderer3d.cpp

        # utils
        ${MOD_DIR}/utils/colorcodingability.cpp
        ${MOD_DIR}/utils/flowmath.cpp
    )
ENDIF()
    
SET(MOD_CORE_HEADERS
    ${MOD_DIR}/flowreenmodule.h

    # datastructures
    ${MOD_DIR}/datastructures/streamline.h
    ${MOD_DIR}/datastructures/streamlinebundle.h
    ${MOD_DIR}/datastructures/streamlinelist.h
    ${MOD_DIR}/datastructures/streamlinelistbase.h
    ${MOD_DIR}/datastructures/streamlinelistdecorator.h
    ${MOD_DIR}/datastructures/streamlinelistobserver.h

    # ports
    ${MOD_DIR}/ports/streamlinelistport.h

    # processors
    ${MOD_DIR}/processors/flowdirectionoverlay.h
    ${MOD_DIR}/processors/streamlinerenderer3d.h
    ${MOD_DIR}/processors/streamline/streamlinecombine.h
    ${MOD_DIR}/processors/streamline/streamlinecreator.h
    ${MOD_DIR}/processors/streamline/streamlinerotation.h
    ${MOD_DIR}/processors/streamline/streamlinesave.h
    ${MOD_DIR}/processors/streamline/streamlineselector.h
    ${MOD_DIR}/processors/streamline/streamlinesource.h
    ${MOD_DIR}/processors/streamline/streamlinetoboundingbox.h
    ${MOD_DIR}/processors/streamline/pathlinecreator.h

    # utils
    ${MOD_DIR}/utils/streamlinebundledetectorbackgroundthread.h
    ${MOD_DIR}/utils/streamlinecreatorbackgroundthread.h
    ${MOD_DIR}/utils/pathlinecreatorbackgroundthread.h
)
 
IF(VRN_OPENGL_COMPATIBILITY_PROFILE)
    #check if core profile compatible or port
    LIST(APPEND MOD_CORE_HEADERS
        # # datastructures
        ${MOD_DIR}/datastructures/deprecated/flow2d.h
        ${MOD_DIR}/datastructures/deprecated/flow3d.h
        ${MOD_DIR}/datastructures/deprecated/simpletexture.h
        ${MOD_DIR}/datastructures/deprecated/streamlinetexture.h
        ${MOD_DIR}/datastructures/deprecated/volumeflow3d.h
        ${MOD_DIR}/datastructures/deprecated/volumeoperatorflowmagnitude.h
        ${MOD_DIR}/datastructures/deprecated/volumeoperatorintensitymask.h

        # io
        ${MOD_DIR}/io/flowreader.h

        # processors
        ${MOD_DIR}/processors/flowarrowrenderer2d.h
        ${MOD_DIR}/processors/flowarrowrenderer3d.h
        ${MOD_DIR}/processors/flowmagnitudes3d.h
        ${MOD_DIR}/processors/floworthogonalslicerenderer.h
        ${MOD_DIR}/processors/flowreenadapter.h
        ${MOD_DIR}/processors/flowreenprocessor.h
        ${MOD_DIR}/processors/flowslicerenderer.h
        ${MOD_DIR}/processors/flowslicerenderer2d.h
        ${MOD_DIR}/processors/flowslicerenderer3d.h
        ${MOD_DIR}/processors/flowstreamlinestexture3d.h
        ${MOD_DIR}/processors/pathlinerenderer3d.h

        # utils
        ${MOD_DIR}/utils/colorcodingability.h
        ${MOD_DIR}/utils/flowmath.h
    )
ENDIF()

# Deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
    ${MOD_DIR}/data
    ${MOD_DIR}/transferfuncs
)
