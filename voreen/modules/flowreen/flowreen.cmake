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
    ${MOD_DIR}/processors/streamline/streamlinetogeometry.cpp

    # utils
    ${MOD_DIR}/utils/streamlinebundledetectorbackgroundthread.cpp
    ${MOD_DIR}/utils/streamlinecreatorbackgroundthread.cpp
)
    
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
    ${MOD_DIR}/processors/streamline/streamlinetogeometry.h

    # utils
    ${MOD_DIR}/utils/streamlinebundledetectorbackgroundthread.h
    ${MOD_DIR}/utils/streamlinecreatorbackgroundthread.h
)

# Deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
)
