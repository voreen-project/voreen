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

    # ports
    ${MOD_DIR}/ports/streamlinelistport.cpp

    # processors
    # geometry
    ${MOD_DIR}/processors/geometry/streamlinetoboundingbox.cpp
    ${MOD_DIR}/processors/geometry/streamlinetogeometry.cpp

    # render
    ${MOD_DIR}/processors/render/flowdirectionoverlay.cpp
    ${MOD_DIR}/processors/render/streamlinerenderer3d.cpp

    # streamline
    ${MOD_DIR}/processors/streamline/pathlinecreator.cpp
    ${MOD_DIR}/processors/streamline/streamlinebundledetector.cpp
    ${MOD_DIR}/processors/streamline/streamlinecombine.cpp
    ${MOD_DIR}/processors/streamline/streamlinecreator.cpp
    ${MOD_DIR}/processors/streamline/streamlinerotation.cpp
    ${MOD_DIR}/processors/streamline/streamlinesave.cpp
    ${MOD_DIR}/processors/streamline/streamlineselector.cpp
    ${MOD_DIR}/processors/streamline/streamlinesource.cpp

    # utils
    ${MOD_DIR}/utils/flowutils.cpp
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
    # geometry
    ${MOD_DIR}/processors/geometry/streamlinetoboundingbox.h
    ${MOD_DIR}/processors/geometry/streamlinetogeometry.h

    # render
    ${MOD_DIR}/processors/render/flowdirectionoverlay.h
    ${MOD_DIR}/processors/render/streamlinerenderer3d.h

    # streamline
    ${MOD_DIR}/processors/streamline/pathlinecreator.h
    ${MOD_DIR}/processors/streamline/streamlinebundledetector.h
    ${MOD_DIR}/processors/streamline/streamlinecombine.h
    ${MOD_DIR}/processors/streamline/streamlinecreator.h
    ${MOD_DIR}/processors/streamline/streamlinerotation.h
    ${MOD_DIR}/processors/streamline/streamlinesave.h
    ${MOD_DIR}/processors/streamline/streamlineselector.h
    ${MOD_DIR}/processors/streamline/streamlinesource.h

    # utils
    ${MOD_DIR}/utils/flowutils.h
)

# Deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
)
