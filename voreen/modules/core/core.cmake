################################################################################
# Core module resources 
################################################################################
SET(MOD_CORE_MODULECLASS CoreModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/input/geometrysource.cpp
    ${MOD_DIR}/processors/input/imageselector.cpp
    ${MOD_DIR}/processors/input/imagesequencesource.cpp
    ${MOD_DIR}/processors/input/imagesource.cpp
    ${MOD_DIR}/processors/input/octreecreator.cpp
    ${MOD_DIR}/processors/input/textsource.cpp
    ${MOD_DIR}/processors/input/volumelistsource.cpp
    ${MOD_DIR}/processors/input/volumeselector.cpp
    ${MOD_DIR}/processors/input/volumesource.cpp
    
    ${MOD_DIR}/processors/output/canvasrenderer.cpp
    ${MOD_DIR}/processors/output/geometrysave.cpp
    ${MOD_DIR}/processors/output/imagesequencesave.cpp
    ${MOD_DIR}/processors/output/textsave.cpp
    ${MOD_DIR}/processors/output/volumelistsave.cpp
    ${MOD_DIR}/processors/output/volumesave.cpp
    ${MOD_DIR}/processors/output/octreesave.cpp

    ${MOD_DIR}/io/datvolumereader.cpp
    ${MOD_DIR}/io/datvolumewriter.cpp
    ${MOD_DIR}/io/rawvolumereader.cpp
    ${MOD_DIR}/io/vvdformat.cpp
    ${MOD_DIR}/io/vvodformat.cpp
    ${MOD_DIR}/io/vvdvolumereader.cpp
    ${MOD_DIR}/io/vvdvolumewriter.cpp
    ${MOD_DIR}/io/vvodvolumereader.cpp
    ${MOD_DIR}/io/vvodvolumewriter.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/input/geometrysource.h
    ${MOD_DIR}/processors/input/imageselector.h
    ${MOD_DIR}/processors/input/imagesequencesource.h
    ${MOD_DIR}/processors/input/imagesource.h
    ${MOD_DIR}/processors/input/octreecreator.h
    ${MOD_DIR}/processors/input/textsource.h
    ${MOD_DIR}/processors/input/volumelistsource.h
    ${MOD_DIR}/processors/input/volumeselector.h
    ${MOD_DIR}/processors/input/volumesource.h
    
    ${MOD_DIR}/processors/output/canvasrenderer.h
    ${MOD_DIR}/processors/output/geometrysave.h
    ${MOD_DIR}/processors/output/imagesequencesave.h
    ${MOD_DIR}/processors/output/textsave.h
    ${MOD_DIR}/processors/output/volumelistsave.h
    ${MOD_DIR}/processors/output/volumesave.h
    ${MOD_DIR}/processors/output/octreesave.h

    ${MOD_DIR}/io/datvolumereader.h
    ${MOD_DIR}/io/datvolumewriter.h
    ${MOD_DIR}/io/rawvolumereader.h
    ${MOD_DIR}/io/vvdformat.h
    ${MOD_DIR}/io/vvodformat.h
    ${MOD_DIR}/io/vvdvolumereader.h
    ${MOD_DIR}/io/vvdvolumewriter.h
    ${MOD_DIR}/io/vvodvolumereader.h
    ${MOD_DIR}/io/vvodvolumewriter.h
)
   

################################################################################
# Qt module resources 
################################################################################
SET(MOD_QT_MODULECLASS CoreModuleQt)

SET(MOD_QT_SOURCES
    ${MOD_DIR}/qt/processor/coreprocessorwidgetfactory.cpp
    ${MOD_DIR}/qt/processor/canvasrendererwidget.cpp
)  
    
SET(MOD_QT_HEADERS
    ${MOD_DIR}/qt/processor/canvasrendererwidget.h
)

SET(MOD_QT_HEADERS_NONMOC
    ${MOD_DIR}/qt/processor/coreprocessorwidgetfactory.h
)


################################################################################
# VoreenVE module resources (currently not present)
################################################################################
# SET(MOD_VE_MODULECLASS CoreModuleVE)

# SET(MOD_VE_SOURCES
# )  
    
# SET(MOD_VE_HEADERS
# )

# SET(MOD_VE_HEADERS_NONMOC
# )


################################################################################
# Deployment
################################################################################
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
)
