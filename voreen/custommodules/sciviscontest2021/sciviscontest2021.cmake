# module class must reside in sample/samplemodule.h + sample/samplemodule.cpp
SET(MOD_CORE_MODULECLASS SciVisContest2021Module)

# module's core source files, path relative to module dir
SET(MOD_CORE_SOURCES
    ${MOD_DIR}/sciviscontest2021module.cpp

    #datastructures
    ${MOD_DIR}/datastructures/sphericalvolumeramproxy.cpp
    ${MOD_DIR}/datastructures/timeserieslist.cpp

    #ports
    ${MOD_DIR}/ports/timeserieslistport.cpp

    #processors
    ${MOD_DIR}/processors/componentidselector.cpp
    ${MOD_DIR}/processors/connectedcomponenttracker.cpp
    ${MOD_DIR}/processors/diskseedpointcreator.cpp
    ${MOD_DIR}/processors/distancetooriginplot.cpp
    ${MOD_DIR}/processors/featureextractor.cpp
    ${MOD_DIR}/processors/sphericalcuttingplaneselector.cpp
    ${MOD_DIR}/processors/sphericalfakevolume.cpp
    ${MOD_DIR}/processors/sphericalgeometrytransformation.cpp
    ${MOD_DIR}/processors/sphericalstreamlinetransformation.cpp
    ${MOD_DIR}/processors/sphericalraycaster.cpp
    #${MOD_DIR}/processors/sphericalsliceviewer.cpp
    ${MOD_DIR}/processors/sphericalvolumelistproxy.cpp
    ${MOD_DIR}/processors/sphericalvolumeproxy.cpp
    ${MOD_DIR}/processors/starcoordinates.cpp
    ${MOD_DIR}/processors/timeserieslistcreator.cpp
    ${MOD_DIR}/processors/timeserieslistsave.cpp
    ${MOD_DIR}/processors/timeserieslistsource.cpp
    ${MOD_DIR}/processors/timeseriesplot.cpp
    ${MOD_DIR}/processors/volumeoverlaptracker.cpp
)

# module's core header files, path relative to module dir
SET(MOD_CORE_HEADERS
    ${MOD_DIR}/sciviscontest2021module.h

    #datastructures
    ${MOD_DIR}/datastructures/sphericalvolumeramproxy.h
    ${MOD_DIR}/datastructures/timeserieslist.h

    #ports
    ${MOD_DIR}/ports/timeserieslistport.h

    #processors
    ${MOD_DIR}/processors/componentidselector.h
    ${MOD_DIR}/processors/connectedcomponenttracker.h
    ${MOD_DIR}/processors/diskseedpointcreator.h
    ${MOD_DIR}/processors/distancetooriginplot.h
    ${MOD_DIR}/processors/featureextractor.h
    ${MOD_DIR}/processors/sphericalcuttingplaneselector.h
    ${MOD_DIR}/processors/sphericalfakevolume.h
    ${MOD_DIR}/processors/sphericalgeometrytransformation.h
    ${MOD_DIR}/processors/sphericalstreamlinetransformation.h
    ${MOD_DIR}/processors/sphericalraycaster.h
    #${MOD_DIR}/processors/sphericalsliceviewer.h
    ${MOD_DIR}/processors/sphericalvolumelistproxy.h
    ${MOD_DIR}/processors/sphericalvolumeproxy.h
    ${MOD_DIR}/processors/starcoordinates.h
    ${MOD_DIR}/processors/timeserieslistcreator.h
    ${MOD_DIR}/processors/timeserieslistsave.h
    ${MOD_DIR}/processors/timeserieslistsource.h
    ${MOD_DIR}/processors/timeseriesplot.h
    ${MOD_DIR}/processors/volumeoverlaptracker.h
)

# deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
    ${MOD_DIR}/workspaces
)