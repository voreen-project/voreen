################################################################################
# Core module resources
################################################################################

# Dependencies
IF(NOT VRN_MODULE_ENSEMBLEANALYSIS) # for plotting
    MESSAGE(FATAL_ERROR "EnsembleAnalysisExtra Module requires EnsembleAnalysis Module")
ENDIF()

SET(MOD_CORE_MODULECLASS EnsembleAnalysisExtraModule)

# deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
    ${MOD_DIR}/workspaces
)

SET(MOD_CORE_SOURCES
    #Datastructures
    ${MOD_DIR}/datastructures/fieldplotdata.cpp

    #IO
    ${MOD_DIR}/io/fieldplotsave.cpp
    ${MOD_DIR}/io/fieldplotsource.cpp

    #Processors
    ${MOD_DIR}/processors/fieldparallelplotcreator.cpp
    ${MOD_DIR}/processors/fieldparallelplothistogram.cpp
    ${MOD_DIR}/processors/fieldparallelplotviewer.cpp
    ${MOD_DIR}/processors/physicalclippinglinker.cpp

    #Ports
    ${MOD_DIR}/ports/fieldplotdataport.cpp
)

SET(MOD_CORE_HEADERS
    #Datastructures
    ${MOD_DIR}/datastructures/fieldplotdata.h

    #IO
    ${MOD_DIR}/io/fieldplotsave.h
    ${MOD_DIR}/io/fieldplotsource.h

    #Processors
    ${MOD_DIR}/processors/fieldparallelplotcreator.h
    ${MOD_DIR}/processors/fieldparallelplothistogram.h
    ${MOD_DIR}/processors/fieldparallelplotviewer.h
    ${MOD_DIR}/processors/physicalclippinglinker.h


    #Ports
    ${MOD_DIR}/ports/fieldplotdataport.h

)

