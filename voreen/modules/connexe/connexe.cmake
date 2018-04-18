SET(MOD_CORE_MODULECLASS ConnexeModule)

# external dependency: connexe library
SET(MOD_INCLUDE_DIRECTORIES 
    ${MOD_DIR}/ext/connexe
)
SET(MOD_INSTALL_FILES
    ${MOD_DIR}/ext/connexe/COPYRIGHT
) 

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/connectedcomponents2d.cpp
    ${MOD_DIR}/processors/connectedcomponents3d.cpp
# ext
    ${MOD_DIR}/ext/connexe/connexe.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/connectedcomponents2d.h
    ${MOD_DIR}/processors/connectedcomponents3d.h
)
