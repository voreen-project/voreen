SET(MOD_CORE_MODULECLASS PVMModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/io/ddsbase.cpp
    ${MOD_DIR}/io/pvmvolumereader.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/io/codebase.h
    ${MOD_DIR}/io/ddsbase.h
    ${MOD_DIR}/io/pvmvolumereader.h
)

# deployment
SET(MOD_INSTALL_FILES
    ${MOD_DIR}/LICENSE.txt
)
 
