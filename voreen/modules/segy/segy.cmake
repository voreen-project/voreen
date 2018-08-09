SET(MOD_CORE_MODULECLASS SEGYModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/io/segyvolumereader.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/io/segyvolumereader.h
)

# deployment
SET(MOD_INSTALL_FILES
    ${MOD_DIR}/README.txt
)
 
