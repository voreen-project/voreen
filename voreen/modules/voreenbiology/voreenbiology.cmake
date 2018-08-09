################################################################################
# core module resources
################################################################################
SET(MOD_CORE_MODULECLASS VoreenBiologyModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/startuplinkingprocessor.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/startuplinkingprocessor.h
)

################################################################################
# dongle action
################################################################################
IF(WIN32)
    SET(MOD_INCLUDE_DIRECTORIES "${MOD_DIR}/ext/include")

    LIST(APPEND MOD_CORE_SOURCES ${MOD_DIR}/ext/src/DongleInterface.cpp)
    
    SET(MOD_LIBRARIES 
        "${MOD_DIR}/ext/lib/hasp_windows_x64.lib"
    )

    SET(MOD_DEBUG_DLLS
        "${MOD_DIR}/ext/lib/hasp_windows_x64.dll"
    )
    SET(MOD_RELEASE_DLLS ${MOD_DEBUG_DLLS})
ENDIF()

