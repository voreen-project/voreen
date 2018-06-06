# module class must reside in sample/samplemodule.h + sample/samplemodule.cpp
SET(MOD_CORE_MODULECLASS SampleModule)

# module's core source files, path relative to module dir
SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/sampleprocessor.cpp
    ${MOD_DIR}/processors/samplerenderprocessor.cpp
)

# module's core header files, path relative to module dir
SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/sampleprocessor.h
    ${MOD_DIR}/processors/samplerenderprocessor.h
)