################################################################################
# Core module resources 
################################################################################
SET(MOD_CORE_MODULECLASS TwoManifoldGeometryClippingModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/twomanifoldgeometryclipping.cpp
    ${MOD_DIR}/util/twomanifoldclippingutil.cpp
    ${MOD_DIR}/util/polygoncreator.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/twomanifoldgeometryclipping.h
    ${MOD_DIR}/util/twomanifoldclippingutil.h
    ${MOD_DIR}/util/loopconstructor.h
    ${MOD_DIR}/util/polygoncreator.h
    ${MOD_DIR}/util/polygontriangulator.h
)

# deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
)
   
