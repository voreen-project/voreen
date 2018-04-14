
################################################################################
# Core module resources 
################################################################################
SET(MOD_CORE_MODULECLASS RandomWalkerModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/randomwalker.cpp
    ${MOD_DIR}/processors/randomwalkeranalyzer.cpp
    ${MOD_DIR}/processors/rwmultilabelloopinitializer.cpp
    ${MOD_DIR}/processors/rwmultilabelloopfinalizer.cpp
    
    ${MOD_DIR}/solver/randomwalkersolver.cpp
    ${MOD_DIR}/solver/randomwalkerseeds.cpp
    ${MOD_DIR}/solver/randomwalkerweights.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/randomwalker.h
    ${MOD_DIR}/processors/randomwalkeranalyzer.h
    ${MOD_DIR}/processors/rwmultilabelloopinitializer.h
    ${MOD_DIR}/processors/rwmultilabelloopfinalizer.h
    
    ${MOD_DIR}/solver/randomwalkersolver.h
    ${MOD_DIR}/solver/randomwalkerseeds.h
    ${MOD_DIR}/solver/randomwalkerweights.h
)
   

################################################################################
# Qt module resources 
################################################################################
SET(MOD_QT_MODULECLASS RandomWalkerModuleQt)

SET(MOD_QT_SOURCES
    ${MOD_DIR}/qt/randomwalkeranalyzerwidget.cpp
    ${MOD_DIR}/qt/randomwalkerprocessorwidgetfactory.cpp 
)  
    
SET(MOD_QT_HEADERS
    ${MOD_DIR}/qt/randomwalkeranalyzerwidget.h 
)

SET(MOD_QT_HEADERS_NONMOC
    ${MOD_DIR}/qt/randomwalkerprocessorwidgetfactory.h
)
