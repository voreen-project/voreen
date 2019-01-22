# Try to find Boost library and include path. Once done this will define:
# Boost_FOUND
# Boost_DEFINITIONS
# Boost_INCLUDE_DIRS
# Boost_LIBRARIES (containing both debug and release libraries on win32)
# win32: Boost_DEBUG_LIBRARIES, Boost_RELEASE_LIBRARIES, Boost_DEBUG_DLLS, Boost_RELEASE_DLLS

# Defines the Boost version on windows. Must be manually updated. 
SET(VRN_Boost_VERSION_WIN32 1_65_1)
SET(VRN_Boost_VERSION_UNIX 1.54.0)

# Define used dynamic libs: "container" not included, as it does not seem to be present on Ubuntu linux...
LIST(APPEND Boost_DYNAMIC_NAMES 
   "atomic" "chrono" "context"
   "coroutine" "date_time" "filesystem" "graph" "iostreams"
   "locale" "log_setup" "log" "math_c99f" "math_c99l"
   "math_c99" "math_tr1f" "math_tr1l" "math_tr1" "prg_exec_monitor"
   "program_options" "random" "regex" "serialization"
   "system" "thread" "timer" "unit_test_framework"
   "wave" "wserialization"
)

# define used static libs
LIST(APPEND Boost_STATIC_NAMES 
    "exception"
)

IF (WIN32)

    # explicitly add "container" and "zlib" for windows
    LIST(APPEND Boost_DYNAMIC_NAMES "container" "zlib")

    # used to get the proper library endings 
    SET(VRN_Boost_LIB_GD_STR "")
    SET(VRN_Boost_LIB_STR "")
    SET(Boost_DEFINITIONS "-DBOOST_ALL_NO_LIB")
    SET(Boost_DIR "${VRN_HOME}/ext/boost" CACHE PATH "If boost is not found, set this path")
    SET(Boost_INCLUDE_DIRS "${Boost_DIR}/include")
    
    # set debug and release libraries
	IF(VRN_MSVC2012)            
		SET(Boost_LIB_DIR "${Boost_DIR}/lib/msvc2012")
		SET(VRN_Boost_LIB_STR "-vc110-mt-${VRN_Boost_VERSION_WIN32}")
		SET(VRN_Boost_LIB_GD_STR "-vc110-mt-gd-${VRN_Boost_VERSION_WIN32}")
	ELSEIF(VRN_MSVC2013)            
		SET(Boost_LIB_DIR "${Boost_DIR}/lib/msvc2013")
		SET(VRN_Boost_LIB_STR "-vc120-mt-${VRN_Boost_VERSION_WIN32}")
		SET(VRN_Boost_LIB_GD_STR "-vc120-mt-gd-${VRN_Boost_VERSION_WIN32}")
	ELSEIF(VRN_MSVC2015)            
		SET(Boost_LIB_DIR "${Boost_DIR}/lib/msvc2015")
		SET(VRN_Boost_LIB_STR "-vc140-mt-${VRN_Boost_VERSION_WIN32}")
		SET(VRN_Boost_LIB_GD_STR "-vc140-mt-gd-${VRN_Boost_VERSION_WIN32}")
	ELSEIF(VRN_MSVC2017)            
		SET(Boost_LIB_DIR "${Boost_DIR}/lib/msvc2017")
		SET(VRN_Boost_LIB_STR "-vc141-mt-${VRN_Boost_VERSION_WIN32}")
		SET(VRN_Boost_LIB_GD_STR "-vc141-mt-gd-${VRN_Boost_VERSION_WIN32}")
	ELSE()
		MESSAGE(FATAL_ERROR "Unknown 64Bit Windows compiler!")
	ENDIF()
  
    # construct windows file names for dynamic and static libs 
    FOREACH(elem ${Boost_DYNAMIC_NAMES})
        LIST(APPEND Boost_DEBUG_LIBRARIES   "${Boost_LIB_DIR}/boost_${elem}${VRN_Boost_LIB_GD_STR}.lib")
        LIST(APPEND Boost_DEBUG_DLLS        "${Boost_LIB_DIR}/boost_${elem}${VRN_Boost_LIB_GD_STR}.dll")
        LIST(APPEND Boost_RELEASE_LIBRARIES "${Boost_LIB_DIR}/boost_${elem}${VRN_Boost_LIB_STR}.lib")
        LIST(APPEND Boost_RELEASE_DLLS      "${Boost_LIB_DIR}/boost_${elem}${VRN_Boost_LIB_STR}.dll")
    ENDFOREACH()
    FOREACH(elem ${Boost_STATIC_NAMES})
        LIST(APPEND Boost_DEBUG_LIBRARIES   "${Boost_LIB_DIR}/libboost_${elem}${VRN_Boost_LIB_GD_STR}.lib")
        LIST(APPEND Boost_RELEASE_LIBRARIES "${Boost_LIB_DIR}/libboost_${elem}${VRN_Boost_LIB_STR}.lib") 
    ENDFOREACH()

    FOREACH(lib ${Boost_DEBUG_LIBRARIES})
        LIST(APPEND Boost_LIBRARIES debug ${lib})
    ENDFOREACH()
    FOREACH(lib ${Boost_RELEASE_LIBRARIES})
        LIST(APPEND Boost_LIBRARIES optimized ${lib})
    ENDFOREACH()
   
    IF(Boost_INCLUDE_DIRS AND Boost_LIBRARIES)
        SET(Boost_FOUND TRUE)
    ELSE()
        SET(Boost_FOUND FALSE)
    ENDIF()

ELSE(WIN32)
    FIND_PACKAGE(Boost ${VRN_Boost_VERSION_UNIX} REQUIRED ${Boost_DYNAMIC_NAMES} )
ENDIF(WIN32)

MARK_AS_ADVANCED(Boost_DIR Boost_INCLUDE_DIRS)
