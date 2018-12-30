include(FindPackageHandleStandardArgs)

set(STIM_INCLUDE_DIR $ENV{STIMLIB_PATH})

find_package_handle_standard_args(STIM DEFAULT_MSG STIM_INCLUDE_DIR)

if(STIM_FOUND)
    set(STIM_INCLUDE_DIRS ${STIM_INCLUDE_DIR})
endif()