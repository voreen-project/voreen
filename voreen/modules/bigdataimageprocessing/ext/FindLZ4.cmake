# Try to find liblz4
#
# Once done, this will define
#
# LZ4_FOUND
# LZ4_INCLUDE_DIR
# LZ4_LIBRARIES
# LZ4_VERSION_STRING
# LZ4_VERSION_MAJOR
# LZ4_VERSION_MINOR
# LZ4_VERSION_RELEASE

find_path(LZ4_INCLUDE_DIR NAMES lz4.h)

if(LZ4_INCLUDE_DIR AND EXISTS "${LZ4_INCLUDE_DIR}/lz4.h")
  foreach(ver "MAJOR" "MINOR" "RELEASE")
    file(STRINGS "${LZ4_INCLUDE_DIR}/lz4.h" LZ4_VER_${ver}_LINE
      REGEX "^#define[ \t]+LZ4_VERSION_${ver}[ \t]+[0-9]+[ \t]+.*$")
    string(REGEX REPLACE "^#define[ \t]+LZ4_VERSION_${ver}[ \t]+([0-9]+)[ \t]+.*$"
      "\\1" LZ4_VERSION_${ver} "${LZ4_VER_${ver}_LINE}")
    unset(${LZ4_VER_${ver}_LINE})
  endforeach()
  set(LZ4_VERSION_STRING
    "${LZ4_VERSION_MAJOR}.${LZ4_VERSION_MINOR}.${LZ4_VERSION_RELEASE}")
endif()

find_library(LZ4_LIBRARIES NAMES lz4)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LZ4
  REQUIRED_VARS LZ4_LIBRARIES LZ4_INCLUDE_DIR
  VERSION_VAR LZ4_VERSION_STRING)

mark_as_advanced(LZ4_INCLUDE_DIR LZ4_LIBRARIES)
