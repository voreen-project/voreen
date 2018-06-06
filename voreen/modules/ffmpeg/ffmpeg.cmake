
################################################################################
# External dependency: FFmpeg library
################################################################################

IF(WIN32)
    SET(MOD_INCLUDE_DIRECTORIES "${MOD_DIR}/ext/ffmpeg/include")

    SET(MOD_LIBRARIES 
        ${MOD_DIR}/ext/ffmpeg/lib/avcodec.lib
        ${MOD_DIR}/ext/ffmpeg/lib/avformat.lib
        ${MOD_DIR}/ext/ffmpeg/lib/avutil.lib
        ${MOD_DIR}/ext/ffmpeg/lib/swscale.lib
    )

    SET(MOD_DEBUG_DLLS
        ${MOD_DIR}/ext/ffmpeg/lib/avcodec-55.dll
        ${MOD_DIR}/ext/ffmpeg/lib/avformat-55.dll
        ${MOD_DIR}/ext/ffmpeg/lib/avutil-52.dll
        ${MOD_DIR}/ext/ffmpeg/lib/swscale-2.dll
    )

    SET(MOD_RELEASE_DLLS ${MOD_DEBUG_DLLS})

    # deployment
    SET(MOD_INSTALL_FILES
        ${MOD_DIR}/ext/ffmpeg/lgpl-2.1.txt
    )

ELSEIF(UNIX)

    FIND_PATH(FFMPEG_INCLUDE_DIR libavcodec/avcodec.h
        PATHS
        /usr/include
        /usr/local/include
        /opt/local/include
    )

    IF(FFMPEG_INCLUDE_DIR)
        MESSAGE(STATUS "FFmpeg include directory: ${FFMPEG_INCLUDE_DIR}")
    ELSE()
        MESSAGE(FATAL_ERROR " FFmpeg include directory could not be found!")
    ENDIF()

    FIND_LIBRARY( FFMPEG_avcodec_LIBRARY
        NAMES avcodec
        PATHS
        /usr/lib
        /usr/local/lib
        /usr/lib64
        /usr/local/lib64
        /usr/lib/x86_64-linux-gnu
        /opt/local/lib
    )

    IF(FFMPEG_avcodec_LIBRARY)
        MESSAGE(STATUS "FFmpeg libavcodec location: ${FFMPEG_avcodec_LIBRARY}")
    ELSE()
        MESSAGE(FATAL_ERROR " FFmpeg libavcodec library could not be found!")
    ENDIF()

    FIND_LIBRARY( FFMPEG_avformat_LIBRARY
        NAMES avformat
        PATHS
        /usr/lib
        /usr/local/lib
        /usr/lib64
        /usr/local/lib64
        /usr/lib/x86_64-linux-gnu
        /opt/local/lib
    )

    IF(FFMPEG_avformat_LIBRARY)
        MESSAGE(STATUS "FFmpeg libavformat location: ${FFMPEG_avformat_LIBRARY}")
    ELSE()
        MESSAGE(FATAL_ERROR " FFmpeg libavformat library could not be found!")
    ENDIF()

    FIND_LIBRARY( FFMPEG_avutil_LIBRARY
        NAMES avutil
        PATHS
        /usr/lib
        /usr/local/lib
        /usr/lib64
        /usr/local/lib64
        /usr/lib/x86_64-linux-gnu
        /opt/local/lib
    )

    IF(FFMPEG_avutil_LIBRARY)
        MESSAGE(STATUS "FFmpeg libavutil location: ${FFMPEG_avutil_LIBRARY}")
    ELSE()
        MESSAGE(FATAL_ERROR " FFmpeg libavutil library could not be found!")
    ENDIF()

    FIND_LIBRARY( FFMPEG_swscale_LIBRARY
        NAMES swscale
        PATHS
        /usr/lib
        /usr/local/lib
        /usr/lib64
        /usr/local/lib64
        /usr/lib/x86_64-linux-gnu
        /opt/local/lib
    )

    IF(FFMPEG_swscale_LIBRARY)
        MESSAGE(STATUS "FFmpeg libswscale location: ${FFMPEG_swscale_LIBRARY}")
    ELSE()
        MESSAGE(FATAL_ERROR " FFmpeg libswscale library could not be found!")
    ENDIF()

    SET(MOD_LIBRARIES
        #-lbz2
        ${FFMPEG_avcodec_LIBRARY}
        ${FFMPEG_avformat_LIBRARY}
        ${FFMPEG_avutil_LIBRARY}
        ${FFMPEG_swscale_LIBRARY}
    )
    SET(MOD_INCLUDE_DIRECTORIES ${FFMPEG_INCLUDE_DIR})

    MARK_AS_ADVANCED(FFMPEG_INCLUDE_DIR FFMPEG_avcodec_LIBRARY FFMPEG_avformat_LIBRARY FFMPEG_avutil_LIBRARY FFMPEG_swscale_LIBRARY )
ENDIF()


################################################################################
# Core module resources
################################################################################
SET(MOD_CORE_MODULECLASS FFmpegModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/videoencoder/videoencoder.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/videoencoder/videoencoder.h
)

