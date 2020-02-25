# This file contains cmake code that will be executed at the very end of the cmake script processing,
# especially after the install target has finished it's execution.
# Due to technical limitations of CMake, we are forced to include the whole directory for that purpose
# which will guarantee the desired order of execution.

MESSAGE(STATUS "Processing post-install script")

# Execute deployment script for linux at the very end, because everything has to be copied properly before.

IF(UNIX AND VRN_DEPLOYMENT)
    OPTION(VRN_DEPLOYMENT_CREATE_APPIMAGE "Create an AppImage inside the installation directory" OFF)
    IF(VRN_DEPLOYMENT_CREATE_APPIMAGE)
        MESSAGE(STATUS "* Creating AppImage")
        INSTALL(CODE "EXECUTE_PROCESS(COMMAND ${VRN_HOME}/tools/snapshot/deployment-linux/linuxdeployment.sh ${CMAKE_INSTALL_PREFIX} ${VRN_HOME}/tools/snapshot/deployment-linux/)")
    ENDIF()
ENDIF()
