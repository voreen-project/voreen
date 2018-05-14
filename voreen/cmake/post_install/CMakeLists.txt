# This file contains cmake code that will be executed at the very end of the cmake script processing,
# especially after the install target has finished it's execution.

# Execute deployment script for linux at the very end, because everything has to be copied properly before.
IF(UNIX AND VRN_DEPLOYMENT)
    INSTALL(CODE "EXECUTE_PROCESS(COMMAND ${VRN_HOME}/tools/snapshot/deployment-linux/linuxdeployment.sh ${CMAKE_INSTALL_PREFIX} ${VRN_HOME}/tools/snapshot/deployment-linux/)")
ENDIF()