#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8.11)

#Name your project here
project(netmets)

#set the module directory
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${VRN_NETMETS_DIRECTORY}")

#find the STIM library
find_package(STIM)

#if the STIM library isn't found, download it
if(NOT STIM_INCLUDE_DIRS)
	file(REMOVE_RECURSE ${VRN_NETMETS_DIRECTORY}/stimlib)	#remove the stimlib directory if it exists
	set(STIM_GIT "https://git.stim.ee.uh.edu/codebase/stimlib.git")
	execute_process(COMMAND git clone --depth 1 ${STIM_GIT} WORKING_DIRECTORY ${VRN_NETMETS_DIRECTORY})
	set(STIM_INCLUDE_DIRS "${VRN_NETMETS_DIRECTORY}/stimlib" CACHE TYPE PATH)
endif(NOT STIM_INCLUDE_DIRS)

#find BOOST
find_package(Boost REQUIRED)

#find cuda
find_package(CUDA REQUIRED)

#find the GLUT library for visualization
#find_package(OpenGL REQUIRED)

#find_package(GLUT REQUIRED)

#GLEW is required for MS Visual C
#if(MSVC)
#	find_package(GLEW REQUIRED)
#	include_directories(${GLEW_INCLUDE_DIR})
#endif(MSVC)

#find OpenCV
#find_package(OpenCV REQUIRED)
#add_definitions(-DUSING_OPENCV)

#find the pthreads package
#find_package(Threads)

#find the X11 package
#find_package(X11)

#find the Approximate Nearest Neighbor Library
#find_package(ANN REQUIRED)

#build the executable in the binary directory on MS Visual Studio
#if ( MSVC )
#	SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${OUTPUT_DIRECTORY}")
#	SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${OUTPUT_DIRECTORY}")
#endif ( MSVC )

MESSAGE(STATUS "VRN ROOT INCLUDE DIR: "${VRN_NETMETS_VRN_ROOT_DIR})
include_directories(
	#${OPENGL_INCLUDE_DIRS}
	#				${GLUT_INCLUDE_DIR}
					${STIM_INCLUDE_DIRS}
          			#${ANN_INCLUDE_DIR}
					${Boost_INCLUDE_DIR}
					${VRN_NETMETS_VRN_ROOT_DIR}/include
					${VRN_NETMETS_VRN_ROOT_DIR}/ext
					${VRN_NETMETS_VRN_ROOT_DIR}/ext/rapidjson
					#${OpenCV_INCLUDE_DIRS}
					)

ADD_DEFINITIONS("-DVRN_VESSELNETWORKANALYSISEXTRA_MINIMAL_VESSELGRAPH")
#Assign source files to the appropriate variables
#file(GLOB SRC_CPP "*.cu")
#file(GLOB SRC_H "*.h")

#create an executable file
cuda_add_library(netmets SHARED lib.cu)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#set the link libraries
target_link_libraries(netmets
#						${OpenGL_LIBRARIES}
#${OPENGL_gl_LIBRARY}
#${OPENGL_glu_LIBRARY}
#${GLUT_LIBRARIES}
#${GLEW_GLEW_LIBRARY}
#${CMAKE_THREAD_LIBS_INIT}
            			#${ANN_LIBRARY}
				#${X11_LIBRARIES}
				#${OpenCV_LIBS}
		  )
	  #if(MSVC)
	  #	target_link_libraries(netmets ${GLEW_LIBRARY})
	  #endif(MSVC)
