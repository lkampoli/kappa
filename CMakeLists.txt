 cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
 project(kappa++ CXX)

################################################################################
# Build types
################################################################################

# Profile
 set(CMAKE_CXX_FLAGS_PROFILE "-g3 -Wall -O3 -DNDEBUG" CACHE STRING
    "Flags used by the C++ compiler during Profile builds."
    FORCE )
 set(CMAKE_C_FLAGS_PROFILE "-g3 -Wall -pedantic -O3 -DNDEBUG" CACHE STRING
    "Flags used by the C compiler during Profile builds."
    FORCE )
 set(CMAKE_EXE_LINKER_FLAGS_PROFILE
    "" CACHE STRING
    "Flags used for linking binaries during Profile builds."
    FORCE )
 set(CMAKE_SHARED_LINKER_FLAGS_PROFILE
    "" CACHE STRING
    "Flags used by the shared libraries linker during Profile builds."
    FORCE )
 mark_as_advanced(
    CMAKE_CXX_FLAGS_PROFILE
    CMAKE_C_FLAGS_PROFILE
    CMAKE_EXE_LINKER_FLAGS_PROFILE
    CMAKE_SHARED_LINKER_FLAGS_PROFILE)

# Update the documentation string of CMAKE_BUILD_TYPE for GUIs and set the
# default build type as Release

 if (NOT CMAKE_BUILD_TYPE)
     SET (CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: None Debug Release 
        RelWithDebInfo MinSizeRel Profile." FORCE)
 endif ()

################################################################################
# Install prefix settings
################################################################################

 set(CMAKE_INSTALL_PREFIX "../install" CACHE STRING
    "Install path prefix, prepended onto install directories." FORCE)
 mark_as_advanced(CMAKE_INSTALL_PREFIX)

 if (CMAKE_COMPILER_IS_GNUCXX)

    #set (CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wold-style-cast")
    #set (CMAKE_CXX_FLAGS "-g")

    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -lm -g -pg -fopenmp 						" )
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I ${CMAKE_CURRENT_SOURCE_DIR}/kappa_deps/x86_64/include/				" )
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I ${CMAKE_CURRENT_SOURCE_DIR}/kappa_deps/x86_64/include/yaml-cpp/			" )
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I ${CMAKE_CURRENT_SOURCE_DIR}/kappa_deps/x86_64/include/armadillo_bits/		" )
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I ${CMAKE_CURRENT_SOURCE_DIR}/kappa_deps/x86_64/include/boost/			" )
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I ${CMAKE_CURRENT_SOURCE_DIR}/src/approximations/					" )
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I ${CMAKE_CURRENT_SOURCE_DIR}/src/interactions/					" )
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I ${CMAKE_CURRENT_SOURCE_DIR}/src/numerics/					" )
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I ${CMAKE_CURRENT_SOURCE_DIR}/src/mixtures/					" )
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I ${CMAKE_CURRENT_SOURCE_DIR}/src/particles					" )

 endif ()

################################################################################
# Doxygen documentation generation
################################################################################
 option (BUILD_DOCUMENTATION 
    "Use Doxygen to create the HTML based API documentation" OFF)

 if (BUILD_DOCUMENTATION)
    FIND_PACKAGE(Doxygen)

    if (NOT DOXYGEN_FOUND)
        message(FATAL_ERROR 
            "Doxygen is needed to build the documentation. Please install it 
             correctly")
    endif()
    # Configure the Template Doxyfile for our specific project
    configure_file(Doxyfile.in 
                    ${PROJECT_BINARY_DIR}/Doxyfile  @ONLY IMMEDIATE)
    # Add a custom target to run Doxygen when ever the project is built
    add_custom_target (docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
        SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)
 endif()

################################################################################
# Source code
################################################################################
# Descend into the src directory to build all targets and libraries
# get all *.cpp files recursively

# file(GLOB_RECURSE SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
 include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/approximations
    ${CMAKE_CURRENT_SOURCE_DIR}/src/interactions
    ${CMAKE_CURRENT_SOURCE_DIR}/src/mixtures
    ${CMAKE_CURRENT_SOURCE_DIR}/src/numerics
    ${CMAKE_CURRENT_SOURCE_DIR}/src/particles)
 add_subdirectory(src)

#################################################################################
# Find the required packages: OpenBLAS(Blas, Lapack), Armadillo, Yaml-cpp, Boost
#################################################################################

 set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# SET(Boost_LIBRARIES ${Boost_LIBRARIES} "${CMAKE_CURRENT_SOURCE_DIR}/kappa_deps/x86_64/lib")
# find_package(Boost REQUIRED)
# set_property(TARGET boost PROPERTY
# add_library(boost INTERFACE IMPORTED)
# INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})
 
 find_package(Yaml REQUIRED)
  add_library(Yaml INTERFACE IMPORTED)
  set_property(TARGET Yaml PROPERTY
  INTERFACE_INCLUDE_DIRECTORIES ${YAML_INCLUDE_DIR})
  include_directories(${YAML_INCLUDE_DIRS})
  LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/kappa_deps/x86_64/lib)

 find_package(BlasLapack) 

 find_package(Armadillo REQUIRED)
  add_library(Armadillo INTERFACE IMPORTED)
  set_property(TARGET Armadillo PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${ARMADILLO_INCLUDE_DIR})
  include_directories(${ARMADILLO_INCLUDE_DIRS})
  LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/kappa_deps/x86_64/lib64)

################################################################################
# Install the data directory in the install folder
################################################################################
 install(FILES data/interaction.yaml     DESTINATION include/kappa++_data)
 install(FILES data/particles.yaml       DESTINATION include/kappa++_data)
################################################################################
# CTest
################################################################################
 option (ENABLE_TESTING
    "Enable regression testing with CTest (requires Boost.Test library)" OFF)
 if (ENABLE_TESTING)
    enable_testing()
    add_subdirectory(tests)
 endif()
