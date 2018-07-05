#
# this module look for Armadillo support
# it will define the following values
#
# ARMADILLO_INCLUDE_DIR = where armadillo.h can be found
# ARMADILLO_LIBRARY     = the library to link against (armadillo etc)
# KAPPA_HAS_ARMADILLO   = set to true after finding the library
#
 FIND_PATH(ARMADILLO_INCLUDE_DIR armadillo ${CMAKE_CURRENT_SOURCE_DIR}/kappa_deps/x86_64/include )
 FIND_LIBRARY(ARMADILLO_LIBRARY armadillo ${CMAKE_CURRENT_SOURCE_DIR}/kappa_deps/x86_64/lib )
 FIND_LIBRARY(ARMADILLO_LIBRARY armadillo ${CMAKE_CURRENT_SOURCE_DIR}/kappa_deps/x86_64/lib64 )

 IF(ARMADILLO_INCLUDE_DIR AND ARMADILLO_LIBRARY)
  SET(KAPPA_HAS_ARMADILLO 1 CACHE BOOL "Found armadillo library")
 ELSE()
  SET(KAPPA_HAS_ARMADILLO 0 CACHE BOOL "Not fount armadillo library")
 ENDIF()

 MARK_AS_ADVANCED(
  ARMADILLO_INCLUDE_DIR
  ARMADILLO_LIBRARY
  KAPPA_HAS_ARMADILLO
 )

 MESSAGE ( " KAPPA_HAS_ARMADILLO:   [${KAPPA_HAS_ARMADILLO}]"   )
 MESSAGE ( " ARMADILLO_INCLUDE_DIR: [${ARMADILLO_INCLUDE_DIR}]" )
 MESSAGE ( " ARMADILLO_LIBRARY:     [${ARMADILLO_LIBRARY}]"     )
