#
# this module look for YAML support
# it will define the following values
#
# YAML_INCLUDE_DIR  = where yaml-cpp.h can be found
# YAML_LIBRARY      = the library to link against (yaml-cpp etc)
# KAPPA_HAS_YAML    = set to true after finding the library
#
 FIND_PATH(YAML_INCLUDE_DIR yaml.h ${CMAKE_CURRENT_SOURCE_DIR}/kappa_deps/x86_64/include)
 FIND_PATH(YAML_INCLUDE_DIR yaml.h ${CMAKE_CURRENT_SOURCE_DIR}/kappa_deps/x86_64/include/yaml-cpp) 
 FIND_LIBRARY(YAML_LIBRARY yaml-cpp ${CMAKE_CURRENT_SOURCE_DIR}/kappa_deps/x86_64/lib)
 FIND_LIBRARY(YAML_LIBRARY yaml-cpp ${CMAKE_CURRENT_SOURCE_DIR}/kappa_deps/x86_64/lib64)

 IF(YAML_INCLUDE_DIR AND YAML_LIBRARY)
  SET(KAPPA_HAS_YAML 1)
 ELSE()
  SET(KAPPA_HAS_YAML 0 CACHE BOOL "Not found yaml library")
 ENDIF()

 MARK_AS_ADVANCED(
  YAML_INCLUDE_DIR
  YAML_LIBRARY
  KAPPA_HAS_YAML
 )

 MESSAGE ( " KAPPA_HAS_YAML:   [${KAPPA_HAS_YAML}]"   )
 MESSAGE ( " YAML_INCLUDE_DIR: [${YAML_INCLUDE_DIR}]" )
 MESSAGE ( " YAML_LIBRARY:     [${YAML_LIBRARY}]"     )
