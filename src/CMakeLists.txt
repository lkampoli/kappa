cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
# Add sources for a target
#  ADD_SOURCES(<target> <source1> [<source2> ...])
 function (add_sources target)
  # define the <target>_SRCS properties if necessary
  get_property(prop_defined GLOBAL PROPERTY ${target}_SRCS DEFINED)
  if(NOT prop_defined)
    define_property(GLOBAL PROPERTY ${target}_SRCS
      BRIEF_DOCS "Sources for the ${target} target"
      FULL_DOCS "List of source files for the ${target} target")
  endif()
  # create list of sources (absolute paths)
  set(SRCS)
  foreach(src IN LISTS ARGN)
    if(NOT IS_ABSOLUTE "${src}")
      get_filename_component(src "${src}" ABSOLUTE)
    endif()
    list(APPEND SRCS "${src}")
  endforeach()
  # append to global property
  set_property(GLOBAL APPEND PROPERTY "${target}_SRCS" "${SRCS}")
 endfunction()

 add_subdirectory(approximations)
 add_subdirectory(interactions)
 add_subdirectory(mixtures)
 add_subdirectory(numerics)
 add_subdirectory(particles)

 get_property(kappa++_SRCS GLOBAL PROPERTY kappa++_SRCS)
 add_library(kappa++ SHARED ${kappa++_SRCS})
 target_link_libraries(kappa++ ${YAML_LIBRARY})
 target_link_libraries(kappa++ ${ARMADILLO_LIBRARY})
 install(TARGETS kappa++ DESTINATION lib)

 # Evaluate coverage
 if(ENABLE_COVERAGE)
     add_coverage(kappa++)
 endif()
