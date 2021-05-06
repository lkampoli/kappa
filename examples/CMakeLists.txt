cmake_minimum_required(VERSION 2.8.2)

# Append the PROJECT_BINARY_DIR to CMAKE_PREFIX_PATH in order to be able to
# find the configuration files at build time
list(APPEND CMAKE_PREFIX_PATH ${PROJECT_BINARY_DIR}/src)

# c++ examples
#add_subdirectory(c++)

# Fortran examples
#if (BUILD_FORTRAN_WRAPPER)
#    add_subdirectory(fortran)
#endif(BUILD_FORTRAN_WRAPPER)

#include(ListSubDirs)
#ListSubDirs(example_subdirs ${CMAKE_CURRENT_LIST_DIR})

# Add all examples targets and test them
add_subdirectory(basic)
add_subdirectory(particles)
add_subdirectory(approximations)
add_subdirectory(interaction)
add_subdirectory(mixtures)
#foreach(subdir ${example_subdirs})
#    add_subdirectory(${subdir})
#    add_test(
#        NAME example_${subdir}
#        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/${subdir}
#        COMMAND ${subdir} )
#endforeach()
