cmake_policy(SET CMP0048 NEW)
cmake_minimum_required(VERSION 3.0)

get_filename_component(example_name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

project( KAPPAExamples LANGUAGES CXX )
message( STATUS "Examples included" )

add_executable( dumpspectrumTest
  dumpspectrumTest.cpp
)
add_executable( particleTest
  particleTest.cpp
)

if(NOT TARGET kappa++)
    find_package(kappa++ REQUIRED)
endif()

set(ALL_EXAMPLE_TARGETS
    dumpspectrumTest
    particleTest
)

foreach( name ${ALL_EXAMPLE_TARGETS} )
    target_link_libraries( ${name}
        INTERFACE
            Yaml
            Armadillo
        PRIVATE
            kappa++
    )
endforeach()
