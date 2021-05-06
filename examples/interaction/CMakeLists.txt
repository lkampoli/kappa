cmake_policy(SET CMP0048 NEW)
cmake_minimum_required(VERSION 3.0)

get_filename_component(example_name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

project( KAPPAExamples LANGUAGES CXX )
message( STATUS "Examples included" )

add_executable( cs_dissTest
  cs_dissTest.cpp
)
add_executable( k_diss_ILTtest
  k_diss_ILTtest.cpp
)
add_executable( k_dissTest
  k_dissTest.cpp
)
add_executable( k_exchArrhTest
  k_exchArrhTest.cpp
)
add_executable( omegaBCEtest
  omegaBCEtest.cpp
)
add_executable( omegaTest
  omegaTest.cpp
)
add_executable( test_FHO
  test_FHO.cpp
)
add_executable( vss_and_neutral_e
  vss_and_neutral_e.cpp
)
add_executable( vssTest
  vssTest.cpp
)

if(NOT TARGET kappa++)
    find_package(kappa++ REQUIRED)
endif()

set(ALL_EXAMPLE_TARGETS
    cs_dissTest
    k_diss_ILTtest
    k_dissTest
    k_exchArrhTest
    omegaBCEtest
    omegaTest
    test_FHO
    vss_and_neutral_e
    vssTest
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
