cmake_policy(SET CMP0048 NEW)
cmake_minimum_required(VERSION 3.0)

get_filename_component(example_name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

project( KAPPAExamples LANGUAGES CXX )
message( STATUS "Examples included" )

add_executable( cvibr_multit
  cvibr_multit.cpp
)
add_executable( cvibr_sts
  cvibr_sts.cpp
)
add_executable( cvrot_sts
  cvrot_sts.cpp
)
add_executable( cvtr_sts
  cvtr_sts.cpp
)

if(NOT TARGET kappa++)
    find_package(kappa++ REQUIRED)
endif()

set(ALL_EXAMPLE_TARGETS
    cvibr_multit
    cvibr_sts
    cvrot_sts
    cvtr_sts
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
