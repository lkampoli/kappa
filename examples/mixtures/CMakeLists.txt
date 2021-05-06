cmake_policy(SET CMP0048 NEW)
cmake_minimum_required(VERSION 3.0)

get_filename_component(example_name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

project( KAPPAExamples LANGUAGES CXX )
message( STATUS "Examples included" )

#add_executable( ionized_mixture_creation
#  ionized_mixture_creation.cpp
#)
add_executable( mixture_string_creation
  mixture_string_creation.cpp
)
add_executable( mixture-sts-basic
  mixture-sts-basic.cpp
)
add_executable( mixture-sts-basic_boltzman
  mixture-sts-basic_boltzman.cpp
)
add_executable( mixture-sts-bulk
  mixture-sts-bulk.cpp
)
add_executable( mixture-sts-bulk_air5
  mixture-sts-bulk_air5.cpp
)
add_executable( mixture-sts-shear
  mixture-sts-shear.cpp
)
add_executable( mixture-sts-shear-bulk-thermal
  mixture-sts-shear-bulk-thermal.cpp
)
add_executable( mixture-sts-thermal_conductivity
  mixture-sts-thermal_conductivity.cpp
)
add_executable( mixture-sts-thermo-diffusion
  mixture-sts-thermo-diffusion.cpp
)
add_executable( TestDiffusion
  TestDiffusion.cpp
)
add_executable( thermal_conductivity
  thermal_conductivity.cpp
)
add_executable( thermodiffusion
  thermodiffusion.cpp
)
add_executable( transport_coefficients_air5
  transport_coefficients_air5.cpp
)
add_executable( mk_TC_air5_database
  mk_TC_air5_database.cpp
)

if(NOT TARGET kappa++)
    find_package(kappa++ REQUIRED)
endif()

set(ALL_EXAMPLE_TARGETS
    #ionized_mixture_creation
    mixture_string_creation
    mixture-sts-basic
    mixture-sts-basic_boltzman
    mixture-sts-bulk
    mixture-sts-bulk_air5
    mixture-sts-shear
    mixture-sts-shear-bulk-thermal
    mixture-sts-thermal_conductivity
    mixture-sts-thermo-diffusion
    TestDiffusion
    thermal_conductivity
    thermodiffusion
    transport_coefficients_air5
    mk_TC_air5_database
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
