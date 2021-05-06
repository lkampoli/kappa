cmake_minimum_required(VERSION 2.6)

include_directories(${CMAKE_CURRENT_BINARY_DIR})
configure_file(Configuration.h.in Configuration.h)

file(GLOB test_sources "test_*.cpp")
add_executable(run_tests run_tests.cpp ${test_sources})
target_link_libraries(run_tests
    PRIVATE
        kappa++
    PUBLIC
        Catch2::Catch2
)

# Add all Catch TEST_CASEs as individual CTest tests
include(CTest)
include(ParseAndAddCatchTests)
set(PARSE_CATCH_TESTS_VERBOSE ON)
ParseAndAddCatchTests(run_tests)

#add_executable(update_comparison update_comparison.cpp)
#target_link_libraries(update_comparison kappa++)
