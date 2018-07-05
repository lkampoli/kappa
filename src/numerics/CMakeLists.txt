cmake_minimum_required(VERSION 2.8)

add_sources(kappa++
	numeric.cpp)

install(FILES numeric.hpp DESTINATION include/kappa++)
install(FILES constants.h DESTINATION include/kappa++)
