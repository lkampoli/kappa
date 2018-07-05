cmake_minimum_required(VERSION 2.8)

add_sources(kappa++
	approximation.cpp
	approximation_multit.cpp
	approximation_onet.cpp)

install(FILES approximation.hpp 	DESTINATION include/kappa++)
install(FILES approximation_multit.hpp 	DESTINATION include/kappa++)
install(FILES approximation_onet.hpp 	DESTINATION include/kappa++)
install(FILES approximations.hpp 	DESTINATION include/kappa++)
install(FILES exceptions.hpp 		DESTINATION include/kappa++)
install(FILES kappa.hpp 		DESTINATION include/kappa++)
install(FILES logging.hpp 		DESTINATION include/kappa++)
install(FILES models.h 			DESTINATION include/kappa++)
