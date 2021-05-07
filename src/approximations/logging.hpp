/*!
    \file logging.hpp
 */

#ifndef kappa_logging_h
#define kappa_logging_h

#include <fstream>
#include <string>
#include <arma>

namespace kappa {

// 0 - no logging; 
// 1 - temperatures, transport coefficients, summed densities;
// 2 - matrices (summed), RHS (summed), 
// 3 - matrices, RHS, non-summed densities
enum logging_sts {no_logging, logging_basic, logging_detailed_1, logging_detailed_2, logging_detailed_sts}; 

void log_new_step(std::ofstream &logfile, int scalar);
void log_scalar(std::ofstream &logfile, std::string name, double scalar, int offset);
void log_vector(std::ofstream &logfile, std::string name, arma::vec vector, int offset);
void log_matrix(std::ofstream &logfile, std::string name, arma::mat matrix, int offset);
}
#endif /* models_h */
