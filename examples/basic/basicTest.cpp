/*
   \file basicTest.cpp
   \brief some basic tests - YAML parsing, linear algebra, utility numeric functions
*/

#include <iostream>
#include "yaml-cpp/yaml.h"
#include <kappa.hpp>

int main(void) {

  YAML::Node node = YAML::Load("{pi: 3.14159, [0, 1]: integers}");
  std::cout << node["pi"].as<double>() << "\n";

  arma::vec v1 = arma::ones(5);
  arma::vec v2 = arma::ones(5);
  v2[0] = 4;
  std::cout << "(v1, v2) = " << arma::dot(v1, v2) << "\n";
  std::cout << "5! = " << kappa::factorial_table[5] << " 69!=" << kappa::factorial_table[69] << "\n";

  auto test_f = [](int i) {return 0.5 + i;}; // find_max_value should give 6
  auto test_f_integration = [](double x) {return x; }; // integral from 0 to 1 should be 0.5

  std::cout << "int(f)_0^1 = " << kappa::integrate_interval(test_f_integration, 0, 1.) << "\n";
  std::cout << "Max i = " << kappa::find_max_value(test_f, 7.) << "\n";

  return 0;
}
