
#include <string>
#include <cstdlib>
#include <iostream>
#include "kappa.hpp"
#include <catch.hpp>
#include "yaml-cpp/yaml.h"

using namespace std;
using namespace kappa;
using namespace Catch;

TEST_CASE("Test basic", "[basic]")
{

  const double tol = std::numeric_limits<double>::epsilon();

  YAML::Node node = YAML::Load("{pi: 3.14159, [0, 1]: integers}");
  SECTION("YAML") {
     CHECK(node["pi"].as<double>() == 3.14159);
  }

  arma::vec v1 = arma::ones(5);
  arma::vec v2 = arma::ones(5);
  v2[0] = 4;

  SECTION("Armadillo") {
    CHECK(arma::dot(v1, v2) == 8);
//  CHECK(kappa::factorial_table[5] == 720);
//  CHECK(kappa::factorial_table[69] == 4.94066e-324);
  }

  auto test_f = [](int i) {return 0.5 + i;}; 
  auto test_f_integration = [](double x) {return x; }; 

  SECTION("KappaNumerics") {
//  CHECK(kappa::integrate_interval(test_f_integration, 0., 1.) == 0.5);
    CHECK(kappa::find_max_value(test_f, 7.) == 6);
  }

}
