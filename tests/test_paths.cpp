
#include <string>
#include <cstdlib>
#include <iostream>
#include "kappa.hpp"
#include <catch.hpp>
#include "yaml-cpp/yaml.h"

using namespace std;
using namespace kappa;
using namespace Catch;

TEST_CASE("Test database path", "[database]")
{

  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
  std::string particle_source    = m_source + "particles.yaml";
  std::string interaction_source = m_source + "interaction.yaml";

  std::size_t part = particle_source.find("particles");
  std::size_t inte = interaction_source.find("interaction");

  SECTION("InvalidSourcePaths") {
    CHECK(particle_source.substr(part) == "particles.yaml");
    CHECK(interaction_source.substr(inte) == "interaction.yaml");
  }

}
