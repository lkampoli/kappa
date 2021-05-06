#include <string>
#include <cstdlib>
#include <iostream>
#include "kappa.hpp"
#include "yaml-cpp/yaml.h"
#include <catch.hpp>

using namespace std;
using namespace kappa;
using namespace Catch;

TEST_CASE("Initialization of particles objects", "[particles]")
{

  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
  std::string particle_source    = m_source + "particles.yaml";
  std::string interaction_source = m_source + "interaction.yaml";

  Particle Ar(    "Ar", particle_source);
  Particle C(     "C",  particle_source);
  Particle e(     "e-", particle_source);
  Particle N(     "N",  particle_source);
  Particle Nplus( "N+", particle_source);
  Particle O(     "O",  particle_source);
  Particle Oplus( "O+", particle_source);
  Particle Ominus("O-", particle_source);

  Atom Ar_atom(    "Ar", particle_source);
  Atom C_atom(     "C",  particle_source);
  Atom N_atom(     "N",  particle_source);
  Atom Nplus_atom( "N+", particle_source);
  Atom O_atom(     "O",  particle_source);
  Atom Oplus_atom( "O+", particle_source);
  Atom Ominus_atom("O-", particle_source);

  Molecule CO(    "CO",  true, false, particle_source);
  Molecule C2(    "C2",  true, false, particle_source);
  Molecule N2(    "N2",  true, false, particle_source);
  Molecule N2plus("N2+", true, false, particle_source);
  Molecule NO(    "NO",  true, false, particle_source);
  Molecule O2(    "O2",  true, false, particle_source);
  Molecule O2plus("O2+", true, false, particle_source);
   
//  SECTION("Particles") {
//    CHECK(Ar.name() == "Ar");
//  }

  CHECK(CO.num_electron_levels == 28);
  CHECK(CO.anharmonic_spectrum == true);
  CHECK(CO.rigid_rotator == false);

}
