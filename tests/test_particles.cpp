/* 
 * \file particleTest.cpp
 */

#include <cstdlib>
#include <iostream>

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#include "kappa.hpp"

#include <catch.hpp>
/*
std::string GetCurrentWorkingDir( void ) {
  char buff[FILENAME_MAX];
  GetCurrentDir( buff, FILENAME_MAX );
  std::string current_working_dir(buff);
  return current_working_dir;
}
*/
using namespace std;
using namespace kappa;
using namespace Catch;

TEST_CASE("Initialization of particles objects", "[particles]")
{

  cout << "Start Test for Particle classes" << endl;
  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
  std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
  std::string particle_source    = m_source + "particles.yaml";
  std::string interaction_source = m_source + "interaction.yaml";

  Molecule CO("CO", true, false, particle_source);

  CHECK(CO.num_electron_levels == 28);
  CHECK(CO.anharmonic_spectrum == true);
  CHECK(CO.rigid_rotator == false);

}
