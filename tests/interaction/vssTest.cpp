/*
  \file vssTest.cpp
*/

#include <iostream>
#include <fstream>

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#include "kappa.hpp"

std::string GetCurrentWorkingDir( void ) {
  char buff[FILENAME_MAX];
  GetCurrentDir( buff, FILENAME_MAX );
  std::string current_working_dir(buff);
  return current_working_dir;
}

using namespace std;
using namespace kappa;

int main(int argc, char** argv) {

  std::cout << "Start Test for Omega integrals, loading particle data" << endl;
  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
  std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << std::endl;
  std::string particle_source    = m_source + "particles.yaml";
  std::string interaction_source = m_source + "interaction.yaml";
  std::string output_dir = GetCurrentWorkingDir();
  std::cout << "Current directory is: " << output_dir << std::endl;

  Approximation ApproximationTest{};

  Molecule N2("N2", true, true, particle_source);
  Atom N("N", particle_source);
  // Atom H("H", particle_source);

  Interaction N2N2(N2, N2, interaction_source);
  Interaction N2N(N2, N, interaction_source);
  // Interaction HH(H, H, interaction_source);

  std::cout << "ESA: " << ApproximationTest.omega_integral(2000., N2N2, 1, 1, models_omega::model_omega_esa) << " " << ApproximationTest.omega_integral(10000., N2N2, 1, 1, models_omega::model_omega_esa) << std::endl;
  std::cout << "RS: " << ApproximationTest.omega_integral(2000., N2N2, 1, 1, models_omega::model_omega_rs) << " " << ApproximationTest.omega_integral(10000., N2N2, 1, 1, models_omega::model_omega_rs) << std::endl;
  std::cout << "VSS: " << ApproximationTest.omega_integral(2000., N2N2, 1, 1, models_omega::model_omega_vss) << " " << ApproximationTest.omega_integral(10000., N2N2, 1, 1, models_omega::model_omega_vss) << std::endl << std::endl;

  std::cout << "ESA: " << ApproximationTest.omega_integral(2000., N2N2, 2, 3, models_omega::model_omega_esa) << " " << ApproximationTest.omega_integral(10000., N2N2, 2, 3, models_omega::model_omega_esa) << std::endl;
  std::cout << "RS: " << ApproximationTest.omega_integral(2000., N2N2, 2, 3, models_omega::model_omega_rs) << " " << ApproximationTest.omega_integral(10000., N2N2, 2, 3, models_omega::model_omega_rs) << std::endl;
  std::cout << "VSS: " << ApproximationTest.omega_integral(2000., N2N2, 2, 3, models_omega::model_omega_vss) << " " << ApproximationTest.omega_integral(10000., N2N2, 2, 3, models_omega::model_omega_vss) << std::endl << std::endl;

  std::cout << "ESA: " << ApproximationTest.omega_integral(2000., N2N, 1, 1, models_omega::model_omega_esa) << " " << ApproximationTest.omega_integral(10000., N2N, 1, 1, models_omega::model_omega_esa) << std::endl;
  std::cout << "RS: " << ApproximationTest.omega_integral(2000., N2N, 1, 1, models_omega::model_omega_rs) << " " << ApproximationTest.omega_integral(10000., N2N, 1, 1, models_omega::model_omega_rs) << std::endl;
  std::cout << "VSS: " << ApproximationTest.omega_integral(2000., N2N, 1, 1, models_omega::model_omega_vss) << " " << ApproximationTest.omega_integral(10000., N2N, 1, 1, models_omega::model_omega_vss) << std::endl << std::endl;

  std::cout << "ESA: " << ApproximationTest.omega_integral(2000., N2N, 2, 3, models_omega::model_omega_esa) << " " << ApproximationTest.omega_integral(10000., N2N, 2, 3, models_omega::model_omega_esa) << std::endl;
  std::cout << "RS: " << ApproximationTest.omega_integral(2000., N2N, 2, 3, models_omega::model_omega_rs) << " " << ApproximationTest.omega_integral(10000., N2N, 2, 3, models_omega::model_omega_rs) << std::endl;
  std::cout << "VSS: " << ApproximationTest.omega_integral(2000., N2N, 2, 3, models_omega::model_omega_vss) << " " << ApproximationTest.omega_integral(10000., N2N, 2, 3, models_omega::model_omega_vss) << std::endl << std::endl;

  // try {
  //   std::cout << ApproximationTest.omega_integral(2000., HH, 2, 3, models_omega::model_omega_vss);
  // } catch (const DataNotFoundException &e) {
  //   std::cout << e.what() << endl;
  // }

  return 0;
}
