/* 
  \file vss_and_neutral_e.cpp
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

  Approximation ApproximationTest{};

  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
  std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
  std::string particle_source    = m_source + "particles.yaml";
  std::string interaction_source = m_source + "interaction.yaml";

  Molecule N2("N2", true, true, particle_source);
  Particle e("e-", particle_source);
  Atom Np("N+", particle_source);
  Atom N("N", particle_source);

  Interaction N2Np(N2, Np, interaction_source);
  Interaction ee(e, e, interaction_source);
  Interaction Ne(N2, e, interaction_source);
  // Interaction N2N(N2, N);
  // Interaction HH(H, H);
  std::vector<int> l_arr = {1, 1, 2, 3};
  std::vector<int> r_arr = {1, 4, 2, 3};
  int l, r;

  std::cout << "collision-reduced mass of two electrons: " << ee.collision_mass << ", e- mass:" << e.mass << std::endl;
  std::cout << "what is it? " << Ne["_Om22_0"] << std::endl;

  for (int i=0;i<4; i++) {
    l = l_arr[i];
    r = r_arr[i];

    std::cout << N2.name << " + " << Np.name << std::endl;
    std::cout << "ESA: " << ApproximationTest.omega_integral(2000., N2Np, l, r, models_omega::model_omega_esa) << " " << ApproximationTest.omega_integral(10000., N2Np, l, r, models_omega::model_omega_esa) << std::endl;
    std::cout << "RS: " << ApproximationTest.omega_integral(2000., N2Np, l, r, models_omega::model_omega_rs) << " " << ApproximationTest.omega_integral(10000., N2Np, l, r, models_omega::model_omega_rs) << std::endl;
    std::cout << "VSS: " << ApproximationTest.omega_integral(2000., N2Np, l, r, models_omega::model_omega_vss) << " " << ApproximationTest.omega_integral(10000., N2Np, l, r, models_omega::model_omega_vss) << std::endl << std::endl;

    std::cout << N.name << " + " << e.name << std::endl;
    std::cout << "ESA: " << ApproximationTest.omega_integral(2000., Ne, l, r, models_omega::model_omega_esa) << " " << ApproximationTest.omega_integral(10000., Ne, l, r, models_omega::model_omega_esa) << std::endl;
    std::cout << "RS: " << ApproximationTest.omega_integral(2000., Ne, l, r, models_omega::model_omega_rs) << " " << ApproximationTest.omega_integral(10000., Ne, l, r, models_omega::model_omega_rs) << std::endl;
    std::cout << std::endl;
  }

  return 0;
}
