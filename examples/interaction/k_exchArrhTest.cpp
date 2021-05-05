/*
  \file k_exchArrhTest.cpp  
  \brief Computing Arrhenius law according to different models.
*/

#include <iomanip>
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

  std::cout << "Start Test for k_exch coefficients, loading particle data" << endl;

  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
  std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
  std::string particle_source    = m_source + "particles.yaml";
  std::string interaction_source = m_source + "interaction.yaml";
  std::string output_dir = GetCurrentWorkingDir();
  std::cout << "Current directory is: " << output_dir << std::endl;

  // Molecule and Atom creation
  Molecule Molecule1("N2", false, true, particle_source);
  Atom Atom1("O", particle_source);

  std::cout << "Loading interaction parameters" << endl;

  // Molecule + Atom interaction definition
  Interaction inter(Molecule1, Atom1, interaction_source);

  // Approximation level selection
  Approximation ApproximationTest{};

  // Output file
  ofstream outf;

  std::vector<double> T_vals = {500., 1000., 2000., 5000., 10000., 15000., 20000., 25000., 30000., 40000.};
  std::vector<int> i_vals = {0, 10, 15, 20, 25};
  std::vector<std::string> model_names = {"ArrheniusScanlon", "ArrheniusPark"};
  std::vector<kappa::models_k_exch> model_vals = {models_k_exch::model_k_exch_arrh_scanlon, models_k_exch::model_k_exch_arrh_park};
  double res;
  int j =0;
  outf.open(output_dir + "/EXCH/" + Molecule1.name + "_" + Atom1.name + ".txt");
  for (auto exch_model: model_vals) {
    for (auto T : T_vals) {
      for (auto i : i_vals) {
        outf << std::setw(15) << T << std::setw(20);
        outf << std::setw(15) << i << std::setw(20);
        outf << std::setw(15) << ApproximationTest.k_exch(T, Molecule1, Atom1, inter, i, models_k_exch::model_k_exch_arrh_scanlon) << std::setw(20);
        outf << std::setw(15) << ApproximationTest.k_exch(T, Molecule1, Atom1, inter, i, models_k_exch::model_k_exch_arrh_park) << std::setw(20);
        outf << std::endl;
      }
    }
    outf.close();
    j++;
  }
  return 0;
}
