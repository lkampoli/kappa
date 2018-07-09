/*! 
    \file: mixture-sts-shear.cpp 
    \brief Test for shear viscosity in binary mixtures.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip> 

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#include "kappa.hpp"

using namespace kappa;

std::string GetCurrentWorkingDir( void ) {
  char buff[FILENAME_MAX];
  GetCurrentDir( buff, FILENAME_MAX );
  std::string current_working_dir(buff);
  return current_working_dir;
}

int main(int argc, char** argv) {

  std::cout << "Start test: computation of shear viscosity" << std::endl;

  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
  std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
  std::string particle_source    = m_source + "particles.yaml";
  std::string interaction_source = m_source + "interaction.yaml";
  std::string output_dir = GetCurrentWorkingDir();
  std::cout << "Current directory is: " << output_dir << std::endl;

  std::vector<kappa::Molecule> molecules;
  std::vector<kappa::Atom> atoms;

  std::cout << "Loading particles data" << std::endl;

  kappa::Molecule mol("NO", false, true, particle_source);
  // kappa::Molecule mol("N2", false, true, particle_source);
  // kappa::Atom at("N", particle_source);
  // kappa::Molecule mol("O2", false, true, particle_source);
  // kappa::Atom at("O", particle_source);
  kappa::Atom N("N", particle_source);
  atoms.push_back(N);
  kappa::Atom O("O", particle_source);
  atoms.push_back(O);
  molecules.push_back(mol);
  // atoms.push_back(at);

  kappa::Mixture mixture(molecules, atoms, interaction_source, particle_source);
  std::cout << "Mixture created" << std::endl;

  // std::vector<double> T_vals = { 500., 1000., 5000., 10000., 20000., 40000. };
  std::vector<double> T_vals;
  // double x_atom = 0.0;
  double x_N2 = 1.0;
  double x_NO = 1.0;
  double x_N = 0.0;
  double x_O = 0.0;
  double pressure = 101325.;
  double n_tot = pressure / (K_CONST_K * 500);

  // for (int i=0; i<80; i++) { //500-40000K
  //   T_vals.push_back(500 + i * 500);
  // }

  for (int i=0; i<50; i++) { 
    T_vals.push_back(100 + i * 50);
  }

  std::vector<arma::vec> mol_ndens;
  mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], x_NO * 101325.0 / (K_CONST_K * T_vals[0]), mol)); // NO

  std::ofstream outf;
  outf.open(output_dir + "/TRANSPORT_COEFFICIENTS/shear_viscosity/" + mol.name + "_xat_" + std::to_string(x_NO) + ".txt");
  outf << std::setw(20) << "Temperature [K]";
  outf << std::setw(20) << "Eta";
  outf << std::endl;

  arma::vec atom_ndens(2);
  atom_ndens[0] = 0.;
  atom_ndens[1] = 0.;

  double tot_ndens;

  for (auto T : T_vals) {
    tot_ndens =  101325.0 / (K_CONST_K * T);
    mol_ndens[0] = mixture.Boltzmann_distribution(T, x_NO * 101325.0 / (K_CONST_K * T), mol);

    // mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens);
    mixture.compute_transport_coefficients(T, mol_ndens);

    outf << std::setw(20) << T;
    outf << std::setw(20) << mixture.get_shear_viscosity();
    outf << std::endl;
  }

  outf.close();
  return 0;
}
