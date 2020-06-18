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

  //kappa::Molecule mol("NO", false, true, particle_source);
  kappa::Molecule mol("N2", true, false, particle_source);
  kappa::Atom at("N", particle_source);
  // kappa::Molecule mol("O2", false, true, particle_source);
  // kappa::Atom at("O", particle_source);
  //kappa::Atom N("N", particle_source);
  //atoms.push_back(N);
  //kappa::Atom O("O", particle_source);
  //atoms.push_back(O);
  molecules.push_back(mol);
  atoms.push_back(at);

  kappa::Mixture mixture(molecules, atoms, interaction_source, particle_source);
  std::cout << "Mixture created" << std::endl;

  //std::vector<double> T_vals = { 500., 1000., 5000., 10000., 20000., 40000. };
  //std::vector<double> T_vals = { 500. };
  std::vector<double> T_vals = { 9971.900000 };
  //std::vector<double> T_vals;
  // double x_atom = 0.0;
  double x_N2 = 0.90;
  double x_NO = 0.0;
  double x_N = 0.10;
  double x_O = 0.0;
  double pressure = 27358.000000; //101325.;
  //double n_tot = pressure / (K_CONST_K * 500);

  arma::vec x_mass_fractions = arma::zeros(molecules.size()+atoms.size());
  x_mass_fractions[0] = 9.991800E-01; // N2
  x_mass_fractions[1] = 1.252500E-03; // N

  // mass to molar conversion
  arma::vec x_molar_fractions = arma::zeros(molecules.size()+atoms.size());
  x_molar_fractions = mixture.convert_mass_frac_to_molar(x_mass_fractions);
  std::cout << std::setw(20) << x_molar_fractions << "\n";

  x_N2 = x_molar_fractions[0]; //0.9975
  x_N = x_molar_fractions[1]; //0.0025
  std::cout << x_N2 << " " << x_N << "\n";

  // for (int i=0; i<80; i++) { //500-40000K
  //   T_vals.push_back(500 + i * 500);
  // }

  //for (int i=0; i<50; i++) { 
  //  T_vals.push_back(100 + i * 50);
  //}

  std::vector<arma::vec> mol_ndens;
  mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], x_N2 * pressure / (K_CONST_K * T_vals[0]), mol)); // N2

  std::ofstream outf;
  outf.open(output_dir + "/TRANSPORT_COEFFICIENTS/shear_viscosity/" + mol.name + "_xat_" + std::to_string(x_N2) + ".txt");
  outf << std::setw(20) << "Temperature [K]";
  outf << std::setw(20) << "Eta";
  outf << std::endl;

  arma::vec atom_ndens(2);
  atom_ndens[0] = 0.;
  atom_ndens[1] = 0.;

  double tot_ndens;

//  std::vector<models_omega> omega_integral_models = {models_omega::model_omega_rs, models_omega::model_omega_vss,
//                                                     models_omega::model_omega_bornmayer, models_omega::model_omega_lennardjones,
//                                                     models_omega::model_omega_esa};

  for (auto T : T_vals) {
    tot_ndens =  pressure / (K_CONST_K * T);
    mol_ndens[0] = mixture.Boltzmann_distribution(T, x_N2 * pressure / (K_CONST_K * T), mol);
    atom_ndens[0] = x_N * pressure / (K_CONST_K * T);

    // mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens);
    mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens, 0, models_omega::model_omega_rs);
    //mixture.compute_transport_coefficients(T, mol_ndens);

    outf << std::setw(20) << T;
    outf << std::setw(20) << mixture.get_shear_viscosity();
    outf << std::endl;
  }

  std::cout << std::setw(20) << " shear viscosity = " << mixture.get_shear_viscosity() << std::endl;
  std::cout << std::setw(20) << " bulk viscosity = " << mixture.get_bulk_viscosity() << std::endl;
  std::cout << std::setw(20) << " thermal conductivity = " << mixture.get_thermal_conductivity() << std::endl;
  arma::vec thd; // thermo-diffusion arma vector
  thd = mixture.get_thermodiffusion();
  std::cout << " thermal diffusion = " << mixture.get_thermodiffusion() << std::endl;

  for (auto moles : molecules) {
    std::cout << std::setw(20) << mol_ndens[0];
    std::cout << std::endl;
  }

  outf.close();
  return 0;
}
