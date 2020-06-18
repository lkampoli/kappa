/*!
    \file: mixture-sts-bulk.cpp
    \brief Test for bulk viscosity in binary mixtures.
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

  std::cout << "Start test: computation of bulk viscosity" << std::endl;

  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
  std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
  std::string particle_source    = m_source + "particles.yaml";
  std::string interaction_source = m_source + "interaction.yaml";
  std::string output_dir = GetCurrentWorkingDir();
  std::cout << "Current directory is: " << output_dir << std::endl;

  std::vector<kappa::Molecule> molecules;
  std::vector<kappa::Atom> atoms;

  std::cout << "Loading particles data" << std::endl;

  //kappa::Molecule mol("N2", true, false, particle_source);
  kappa::Molecule mol("N2", true, true, particle_source);
  kappa::Atom at("N", particle_source);

  // kappa::Molecule mol("NO", false, true, particle_source);

  //kappa::Molecule mol("O2", true, false, particle_source);
  //kappa::Atom at("O", particle_source);

  // kappa::Atom N("N", particle_source);
  // atoms.push_back(N);
  // kappa::Atom O("O", particle_source);
  // atoms.push_back(O);
  molecules.push_back(mol);
  atoms.push_back(at);

  kappa::Mixture mixture(molecules, atoms, interaction_source, particle_source);
  std::cout << "Mixture created" << std::endl;

  //std::vector<double> T_vals = { 500., 1000., 5000., 10000., 20000., 40000. };
  //std::vector<double> T_vals = { 500. };
  std::vector<double> T_vals;
  // double x_atom = 0.0;
  double x_N2 = 0.90;
  double x_O2 = 0.90;
  double x_NO = 0.0;
  double x_N = 0.10;
  double x_O = 0.10;
  double pressure = 101325.;
  //double n_tot = pressure / (K_CONST_K * 500);

   for (int i=0; i<195; i++) { //500-40000K
     T_vals.push_back(500 + i * 100);
   }

  //for (int i=0; i<50; i++) {
  //  T_vals.push_back(100 + i * 50);
  //}

  std::vector<arma::vec> mol_ndens;
  mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], x_N2 * 101325.0 / (K_CONST_K * T_vals[0]), mol)); // N2
  //mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], x_O2 * 101325.0 / (K_CONST_K * T_vals[0]), mol)); // O2

  std::ofstream outf;
  outf.open(output_dir + "/TRANSPORT_COEFFICIENTS/bulk_viscosity/" + mol.name + "_xat_" + std::to_string(x_N2) + ".txt");
  //outf.open(output_dir + "/TRANSPORT_COEFFICIENTS/bulk_viscosity/" + mol.name + "_xat_" + std::to_string(x_O2) + ".txt");
  outf << std::setw(20) << "Temperature [K]";
  outf << std::setw(20) << "Eta";
  outf << std::endl;

  arma::vec atom_ndens(2);
  atom_ndens[0] = 0.;
  //atom_ndens[1] = 0.;

  double tot_ndens;

//  std::vector<models_omega> omega_integral_models = {models_omega::model_omega_rs, models_omega::model_omega_vss,
//                                                     models_omega::model_omega_bornmayer, models_omega::model_omega_lennardjones,
//                                                     models_omega::model_omega_esa};

  for (auto T : T_vals) {
    tot_ndens =  101325.0 / (K_CONST_K * T);
    mol_ndens[0] = mixture.Boltzmann_distribution(T, x_N2 * 101325.0 / (K_CONST_K * T), mol);
    atom_ndens[0] = x_N * 101325.0 / (K_CONST_K * T);
    //mol_ndens[0] = mixture.Boltzmann_distribution(T, x_O2 * 101325.0 / (K_CONST_K * T), mol);
    //atom_ndens[0] = x_O * 101325.0 / (K_CONST_K * T);

    // Attention: looks like the perturbation is non-zero is not defined! Better to pass it explicitely
    mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens, 0, models_omega::model_omega_rs, 0.0);

    std::cout << T << "     " << mixture.get_bulk_viscosity() << std::endl;

    outf << std::setw(20) << T;
    outf << std::setw(20) << mixture.get_bulk_viscosity();
    outf << std::endl;
  }

  //std::cout << std::setw(20) << T_vals[0] << std::endl;
  //std::cout << std::setw(20) << mixture.get_bulk_viscosity() << std::endl;
  //for (auto moles : molecules) {
  //  std::cout << std::setw(20) << mol_ndens[0];
  //  std::cout << std::endl;
  //}

  //std::cout << "MOLECULAR NUMBER DENSITIES"  << std::endl;
  //std::cout << mol_ndens[0]  << std::endl;
  //std::cout << "ATOMIC NUMBER DENSITIES"  << std::endl;
  //std::cout << atom_ndens[0]  << std::endl;
  //std::cout << "MOLECULAR MOLAR FRACTIONS"  << std::endl;
  //std::cout << mol_ndens[0] / 101325.0 * (K_CONST_K * T_vals[0]) << std::endl;
  //std::cout << "ATOMIC MOLAR FRACTIONS"  << std::endl;
  //std::cout << atom_ndens[0] / 101325.0 * (K_CONST_K * T_vals[0])  << std::endl;
  //std::cout << "MOLECULAR MASS FRACTIONS"  << std::endl;
  //std::cout << mol_ndens[0] / 101325.0 * (K_CONST_K * T_vals[0])  << std::endl;
  //std::cout << "ATOMIC MASS FRACTIONS"  << std::endl;
  //std::cout << atom_ndens[0] / 101325.0 * (K_CONST_K * T_vals[0])  << std::endl;

  //double density = mixture.compute_density(mol_ndens, atom_ndens, 0);
  //std::cout << "density" << density << std::endl;

  // molar to mass conversion
  //double rho = 0.0;
  //for (int i=0; i<48; i++) {
  //  rho += mol_ndens[0].at(i) * mol.mass;
  //}
  //rho += atom_ndens[0] * at.mass;
  //std::cout << "rho" << rho << std::endl;

  //arma::vec MassFractions = arma::zeros(49);
  //for (int i=0; i<48; i++) {
  //  MassFractions(i) = mol_ndens[0].at(i) * mol.mass / rho;
  //}
  //MassFractions(48) = atom_ndens[0] * at.mass / rho;
  //std::cout << MassFractions << std::endl;

//  arma::vec x_molar_fractions = arma::zeros(49);
//  arma::vec x_mass_fractions = arma::zeros(49);
//
//  x_molar_fractions =  mol_ndens[0] / 101325.0 * (K_CONST_K * T_vals[0]);
//  x_mass_fractions = mixture.convert_molar_frac_to_mass(x_molar_fractions);
//  std::cout << std::setw(20) << x_mass_fractions << "\n";
//
//  std::vector<double> mass_fractions;
//  for ( int i = 0; i < x_mass_fractions.n_rows; ++i) {
//     mass_fractions.push_back(x_mass_fractions(i));
//     std::cout << std::setw(20) << mass_fractions[i] << "\n";
//  };


  //arma::vec thd; // thermo-diffusion arma vector
  //thd = mixture.get_thermodiffusion();
  //std::cout << "mixture.get_thermodiffusion()" << std::endl;
  //std::cout << thd << std::endl;

  outf.close();
  return 0;
}
