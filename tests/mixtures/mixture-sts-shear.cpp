/*! \file: mixture-sts-shear.cpp 
 *  \brief Test for shear viscosity in binary mixtures
 */

#include <iostream>
#include "kappa.hpp"
#include <fstream>
#include <string>
#include <iomanip> // std::setw
#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

using namespace kappa;

 std::string GetCurrentWorkingDir( void ) 
 {
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
    //kappa::Molecule mol("N2", false, true, particle_source);
    //kappa::Atom at("N", particle_source);
    //kappa::Molecule mol("O2", false, true, particle_source);
    //kappa::Atom at("O", particle_source);
    kappa::Atom N("N", particle_source);
    atoms.push_back(N);
    kappa::Atom O("O", particle_source);
    atoms.push_back(O);
    molecules.push_back(mol);
    //atoms.push_back(at);

// George stuff
//  kappa::Molecule mol("N2", true, true, particle_source);
//  kappa::Mixture mix("N2", interaction_source, particle_source);
//  std::vector<double> T_arr = {500., 1000., 2000.};
//  arma::vec n = arma::zeros(1);
//  std::vector<arma::vec> nmol;
//  nmol.push_back(mix.Boltzmann_distribution(500., 1e23, mol));
//  n[0] = 1e23;
//  std::cout << mix.get_names() << std::endl;
//  for (auto T: T_arr) {
//    nmol[0] = mix.Boltzmann_distribution(T, 1e23, mol);
//    mix.compute_transport_coefficients(T, nmol);
////  std::cout << T << " " << mix.get_shear_viscosity() << std::endl;
//    std::cout << T << " " << nmol[0] << std::endl;
//  }

//    kappa::Molecule N2("N2", true, true, particle_source);
//    kappa::Molecule O2("O2", true, true, particle_source);
//    //kappa::Molecule NO("NO", true, true, particle_source);
//    kappa::Atom N("N",                   particle_source);
//    kappa::Atom O("O",                   particle_source);
//    molecules.push_back(N2);
//    molecules.push_back(O2);
//    //molecules.push_back(NO);
//    atoms.push_back(N);
//    atoms.push_back(O);
//
//    std::cout << "Finished loading particles data" << std::endl;
//
//    int vibr_l = 0.; // for N2 d0
//    double Re = 1.097E-10; // N2
//    double be = 2.25E-10;
//    double omega_e = 235860.; // m^-1 (N2)
//    double mu = 0.028; // kg/mol (N2)
//    double l_alpha = sqrt(16.863/(omega_e * mu));
//    double beta = 2.6986E+10; // N2
//    double d0 = Re + be + (9./2.)*beta*l_alpha*l_alpha*exp( 2.*sqrt(beta*l_alpha)*(vibr_l - 1.) );
//
//
//
    kappa::Mixture mixture(molecules, atoms, interaction_source, particle_source);
    std::cout << "Mixture created" << std::endl;

    //std::vector<double> T_vals = { 500., 1000., 5000., 10000., 20000., 40000. };
    std::vector<double> T_vals;
    //double x_atom = 0.0;
    double x_N2 = 1.0;
    double x_NO = 1.0;
    double x_N = 0.0;
    double x_O = 0.0;
    double pressure = 101325.;
    //double n_tot = pressure / (K_CONST_K * T_vals[0]);
    double n_tot = pressure / (K_CONST_K * 500);

    //for (int i = 0; i < 80; i++) { //500-40000K
    // T_vals.push_back(500 + i * 500);
    //}
    for (int i = 0; i < 50; i++) { 
     T_vals.push_back(100 + i * 50);
    }

    std::vector<arma::vec> mol_ndens;
    //mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], n_tot, mol));
    //mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], n_tot, mol1));
    //mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], x_N2 * 101325.0 / (K_CONST_K * T_vals[0]), N2));
    // mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], (1. - x_N2) * 101325.0 / (K_CONST_K * T_vals[0]), O2));
    //mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (K_CONST_K * T_vals[0]), NO));
    mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], x_NO * 101325.0 / (K_CONST_K * T_vals[0]), mol)); // NO

    std::ofstream outf;
    //outf.open(output_dir + "/TRANSPORT_COEFFICIENTS/shear_viscosity/" + mol.name + "_" + at.name + "_xat_" + std::to_string(x_N2) + "_pure_new_diameter.txt");
    //outf.open(output_dir + "/TRANSPORT_COEFFICIENTS/shear_viscosity/" + N2.name + "_" + O2.name + "_xat_" + std::to_string(x_N2) + ".txt");
    outf.open(output_dir + "/TRANSPORT_COEFFICIENTS/shear_viscosity/" + mol.name + "_xat_" + std::to_string(x_NO) + "new_diameter.txt");

    outf << std::setw(20) << "Temperature [K]";
    outf << std::setw(20) << "Eta";
    //outf << std::setw(20) << "H_N2N2";
    //outf << std::setw(20) << "H_N2N";
    //outf << std::setw(20) << "H_NN";
    //outf << std::setw(20) << "omega_11_N2N2";
    //outf << std::setw(20) << "omega_22_N2N2";
    //outf << std::setw(20) << "omega_11_NN";
    //outf << std::setw(20) << "omega_22_NN";
    outf << std::endl;

    arma::vec atom_ndens(2);
    atom_ndens[0] = 0.;
    atom_ndens[1] = 0.;

    double tot_ndens;

    for (auto T : T_vals) {
     tot_ndens =  101325.0 / (K_CONST_K * T);
     //mol_ndens[0] = mixture.Boltzmann_distribution(T, x_N2 * 101325.0 / (K_CONST_K * T), N2);
     //mol_ndens[1] = mixture.Boltzmann_distribution(T, (1. - x_N2) * 101325.0 / (K_CONST_K * T), O2);
      mol_ndens[0] = mixture.Boltzmann_distribution(T, x_NO * 101325.0 / (K_CONST_K * T), mol);
     //mol_ndens[2] = mixture.Boltzmann_distribution(T, 101325.0 / (K_CONST_K * T), NO);
     //mol_ndens[0] = mixture.Boltzmann_distribution(T, tot_ndens, mol);
     //atom_ndens[0] = x_N2 * tot_ndens;
     //std::cout << atom_ndens[0] << std::endl;

     //double eta0 = (5./(16.*d0*d0))*sqrt(mol.mass * K_CONST_K * T/K_CONST_PI);

//   arma::vec molarmass;
//   std::cout << mixture.get_n_vibr_levels() << std::endl;
//   molarmass=mixture.convert_molar_frac_to_mass(atom_ndens);
//   mixture.convert_molar_frac_to_mass(atom_ndens);
//   molarmass=mixture.convert_mass_frac_to_molar(atom_ndens);
//   std::cout << "molar mass: " << molarmass << std::endl;

     //mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens);
     mixture.compute_transport_coefficients(T, mol_ndens);

     outf << std::setw(20) << T;
     outf << std::setw(20) << mixture.get_shear_viscosity();
     //outf << std::setw(20) << mixture.get_shear_viscosity()/eta0;
     outf << std::endl;
    }

    outf.close();
    return 0;
}
