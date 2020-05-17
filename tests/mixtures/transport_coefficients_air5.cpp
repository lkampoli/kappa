/*
   \file transport_coefficients_air5.cpp
   \brief Computation of transport coefficients for air5 mixture in STS approach.
*/

#include <iostream>
#include <fstream>
#include <iomanip> 

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

int main(int argc, char** argv) {
   
  std::cout << "Start computation of transport coefficients" << std::endl;
  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
  std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
  std::string particle_source    = m_source + "particles.yaml";
  std::string interaction_source = m_source + "interaction.yaml";
  std::string output_dir = GetCurrentWorkingDir();
  std::cout << "Current directory is: " << output_dir << std::endl;

  std::cout << "Loading particles data" << std::endl;
  
  kappa::Molecule N2("N2", true, false, particle_source);
  kappa::Molecule O2("O2", true, false, particle_source);
  kappa::Molecule NO("NO", true, false, particle_source);
  kappa::Atom N("N",                    particle_source);
  kappa::Atom O("O",                    particle_source);
  
  std::vector<kappa::Molecule> mol;
  std::vector<kappa::Atom> at;

  mol.push_back(N2);
  mol.push_back(O2);
  mol.push_back(NO);
  at.push_back(N);
  at.push_back(O);

  std::cout << "Finished loading particles data" << std::endl;

  //std::cout << "Molecule's name " << mol.name << std::endl;
  //std::cout << "Atom's name " << at.name << std::endl;
  //std::cout << "Molecule vibrational levels " << mol.num_vibr_levels[0] << std::endl;

  // air5 mixture
  kappa::Mixture mixture(mol, at, interaction_source, particle_source);

  // some check print
  std::cout << "particles: " << mixture.get_n_particles() << std::endl;
  std::cout << "names: " << mixture.get_names() << std::endl;

  // spectroscopic constant for diatomic molecules in ground electronic state
  int vibr_l = 0; // for N2 d0
  double Re = 1.097E-10; // N2
  double be = 2.25E-10;
  double omega_e = 235860; // m^-1 (N2)
  double mu = 0.028; // kg/mol (N2)
  double l_alpha = sqrt( 16.863/(omega_e * mu) );
  double beta = 2.6986E+10; // N2
  double d0 = Re + be + (9./2.)*beta*l_alpha*l_alpha*exp( 2*sqrt(beta*l_alpha)*(vibr_l - 1) );
   
  // set a range for temperature
  std::vector<double> T_vals;
  std::vector<double> X_vals;
  double tot_ndens;

  // vibrational levels
  int i;

  for (i=0; i<80; i++) { //500-40000K
    T_vals.push_back(500 + i * 500);
  }

  for (i=0; i<=20; i++) { //0.-1.
    X_vals.push_back(i/20.);
  }

  // arma vector for atom number density
  arma::vec atom_ndens(2);

  // vector of arma vector for molecular number density
  std::vector<arma::vec> mol_ndens;

  // assume an initial boltzmann ditribution 
  mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (kappa::K_CONST_K * T_vals[0]), N2));
  mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (kappa::K_CONST_K * T_vals[0]), O2));
  mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (kappa::K_CONST_K * T_vals[0]), NO));

  // output files
  std::ofstream outf, outf2, outf3, outf4;
  
  outf.open(output_dir + "/TCs_air5" + ".txt");

  outf << std::setw(20) << "Temperature [K]"; 
  outf << std::setw(20) << "Thermal conduction."; 
  outf << std::setw(20) << "Shear viscosity.";
  outf << std::setw(20) << "Bulk viscosity.";
  outf << std::endl;

  // transport coefficients
  double th_c; // thermal conductivity
  double sh_v; // shear viscosity
  double bk_v; // bulk viscosity
  double bk_o_sh; // ratio bulk/shear viscosity
  arma::vec th_d; // thermal diffusion
  arma::mat diff; // diffusion coeffs.

  double x_N2 = 1.0, x_O2 = 0.0, x_NO = 0.0, x_N = 0.0, x_O = 0.0;

   for (auto T : T_vals) {

     // total number density
     tot_ndens = 101325.0 / (kappa::K_CONST_K * T);

     mol_ndens[0] = mixture.Boltzmann_distribution(T, (1 - x_N2) * tot_ndens, N2);
     mol_ndens[1] = mixture.Boltzmann_distribution(T, (1 - x_O2) * tot_ndens, O2);
     mol_ndens[2] = mixture.Boltzmann_distribution(T, (1 - x_NO) * tot_ndens, NO);
     atom_ndens[0] = x_N * tot_ndens;
     atom_ndens[1] = x_O * tot_ndens;

     // computation of transport coefficients
     mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens);

     // retrieve thermal conductivity coefficients
     th_c = mixture.get_thermal_conductivity();

     // retrieve shear viscosity coefficients
     sh_v = mixture.get_shear_viscosity();

     // retrieve bulk viscosity coefficients
     bk_v = mixture.get_bulk_viscosity();

     outf << std::setw(20) << T; 
     outf << std::setw(20) << th_c; 
     outf << std::setw(20) << sh_v;
     outf << std::setw(20) << bk_v;
     outf << std::endl;
   }
     
   outf.close();
   outf2.close();
   outf3.close();
   outf4.close();

   return 0;
}
