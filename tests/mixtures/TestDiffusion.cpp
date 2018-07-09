/*! 
   \file TestDiffusion.cpp
   \brief Test for multi-component diffusion coefficients' computation in the StS approach.
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
   
  std::cout << "Start test: computation of thermal diffusion coefficients" << std::endl;
  
  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
  std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
  std::string particle_source    = m_source + "particles.yaml";
  std::string interaction_source = m_source + "interaction.yaml";
  std::string output_dir = GetCurrentWorkingDir();
  std::cout << "Current directory is: " << output_dir << std::endl;

  std::cout << "Loading particles data" << std::endl;

  // N2 molecule (non-rigid model)
  kappa::Molecule mol("N2", true, false, particle_source);
 
  // N atom
  kappa::Atom at("N", particle_source);

  std::cout << "Finished loading particles data" << std::endl;

  std::cout << "Molecule's name " << mol.name << std::endl;
  std::cout << "Atom's name " << at.name << std::endl;
  std::cout << "Molecule vibrational levels " << mol.num_vibr_levels[0] << std::endl;

  // N2/N binary mixture creation
  kappa::Mixture mixture(mol, at, interaction_source, particle_source);

  // some check print
  std::cout << "particles: " << mixture.get_n_particles() << std::endl;
  std::cout << "names: " << mixture.get_names() << std::endl;

  // set a range for temperature
  std::vector<double> T_vals = {10000.0};
  // std::vector<double> T_vals = {2500.0, 5000.0, 20000.0, 50000.0};
  // std::vector<double> T_vals = {15000.0};
  // std::vector<double> T_vals = {30000.0};
  // std::vector<double> T_vals = {1.1986210569919167e+04};

  // double p = 2.6298517223967945e+04; //101325.;
  double p = 101325.;

  // vibrational levels
  int i;

  // for (i = 0; i < 80; i++) { // assume a max num. of vibr. levels a priori
  //   T_vals.push_back(500 + i * 500);
  // }

  // set an arbitrary atom mass fraction (0.20)
  int x_atom_perc = 0.;
  double x_atom = x_atom_perc / 100.;

  // arma vector for atom number density
  arma::vec atom_ndens(1);

  // vector of arma vector for molecular number density
  std::vector<arma::vec> mol_ndens;

  double tot_ndens;
  tot_ndens =  p / (kappa::K_CONST_K * T_vals[0]);
  std::cout << " tot_ndens " << tot_ndens << std::endl;

  mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], (1. - x_atom) * tot_ndens, mol));

  // atom number density
  atom_ndens = x_atom * tot_ndens;

  for (int at=0; at<atom_ndens.size(); at++) {
    std::cout << " atom_ndens " << atom_ndens.at(at) << std::endl;
  }

  std::ofstream outf;
  outf.open(output_dir + "/diff_" + mol.name + "_" + at.name + "_xat" + std::to_string(x_atom_perc) + ".txt");

  // output files' header
  outf << std::setw(20) << "Temperature [K]"; 
  outf << std::setw(20) << "Diff. coeffs."; 
  outf << std::endl;

  std::cout << std::setw(20) << "Temperature [K]" << std::endl;

  mixture.compute_transport_coefficients(T_vals[0], mol_ndens, atom_ndens);
  mixture.get_diffusion();
    
 // simplified binary diffusion coefficients (for checking)
 // mixture.binary_diffusion(T_vals[0]);
 // mixture.binary_diffusion(T);

  return 0;
}
