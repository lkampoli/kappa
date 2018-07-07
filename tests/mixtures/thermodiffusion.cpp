/* 
 * \file thermodiffusion.cpp
 * \brief Test for thermal-diffusion coefficients' computation in the StS approach.
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

  // N2 molecule (rigid and non rigid rotator model)
  kappa::Molecule mol("N2", true, true,         particle_source);
  kappa::Molecule mol_nonrig("N2", true, false, particle_source);
 
  // N atom
  kappa::Atom at("N", particle_source);

  std::cout << "Finished loading particles data" << std::endl;

  std::cout << "Molecule's name " << mol.name << std::endl;
  std::cout << "Atom's name " << at.name << std::endl;
  std::cout << "Molecule vibrational levels " << mol.num_vibr_levels[0] << std::endl;

  // N2/N binary mixture creation
  kappa::Mixture mixture(mol, at,               interaction_source, particle_source);
  kappa::Mixture mixture_nonrig(mol_nonrig, at, interaction_source, particle_source);

  // some check print
  std::cout << "particles: " << mixture.get_n_particles() << std::endl;
  std::cout << "names: " << mixture.get_names() << std::endl;

  // set a range for temperature
  std::vector<double> T_vals = {2500.0, 5000.0, 20000.0, 50000.0};

  // vibrational levels
  int i;

  // spectroscopic constant for diatomic molecules in ground electronic state 
  int vibr_l = 0; // for N2 d0
  double Re = 1.097E-10; // N2
  double be = 2.25E-10;
  double omega_e = 235860; // m^-1 (N2)
  double mu = 0.028; // kg/mol (N2)
  double l_alpha = sqrt( 16.863/(omega_e * mu) );
  double beta = 2.6986E+10; // N2
  double d0 = Re + be + (9./2.)*beta*l_alpha*l_alpha*exp( 2*sqrt(beta*l_alpha)*(vibr_l - 1) );

  for (i=0; i<159; i++) { // assume a max num. of vibr. levels a priori
    T_vals.push_back(500 + i * 500);
  }

  // set an arbitrary atom mass fraction
  int x_atom_perc = 95.0;
  double x_atom = x_atom_perc / 100.;

  // arma vector for atom number density
  arma::vec atom_ndens(1);

  // vector of arma vector for molecular number density
  std::vector<arma::vec> mol_ndens;

  // assume an initial boltzmann ditribution 
  mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (kappa::K_CONST_K * T_vals[0]), mol));

  // output files
  std::ofstream outf, outf2;
  outf.open(output_dir + "/thd_" + mol.name + "_" + at.name + "_xat" + std::to_string(x_atom_perc) + ".txt");
  outf2.open(output_dir + "/thd_sts_" + mol.name + "_" + at.name + "_xat" + std::to_string(x_atom_perc) + ".txt");

  // output files' header
  outf << std::setw(20) << "Temperature [K]"; 
  outf << std::setw(20) << "Lambda"; 
  outf << std::endl;

  outf2 << std::setw(20) << "Temperature [K]"; 
  outf2 << std::setw(20) << "Lambda"; 
  outf2 << std::endl;
  
  double tot_ndens;
  arma::vec thd; // thermo-diffusion arma vector

  std::cout << std::setw(20) << "Temperature [K]" << std::endl;

  // main loop on temperatures
  for (auto T : T_vals) {

    tot_ndens =  101325.0 / (kappa::K_CONST_K * T);
    mol_ndens[0] = mixture.Boltzmann_distribution(T, (1 - x_atom) * tot_ndens, mol);
    atom_ndens[0] = x_atom * tot_ndens;

    std::cout << std::setw(20) << T << std::endl;

    // computation of transport coefficients
    mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens);
    mixture_nonrig.compute_transport_coefficients(T, mol_ndens, atom_ndens);

    double D0 = (3./(8.* tot_ndens * d0 * d0) ) * sqrt( (kappa::K_CONST_K * T)/(kappa::K_CONST_PI*mol.mass) );

    // retrieve thermal-diffusion coefficients (Rigid Rotator)
    thd = mixture.get_thermodiffusion();


    for (i=0; i<thd.n_elem; i++) {
      outf << std::setw(20) << T; 
      outf << std::setw(25) << std::setprecision(18) << thd[i]/D0; 
      outf << std::endl;
    }

    // retrieve thermal-diffusion coefficients (Non-Rigid Rotator)
    thd = mixture_nonrig.get_thermodiffusion();
    for (i=0; i<thd.n_elem; i++) {
      outf2 << std::setw(20) << T; 
      outf2 << std::setw(25) << std::setprecision(18) << thd[i]/D0; 
      outf2 << std::endl;
    }
  }

  outf.close();
  outf2.close();

  return 0;
}
