/*
   \file mixture-sts-shear-bulk-thermal.cpp
*/

#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip> 
#include <string>

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

  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
  std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
  std::string particle_source    = m_source + "particles.yaml";
  std::string interaction_source = m_source + "interaction.yaml";
  std::string output_dir = GetCurrentWorkingDir();
  std::cout << "Current directory is: " << output_dir << std::endl;

  std::cout << "Start test: computation of shear viscosity" << std::endl;

  std::vector<kappa::Molecule> molecules, moleculesRR;
  std::vector<kappa::Atom> atoms;

  std::cout << "Loading particles data" << std::endl;

  kappa::Molecule N2("N2", true, false, particle_source);
  kappa::Molecule N2_RR("N2", true, true, particle_source);
  // kappa::Molecule O2("O2", true, true, particle_source);
  // kappa::Molecule NO("NO", true, true, particle_source);
  kappa::Atom N("N", particle_source);
  // kappa::Atom O("O", particle_source);
  
  molecules.push_back(N2); 
  moleculesRR.push_back(N2_RR); 
  // molecules.push_back(O2);
  // molecules.push_back(NO);
  atoms.push_back(N);
  // atoms.push_back(O);

  std::cout << "Finished loading particles data" << std::endl;
  
  kappa::Mixture mixture(molecules, atoms, interaction_source, particle_source);
  kappa::Mixture mixture_RR(moleculesRR, atoms, interaction_source, particle_source);

  std::vector<double> T_vals = { 500., 1000., 5000., 10000., 15000., 20000., 25000., 30000., 35000., 40000. };
  // double x_N2 = 0.2, x_O2 = 0.2, x_NO = 0.2, x_N = 0.2, x_O = 0.2;
  double x_N2 = 1.0, x_O2 = 0.0, x_NO = 0.0, x_N = 0.0, x_O = 0.0;

  std::vector<arma::vec> mol_ndens;
  arma::vec atom_ndens(2);

  for (auto moles : molecules) {
    mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (K_CONST_K * T_vals[0]), moles));
    std::cout << std::setw(20) << mol_ndens[0];
    std::cout << std::endl;
  }

  // prepare the output files
  std::ofstream outf,outf2;
  outf.open(output_dir + "/coeffs_" + N2.name + "_" + N.name + "_xat" + std::to_string(x_N) + ".txt");
  outf2.open(output_dir + "/coeffs_ratio_" + N2.name + "_" + N.name + "_xat" + std::to_string(x_N) + ".txt");
  outf << std::setw(20) << "Temperature [K]";
  outf << std::setw(20) << "shear visc.";
  outf << std::setw(20) << "bulk visc.";
  outf << std::setw(20) << "thermal-diffusion";
  outf << std::endl;

  outf2 << std::setw(20) << "Temperature [K]";
  outf2 << std::setw(20) << "shear visc.";
  outf2 << std::setw(20) << "bulk visc.";
  outf2 << std::setw(20) << "thermal-diffusion";
  outf2 << std::endl;

  // mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (K_CONST_K * T_vals[0]), N2));
  // mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (K_CONST_K * T_vals[0]), O2));
  // mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (K_CONST_K * T_vals[0]), NO));
    
  double tot_ndens;

  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

  // Write a header line for the table
  std::cout << std::setw(25) << "temperature";
  std::cout << std::setw(25) << "shear viscosity";
  std::cout << std::setw(25) << "bulk viscosity";
  std::cout << std::setw(25) << "thermal conductivity";
  std::cout << std::endl;

  for (auto T : T_vals) {
    tot_ndens =  101325.0 / (K_CONST_K * T);
    mol_ndens[0] = mixture.Boltzmann_distribution(T, (1.-x_N) * tot_ndens, N2);
    mol_ndens[0] = mixture_RR.Boltzmann_distribution(T, (1.-x_N) * tot_ndens, N2_RR);
    // mol_ndens[1] = mixture.Boltzmann_distribution(T, x_O2 * tot_ndens, O2);
    // mol_ndens[2] = mixture.Boltzmann_distribution(T, x_NO * tot_ndens, NO);
    atom_ndens[0] = x_N * tot_ndens;
    // atom_ndens[1] = x_O * tot_ndens;

    mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens);
    mixture_RR.compute_transport_coefficients(T, mol_ndens, atom_ndens);
  
    // std::cout << std::setw(25) << T;
    // std::cout << std::setw(25) << mixture.get_shear_viscosity();
    // std::cout << std::setw(25) << mixture.get_bulk_viscosity();
    // std::cout << std::setw(25) << mixture.get_thermal_conductivity();
    // std::cout << std::endl;

     outf << std::setw(20) << T;
     outf << std::setw(20) << mixture.get_shear_viscosity();
     outf << std::setw(20) << mixture.get_bulk_viscosity();
     outf << std::setw(20) << mixture.get_thermal_conductivity();
     outf << std::endl;

     outf2 << std::setw(20) << T;
     outf2 << std::setw(20) << mixture.get_shear_viscosity()/mixture_RR.get_shear_viscosity();
     outf2 << std::setw(20) << mixture.get_bulk_viscosity()/mixture_RR.get_bulk_viscosity();
     outf2 << std::setw(20) << mixture.get_thermal_conductivity()/mixture_RR.get_thermal_conductivity();
     outf2 << std::endl;
  }

  // timing
  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
  std::cout << " elapsed time: " << duration << " ms "<< std::endl;

  return 0;
}
