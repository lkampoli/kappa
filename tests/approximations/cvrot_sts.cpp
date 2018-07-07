/*
   \file cvrot_sts.cpp
   \brief test for rotational specific heat at constant volume for STS approach
*/

#include <fstream>
#include <iomanip> 
#include <cstdlib>
#include <iostream>

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
using namespace arma; 

int main(){

  // Retrieve the kappa data path from environment variable defined in ~/.bashrc
  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
  std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
  std::string particle_source    = m_source + "particles.yaml";
  std::string interaction_source = m_source + "interaction.yaml";
  std::string output_dir = GetCurrentWorkingDir();
  std::cout << "Current directory is: " << output_dir << std::endl;

  double T = 0.;
  double p = 101325.;
  double tot_ndens = 0.;
  double cv_rot = 0.;
  double rho = 0.;
  double n = 0.;
  double n_at = 0.;
  double n_mol = 0.;
  double pressure = 0.;

  Molecule MoleculeN2("N2", false, true, particle_source);
  Atom AtomN("N", particle_source);

  std::vector<kappa::Molecule> molecules;
  std::vector<kappa::Atom> atoms;

  molecules.push_back(MoleculeN2);
  atoms.push_back(AtomN);
  
  // create mixture
  kappa::Mixture mixture(molecules, atoms, interaction_source, particle_source);

  // select approximation level
  Approximation approx{};

  std::cout << std::setw(20) << "T [K]";
  std::cout << std::setw(25) << "Rot. specific heat capacity, cv_rot";
  std::cout << std::endl;
  
  std::vector<arma::vec> mol_ndens;
  arma::vec atom_ndens(1);
  int x_atom_perc = 50.; // assume a certain atom percentage
  double x_atom = x_atom_perc / 100.;
  atom_ndens[0] = x_atom * tot_ndens;

  // for (int at=0; at<atom_ndens.size(); ++at) {
  //   std::cout << " at_ndens " << atom_ndens.at(at) << std::endl;
  // }

  // compute boltzmann distribution
  mol_ndens.push_back(mixture.Boltzmann_distribution(T, (1.-x_atom) * tot_ndens, MoleculeN2));
  while (T <= 3000.0) {

    // mixture number density
    tot_ndens = p / (kappa::K_CONST_K * T);
 
    mol_ndens[0] = mixture.Boltzmann_distribution(T, (1.-x_atom) * tot_ndens, MoleculeN2);

    // for (auto i=mol_ndens.begin(); i!=mol_ndens.end(); ++i) {
    //   std::cout << " mol_ndens " << (*i) << std::endl;
    // }

    // mixture mass density
    rho = mixture.compute_density(mol_ndens, atom_ndens); 
    // std::cout << " rho " << rho << std::endl;

    // mixture pressure
    pressure = mixture.compute_pressure(T, mol_ndens, atom_ndens); 
    // std::cout << " pressure " << pressure << std::endl;

    n = mixture.compute_n(mol_ndens, atom_ndens); 
    n_at = mixture.compute_n(atom_ndens); 
    n_mol = mixture.compute_n(mol_ndens); 
 
    // cross-check, eq. 1.23
    cv_rot = kappa::K_CONST_K / MoleculeN2.mass;

    int c_rot = 0;
    for (int i=0; i<MoleculeN2.num_vibr_levels[0]; i++) {
      c_rot += approx.c_rot(T, MoleculeN2, i, 0);
    }

    std::cout << std::setw(20) << T;
    std::cout << std::setw(25) << approx.c_rot(T, MoleculeN2) 
              << std::setw(20) << cv_rot;
    std::cout << std::endl;

    // conservation test
    double test=0.;
    for (int i=0; i<mol_ndens[0].size(); i++) {
      test += mol_ndens[0].at(i);
    }
    // std::cout << " conservation tests (expected value = 1) " << (test + atom_ndens[0]) / n << std::endl;

    T += 10;
  }
}
