/*
   \file mixture-sts-basic.cpp
   \brief conversion molar to mass fraction
*/

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip> 

#include "kappa.hpp"

using namespace kappa;

int main(int argc, char** argv) {

  std::cout << "Start Test state-to-state mixture" << std::endl;
  std::string m_libPath = std::getenv("KAPPA_DATA_DIRECTORY");
  std::string particle_source    = m_libPath + "particles.yaml";
  std::string interaction_source = m_libPath + "interaction.yaml";

  // instantiates molecules and atoms kappa vectors
  std::vector<kappa::Molecule> molecules;
  std::vector<kappa::Atom> atoms;

  std::cout << "Loading particles data" << std::endl;

  // create molecules and atoms
  kappa::Molecule N2("N2", true, true, particle_source);
  kappa::Molecule O2("O2", true, true, particle_source);
  kappa::Molecule NO("NO", true, true, particle_source);
  kappa::Atom N("N",                   particle_source);
  kappa::Atom O("O",                   particle_source);

  // populate vectors with particles
  molecules.push_back(N2);
  molecules.push_back(O2);
  molecules.push_back(NO);
  atoms.push_back(N);
  atoms.push_back(O);
 
  std::cout << "Finished loading particles data" << std::endl;

  std::cout << "Size of molecules: " << molecules.size() << std::endl;
  std::cout << "Size of atoms: " << atoms.size() << std::endl;

  // create mixture
  kappa::Mixture mixture(molecules, atoms, interaction_source, particle_source);

  // evaluate particles interactions
  for (auto at: atoms) {
    for (auto mo: molecules) {
      std::cout << at.name << "+" << mo.name << ", interaction=" 
                << mixture.interaction(at, mo).particle1_name << "+" 
                << mixture.interaction(at, mo).particle2_name << std::endl;
    }
  }

  for (auto at: atoms) {
    for (auto at2: atoms) {
      std::cout << at.name << "+" << at2.name << ", interaction=" 
                << mixture.interaction(at, at2).particle1_name << " " 
      	        << mixture.interaction(at, at2).particle2_name << std::endl;
    }
  }

  for (auto mo: molecules) {
    for (auto mo2: molecules) {
      std::cout << mo.name << "+" << mo2.name << ", interaction=" 
                << mixture.interaction(mo, mo2).particle1_name << " " 
     	        << mixture.interaction(mo, mo2).particle2_name << std::endl;
    }
  }

  // compute boltzmann distribution
  std::vector<double> T_vals = { 500., 1000., 5000., 10000., 20000., 40000. };
  std::vector<arma::vec> mol_ndens;
  mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (K_CONST_K * T_vals[0]), N2));
  mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (K_CONST_K * T_vals[0]), O2));
  mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (K_CONST_K * T_vals[0]), NO));

  std::cout << std::setw(20) << " Temp. [K] ";
  std::cout << std::setw(20) << " N2 ";
  std::cout << std::setw(20) << " O2 ";
  std::cout << std::setw(20) << " NO ";
  std::cout << std::endl;
  for (auto T : T_vals) {
    mol_ndens[0] = mixture.Boltzmann_distribution(T, 101325.0 / (K_CONST_K * T), N2);
    mol_ndens[1] = mixture.Boltzmann_distribution(T, 101325.0 / (K_CONST_K * T), O2);
    mol_ndens[2] = mixture.Boltzmann_distribution(T, 101325.0 / (K_CONST_K * T), NO);
    std::cout << std::setw(20) << T;
    std::cout << std::setw(20) << mol_ndens[0][0] / (101325.0 / (K_CONST_K * T));
    std::cout << std::setw(20) << mol_ndens[1][0] / (101325.0 / (K_CONST_K * T));
    std::cout << std::setw(20) << mol_ndens[2][0] / (101325.0 / (K_CONST_K * T));
    std::cout << std::endl;
  }
  arma::vec atom_ndens(2);

  arma::vec x_molar_fractions = arma::zeros(molecules.size()+atoms.size());
  x_molar_fractions[0] = 0.1; // N2
  x_molar_fractions[1] = 0.2; // O2
  x_molar_fractions[2] = 0.3; // NO
  x_molar_fractions[3] = 0.2; // N
  x_molar_fractions[4] = 0.2; // O

  // molar to mass conversion
  arma::vec x_mass_fractions = arma::zeros(molecules.size()+atoms.size());
  x_mass_fractions = mixture.convert_molar_frac_to_mass(x_molar_fractions);
  std::cout << std::setw(20) << x_mass_fractions << "\n";
  std::cout << std::setw(20) << x_mass_fractions.size() << "\n";
  std::cout << std::setw(20) << x_mass_fractions.n_rows << "\n";
  std::cout << std::setw(20) << x_mass_fractions.n_cols << "\n";

  std::vector<double> mass_fractions;
  for ( int i = 0; i < x_mass_fractions.n_rows; ++i) {
     mass_fractions.push_back(x_mass_fractions(i));
     std::cout << std::setw(20) << mass_fractions[i] << "\n";
  };

   return 0;
}
