/*
   \fileionized_mixture_creation.cpp
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

  std::cout << "Loading particles data" << std::endl;

  kappa::Molecule mol("N2", true, true, 	particle_source);
  kappa::Molecule mol_ion("N2+", true, true, 	particle_source);
  kappa::Atom at("N", 				particle_source);
  kappa::Atom at_o("O", 			particle_source);
  kappa::Atom at_ion("N+", 			particle_source);

  std::cout << "particle's name " << at.name << "\n";
  std::cout << "particle's type " << at.particleType << "\n";
  std::cout << "particle's atomic weigth " << at.atomic_weight << "\n";

  for(int i=0; i<at.stoichiometry.size(); i++) {
    std::cout << at.stoichiometry[i].first << " " 
              << at.stoichiometry[i].second << std::endl;
    std::cout << "molar weight: " 
              << at.stoichiometry[i].second * at.atomic_weight << std::endl;
  }

  // create a vector pair of elements
  std::vector<std::pair<std::string,double>> elements;
  elements.push_back(std::make_pair(at.name,at.atomic_weight));
  elements.push_back(std::make_pair(at_o.name,at_o.atomic_weight));
  for(int i=0; i<elements.size(); i++) 
    std::cout << "vector pair: " << elements[i].first << ", " << elements[i].second << std::endl;
  }

  std::cout << "Finished loading particles data" << std::endl;

  kappa::Mixture mixture1(mol_ion, interaction_source, particle_source);
  std::cout << "Finished loading N2+,e- Mixture" << std::endl;

  kappa::Mixture mixture2(at_ion, interaction_source, particle_source);
  std::cout << "Finished loading N+,e- Mixture" << std::endl;

  kappa::Mixture mixture3(mol_ion, at, interaction_source, particle_source);
  std::cout << "Finished loading N2+,N,e- Mixture" << std::endl;

  kappa::Mixture mixture4("N2+, N2", interaction_source, particle_source);
  std::cout << "Finished loading N2+,N2,e- Mixture" << std::endl;

  kappa::Mixture mixture5("N2, N, N+", interaction_source, particle_source);
  std::cout << "Finished loading N2,N,N+,e- Mixture" << std::endl;

  std::cout << mixture1.get_names() << std::endl;
  std::cout << mixture2.get_names() << std::endl;
  std::cout << mixture3.get_names() << std::endl;
  std::cout << mixture4.get_names() << std::endl;
  std::cout << mixture5.get_names() << std::endl << std::endl;

  kappa::Interaction inter(mol, mol, interaction_source);

  inter = mixture5.interaction(mol, mol); // N2 + N2
  std::cout << inter.particle1_name << " " << inter.particle2_name << std::endl;
  inter = mixture5.interaction(mol, at); // N2 + N
  std::cout << inter.particle1_name << " " << inter.particle2_name << std::endl;
  inter = mixture5.interaction(mol, at_ion); // N2 + N+
  std::cout << inter.particle1_name << " " << inter.particle2_name << std::endl << std::endl;
  inter = mixture5.interaction(at, at_ion); // N + N+
  std::cout << inter.particle1_name << " " << inter.particle2_name << std::endl << std::endl;

  return 0;
}
