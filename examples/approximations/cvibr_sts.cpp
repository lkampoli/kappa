/*
   \file cvibr_sts.cpp
   \brief compute specific vibrational heat at constant volume in the STS approach.
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

   
  double T;
  // create molecules and atoms
  Molecule MoleculeN2("N2", false, true, particle_source);
  Atom AtomN("N", particle_source);

  // select approximation level (STS)
  Approximation approx{};

  // select temperature
  T = 300.;
      
  // compute specific vibrational heat at constant volume
  std::cout << std::setw(20) << "T [K]";
  std::cout << std::setw(25) << "Vibr. heat capacity, cv_vib";
  std::cout << std::endl;
  while ( T <= 3000.0) {
    std::cout << std::setw(20) << T;
    std::cout << std::setw(25) << approx.c_vibr_approx(T, MoleculeN2);
    std::cout << std::endl;
    T += 100;
  }
}
