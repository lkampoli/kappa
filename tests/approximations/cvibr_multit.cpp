/*
   \file cvibr_multit.cpp
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

  double T, T1;

  // name, anharmonic_spectrum, rigid_rotator, filename       
  Molecule Molecule1_H(  "N2", false, true,  particle_source);
  Molecule Molecule1_AnH("N2", true,  true,  particle_source);
  Molecule Molecule2_H(  "O2", false, true,  particle_source);
  Molecule Molecule2_AnH("O2", true,  true,  particle_source);

  // Create an instance of ApproximationMultiT model
  ApproximationMultiT approx{};

  // ... just to check the corret name 
  cout << "Molecule1_H" << Molecule1_H.name << endl;
      
  // Initialize temperatures
  T1 = 2000.0; T = 100.0;
      
  cout << "Calculation of the vibrational heat capacity at T1 = " << T1 << " and different values of T [100, 3000]" << endl;
  cout << " - c_vibr_T(T, T1, AnHarmonic)" << endl;

  std::cout << std::setw(20) << "T1 [K]";
  std::cout << std::setw(20) << "T [K]";
  std::cout << std::setw(25) << "Vibr. heat capacity";
  std::cout << std::endl;
  while ( T <= 3000.0) {
    std::cout << std::setw(20) << T1;
    std::cout << std::setw(20) << T;
    std::cout << std::setw(25) << -approx.c_vibr_T(T, T1, Molecule1_AnH);
    std::cout << std::endl;
    T += 100;
  }
  
  // Reset temperatures
  T1 = 2000.0; T = 100.0;

  std::cout << "c_vibr_T1(T, T1, AnHarmonic)" << endl;
  std::cout << std::setw(20) << "T1 [K]";
  std::cout << std::setw(20) << "T [K]";
  std::cout << std::setw(25) << "Vibr. heat capacity";
  std::cout << std::endl;
  while (T <= 3000.0) {
    std::cout << std::setw(20) << T1;
    std::cout << std::setw(20) << T;
    std::cout << std::setw(25) << approx.c_vibr_T1(T, T1, Molecule1_AnH);
    std::cout << std::endl;
    T += 100;
  }

  T1 = 2000.0; T = 100.0;
  cout << "c_vibr_T1(T, T1, Harmonic)" << endl;
  while (T <= 3000.0) {
    cout  << approx.c_vibr_T1(T, T1, Molecule1_H) << endl;
    T += 100;
  }

  std::cout << "My chart" << endl;
  std::cout << "Molecule name = " << Molecule1_AnH.name << endl;
  T = 1000; T1 = 1000;
  cout << "Calculation of the vibrational heat capacity at T1 =" << T1 << " and different values of T [1000, 20000], anharmonic spectrum." << endl;
  cout << "- c_vibr_T(Anharmonic, N2)  " << endl;
  while (T <= 20000) {
    std::cout  << - approx.c_vibr_T(T, T1, Molecule1_AnH) << endl;
    T += 500;
  }

  std::cout << endl;
  T = 1000; T1 = 1000;
  cout << "Calculation of the vibrational heat capacity at T1 =" << T1 << " and different values of T [1000, 20000], harmonic spectrum." << endl;
  cout << " c_vibr_T(Harmonic, N2)  " << endl;
  while (T <= 20000) {
    std::cout  << approx.c_vibr_T(T, T1, Molecule1_H) << endl;
    T += 500;
  }

  std::cout << endl;
  std::cout << "***************************************************************" << endl;
  T = 1000; T1 = 1000;
  cout << "Calculation of the vibrational heat capacity at T1 =" << T1 << " and different values of T [1000, 20000], anharmonic spectrum." << endl;
  cout << " c_vibr_T1(Anharmonic, N2)  " << endl;
  while (T <= 20000) {
    std::cout  << approx.c_vibr_T1(T, T1, Molecule1_AnH) << endl;
    T += 500;
  }

  std::cout << endl;
  T = 1000; T1 = 1000;
  cout << "Calculation of the vibrational heat capacity at T1 =" << T1 << " and different values of T [1000, 20000], harmonic spectrum." << endl;
  cout << " c_vibr_T1(Harmonic, N2)  " << endl;
  while (T <= 20000) {
    std::cout << approx.c_vibr_T1(T, T1, Molecule1_H) << endl;
    T += 500;
  }
  std::cout << endl;
  std::cout << endl;
  std::cout << "***************************************************************" << endl;
  std::cout << "***************************************************************" << endl;

  T = 5000; T1 = 1000;
  cout << "Calculation of the vibrational heat capacity at T =" << T << " and different values of T1 [1000, 20000], anharmonic spectrum." << endl;
  cout << "- c_vibr_T(Anharmonic, N2)  " << endl;
  while (T1 <= 20000) {
    std::cout << -approx.c_vibr_T(T, T1, Molecule1_AnH) << endl;
    T1 += 500;
  }
  std::cout << endl;
  T = 5000; T1 = 1000;
  cout << "Calculation of the vibrational heat capacity at T =" << T << " and different values of T1 [1000, 20000], harmonic spectrum." << endl;
  cout << " c_vibr_T(Harmonic, N2)  " << endl;
  while (T1 <= 20000) {
    std::cout  << approx.c_vibr_T(T, T1, Molecule1_H) << endl;
    T1 += 500;
  }
  std::cout << endl;
  std::cout << "***************************************************************" << endl;
  T = 5000; T1 = 1000;
  cout << "Calculation of the vibrational heat capacity at T =" << T << " and different values of T1 [1000, 20000], anharmonic spectrum." << endl;
  cout << " c_vibr_T1(Anharmonic, N2)  " << endl;
  while (T1 <= 20000) {
    std::cout  << approx.c_vibr_T1(T, T1, Molecule1_AnH) << endl;
    T1 += 500;
  }
  std::cout << endl;
  T = 5000; T1 = 1000;
  cout << "Calculation of the vibrational heat capacity at T =" << T << " and different values of T1 [1000, 20000], harmonic spectrum." << endl;
  cout << " c_vibr_T1(Harmonic, N2)  " << endl;
  while (T1 <= 20000) {
    std::cout  << approx.c_vibr_T1(T, T1, Molecule1_H) << endl;
    T1 += 500;
  }
  std::cout << endl;
  std::cout << endl;
  std::cout << "***************************************************************" << endl;
  std::cout << "***************************************************************" << endl;
  std::cout << endl;
  std::cout << endl;
  std::cout << "Molecule name = " << Molecule2_AnH.name << endl;
  T = 1000; T1 = 1000;
  cout << "Calculation of the vibrational heat capacity at T1 =" << T1 << " and different values of T [1000, 20000], anharmonic spectrum." << endl;
  cout << "- c_vibr_T(Anharmonic, O2)  " << endl;
  while (T <= 20000) {
   std::cout << -approx.c_vibr_T(T, T1, Molecule2_AnH) << endl;
   T += 500;
  }
  std::cout << endl;
  T = 1000; T1 = 1000;
  cout << "Calculation of the vibrational heat capacity at T1 =" << T1 << " and different values of T [1000, 20000], harmonic spectrum." << endl;
  cout << " c_vibr_T(Harmonic, O2)  " << endl;
  while (T <= 20000) {
    std::cout << approx.c_vibr_T(T, T1, Molecule2_H) << endl;
    T += 500;
  }
  std::cout << endl;
  std::cout << "***************************************************************" << endl;
  T = 1000; T1 = 1000;
  cout << "Calculation of the vibrational heat capacity at T1 =" << T1 << " and different values of T [1000, 20000], anharmonic spectrum." << endl;
  cout << " c_vibr_T1(Anharmonic, O2)  " << endl;
  while (T <= 20000) {
    std::cout << approx.c_vibr_T1(T, T1, Molecule2_AnH) << endl;
    T += 500;
  }
  std::cout << endl;
  T = 1000; T1 = 1000;
  cout << "Calculation of the vibrational heat capacity at T1 =" << T1 << " and different values of T [1000, 20000], harmonic spectrum." << endl;
  cout << " c_vibr_T1(Harmonic, O2)  " << endl;
  while (T <= 20000) {
    std::cout  << approx.c_vibr_T1(T, T1, Molecule2_H) << endl;
    T += 500;
  }
  std::cout << endl;
  std::cout << endl;
  std::cout << "***************************************************************" << endl;
  std::cout << "***************************************************************" << endl;
  T = 5000; T1 = 1000;
  cout << "Calculation of the vibrational heat capacity at T =" << T << " and different values of T1 [1000, 20000], anharmonic spectrum." << endl;
  cout << "- c_vibr_T(Anharmonic, O2)  " << endl;
  while (T1 <= 20000) {
    std::cout  <<- approx.c_vibr_T(T, T1, Molecule2_AnH) << endl;
    T1 += 500;
  }
  std::cout << endl;
  T = 5000; T1 = 1000;
  cout << "Calculation of the vibrational heat capacity at T =" << T << " and different values of T1 [1000, 20000], harmonic spectrum." << endl;
  cout << " c_vibr_T(Harmonic, O2)  " << endl;
  while (T1 <= 20000) {
   std::cout << "T1=" << T1 << " c_vibr_T=" << approx.c_vibr_T(T, T1, Molecule2_H) << endl;
   T1 += 500;
  }
  std::cout << endl;
  std::cout << "***************************************************************" << endl;
  T = 5000; T1 = 1000;
  cout << "Calculation of the vibrational heat capacity at T =" << T << " and different values of T1 [0, 20000], anharmonic spectrum." << endl;
  cout << " c_vibr_T1(Anharmonic, O2)  " << endl;
  while (T1 <= 20000) {
    std::cout << approx.c_vibr_T1(T, T1, Molecule2_AnH) << endl;
    T1 += 500;
  }
  std::cout << endl;
  
  T = 5000; T1 = 1000;
  cout << "Calculation of the vibrational heat capacity at T =" << T << " and different values of T1 [0, 20000], harmonic spectrum." << endl;
  cout << " c_vibr_T1(Harmonic, O2)  " << endl;
  while (T1 <= 20000) {
    std::cout  << approx.c_vibr_T1(T, T1, Molecule2_H) << endl;
    T1 += 500;
  }
}
