/*
 * \file test_FHO.cpp
 * \brief Test for FHO model for k_VT and VT probability
 */
 
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>        
#include <stdlib.h>
#include <chrono>
#include <ctime>

#include "kappa.hpp"
 
int main(int argc, char** argv) {
   
  std::cout << "Start test: computation of thermal diffusion coefficients" << std::endl;
  
  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
  std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
  std::string particle_source    = m_source + "particles.yaml";
  std::string interaction_source = m_source + "interaction.yaml";

  std::cout << "Loading particles data" << std::endl;

  kappa::Molecule mol("N2", false, true, particle_source);
  kappa::Atom at("N", particle_source);

  std::cout << "Finished loading particles data" << std::endl;

  std::cout << "Molecule's name " << mol.name << std::endl;
  std::cout << "Atom's name " << at.name << std::endl;
  std::cout << "Molecule vibrational levels " << mol.num_vibr_levels[0] << std::endl;

   kappa::Approximation appr{};
   kappa::Interaction inter_m_a(mol, at, interaction_source);

   // N2(1) + N -> N2(0) + N
//    std::clock_t c_start = std::clock();
//    auto t_start = std::chrono::high_resolution_clock::now();
//  
    std::cout << "probability_VT = " << appr.probability_VT(20000., mol, inter_m_a, 1, -1, 0) << std::endl;
//
//    std::clock_t c_end = std::clock();
//    auto t_end = std::chrono::high_resolution_clock::now();
//
//    std::cout << std::fixed << std::setprecision(2) << "CPU time used: "
//              << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms\n"
//              << "Wall clock time passed: "
//              << std::chrono::duration<double, std::milli>(t_end-t_start).count()
//              << " ms\n";

   double res_kf = 0;
   double res_kb = 0;
   // double T = 10000.0;

   std::cout << "k_VT FHO: " << appr.k_VT(4000., mol, inter_m_a, 1, -1) << std::endl;
   std::cout << "k_VT SSH: " << appr.k_VT(4000., mol, inter_m_a, 1, -1, kappa::models_k_vt::model_k_vt_ssh) << std::endl;
   // std::cout << "k_VT SSH: " << appr.k_VT(3000., mol, inter_m_a, 1, -1, kappa::models_k_vt::model_k_vt_ssh) << std::endl;
   std::cout << "k_VT FHO: " << appr.k_VT(4000., mol, inter_m_a, 2, -1) << std::endl;
   std::cout << "k_VT SSH: " << appr.k_VT(4000., mol, inter_m_a, 2, -1, kappa::models_k_vt::model_k_vt_ssh) << std::endl;
   // std::cout << "k_VT SSH: " << appr.k_VT(3000., mol, inter_m_a, 2, -1, kappa::models_k_vt::model_k_vt_ssh) << std::endl;
   std::cout << "k_VT FHO: " << appr.k_VT(4000., mol, inter_m_a, 10, -1) << std::endl;
   std::cout << "k_VT SSH: " << appr.k_VT(4000., mol, inter_m_a, 10, -1, kappa::models_k_vt::model_k_vt_ssh) << std::endl;
   // std::cout << "k_VT SSH: " << appr.k_VT(3000., mol, inter_m_a, 10, -1, kappa::models_k_vt::model_k_vt_ssh) << std::endl;

   double T = 12000.0;
   double T1 = 10000.0;
   double n = 1e24;
   auto distrib = appr.Boltzmann_distribution(T1, n, mol);

   for (int i=0; i<mol.num_vibr_levels[0]-1; i++) {
     res_kf += i * appr.k_VT(T, mol, inter_m_a, i+1, -1) * distrib[i+1] * n; 					// i+1 -> i
     res_kb += i * appr.k_VT(T, mol, inter_m_a, i+1, -1) * appr.k_bf_VT(T, mol, i+1, -1) * distrib[i] * n; 	// i -> i+1
   }

   for (int i=1; i<mol.num_vibr_levels[0]; i++) {
     res_kf += i * appr.k_VT(T, mol, inter_m_a, i, -1) * appr.k_bf_VT(T, mol, i, -1) * distrib[i-1] * n; 	// i-1 -> i
     res_kb += i * appr.k_VT(T, mol, inter_m_a, i, -1) * distrib[i] * n; 					// i -> i-1
   }
   
   std::cout << "R_VT = " << (res_kf - res_kb) << std::endl;

   return 0;
}
