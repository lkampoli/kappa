/*
  \file basicspeedtest.cpp
*/

#include <iostream>
#include <ctime>
#include <chrono>

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#include "kappa.hpp"

using namespace std;
using namespace kappa;

std::string GetCurrentWorkingDir(void) {
  char buff[FILENAME_MAX];
  GetCurrentDir(buff, FILENAME_MAX);
  std::string current_working_dir(buff);
  return current_working_dir;
}

int main(int argc, char** argv) {

  std::cout << "Start Test for timing ..." << endl;
  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY"); 
  std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
  std::string particle_source    = m_source + "particles.yaml";
  std::string interaction_source = m_source + "interaction.yaml";
  std::string output_dir = GetCurrentWorkingDir();
  std::cout << "Current directory is: " << output_dir << std::endl;

  Approximation ApproximationTest{};

  double n = 1.25 / (2 * 23.2651287e-27); //particle number density of N2

  //  Molecule N2("N2"); FIXME
  //  kappa::Molecule N2("N2", "true", "true", particle_source);
  Molecule N2("N2", "true", "true", particle_source);

  //  Interaction int_m1_m2(N2, N2); FIXME
  Interaction int_m1_m2(N2, N2, interaction_source);

  double res1=0., res2=0;
  // std::vector<double> T_vals = {200., 300., 400., 500., 600., 700., 800., 900., 1000.};
  std::string curr_model;
  ofstream out_m1_m2;
    
  clock_t t;
  // t = clock();
  auto begin = std::chrono::high_resolution_clock::now();
  for (int i=0; i<20000; i++) {
    res1 += ApproximationTest.Z_coll(2000. + i, 1e23, int_m1_m2);
  }
  // t = clock()-t;
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  cout << duration << "ns total, average : " << duration / 20000 << "ns." << std::endl;
  // cout << "Res: " << res1 << "clicks: " << t << ", time: " << float(t) / CLOCKS_PER_SEC << " seconds" << std::endl;

  // t = clock();
  begin = std::chrono::high_resolution_clock::now();
  for (int i=0; i<20000; i++) {
    res2 += ApproximationTest.Z_coll(2000. + i, 1e23, int_m1_m2);
  }
  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  cout << duration << "ns total, average : " << duration / 20000 << "ns." << std::endl;
  // t = clock()-t;
    
  // cout << "New optimised Res: " << res2 << "clicks: " << t << ", time: " << float(t) / CLOCKS_PER_SEC << " seconds" << std::endl;
  // cout << "Z=" <<  ApproximationTest.Z_coll(2500, n, int_m1_m2) << ", Z_opt=" <<  ApproximationTest.Z_coll2(2500, n, int_m1_m2) << std::endl;
  cout << "Z=" <<  ApproximationTest.Z_coll(2500, n, int_m1_m2) << ", Z_opt=" <<  ApproximationTest.Z_coll(2500, n, int_m1_m2) << std::endl;

  return 0;
}
