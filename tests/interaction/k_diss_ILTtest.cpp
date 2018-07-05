/*
   \file 
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

#include "models.h"
#include "kappa.hpp"

std::string GetCurrentWorkingDir( void ) {
  char buff[FILENAME_MAX];
  GetCurrentDir( buff, FILENAME_MAX );
  std::string current_working_dir(buff);
  return current_working_dir;
}


using namespace std;
using namespace kappa;

void calculate(Molecule m, Atom a, Approximation appr) {

  std::cout << "Loading interaction parameters for " << m.name << "+" << a.name << endl;

  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
  std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
  std::string particle_source    = m_source + "particles.yaml";
  std::string interaction_source = m_source + "interaction.yaml";
  std::string output_dir = GetCurrentWorkingDir();
  std::cout << "Current directory is: " << output_dir << std::endl;

  Interaction inter(m, a, interaction_source);

  ofstream outf;
  
  std::cout << "Calculating k_diss" << endl;

  std::vector<double> T_vals = { 500., 1000., 2000., 5000., 10000., 15000., 20000., 25000., 30000., 40000. };
  std::vector<int> i_vals;

  // std::vector<std::string> non_TM_model_names = {"RS, center-of-mass line energy + vibr energy", "RS, full energy + vibr energy",
  //                                                "VSS, center-of-mass line energy + vibr energy", "VSS, full energy + vibr energy",
  //                                                "Phys4Entry approximation", "ILT"};
  //
  // std::vector<kappa::models_k_diss> non_TM_model_vals = {models_k_diss::model_k_diss_rs_thresh_cmass_vibr, models_k_diss::model_k_diss_rs_thresh_vibr,
  //                                                        models_k_diss::model_k_diss_vss_thresh_cmass_vibr, models_k_diss::model_k_diss_vss_thresh_vibr,
  //                                                        models_k_diss::model_k_diss_phys4entry, models_k_diss::model_k_diss_ilt};

  std::vector<std::string> non_TM_model_names = {"ILT"};
  std::vector<kappa::models_k_diss> non_TM_model_vals = {models_k_diss::model_k_diss_ilt};

  // int cc = 0;
  // for (auto diss_model: non_TM_model_vals) {
  //   cout << non_TM_model_names[cc] << endl;
  //   cc++;
  // }

  if (m.name == "O2") {
    i_vals = { 0, 5, 10, 20 };
  } else {
    i_vals = { 0, 5, 10, 20, 30 };
  }

  int counter;
  double res, tmp;
  outf.open(m.name + "_" + a.name + "_k_ILT_test.txt");
  outf << std::setw(20) << "T [K]";
  outf << std::setw(20) << "vibr. level";
  outf << std::setw(20) << "vibr. energy (eV)";
  outf << std::setw(20) << "ILT";
  outf << std::setw(20) << "STS";
  outf << std::endl;

  double c1, c2;
    
  for (auto T : T_vals) {
    for (auto i : i_vals) {

      outf << std::setw(20) << T;
      outf << std::setw(20) << i;
      outf << std::setw(20) << m.vibr_energy[0][i] / K_CONST_EV;
      outf << std::setw(20) << appr.k_diss(T, m, inter, i, models_k_diss::model_k_diss_ilt);

      if (m.name == "O2") {
        c1 = 0.3867 * i * i * i - 2.7425 * i * i - 1901.9 * i + 61696;
        if (i <= 31) {
          c2 = 1.63e-9 * i * i * i - 1.25e-7 * i * i + 3.24e-6 * i + 7.09e-5;
        } else if (i<=37) {
          c2 = - 6.67e-6 * i * i + 4.65e-4 * i - 7.91e-3;
        } else {
          c2 = 7.83e-7 * i * i * i * i - 1.31e-4 * i * i * i + 8.24e-3 * i * i - 0.23 * i + 2.4049;
        }
 
        outf << std::setw(20) << pow(T, -0.1) * 1.53e-10 * c2 * exp(-c1 / T);
        outf << endl;
 
      } else {

        if (i <= 8) {
          c1 = 1.786e-18;
        } else if (i<=34) {
          c1 = 1.71e-18;
        } else if (i<=52) {
          c1 = 1.68e-18;
        } else {
          c1 = 1.66e-18;
        }
        c2 = 4e-19 * i * i * i * i + 5.24e-19 * i * i * i - 7.41e-17 * i * i + 6.42e-15 * i + 7.3e-14;
        outf << std::setw(20) << 7.16e-2 * pow(T, -0.25) * c2 * exp((m.vibr_energy[0][i] - c1) / (K_CONST_K * T));
        outf << endl;
      }
    }
  }
}

int main(int argc, char** argv) {

  std::cout << "Start Test for k_diss coefficients, loading particle data" << endl;

  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
  std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
  std::string particle_source    = m_source + "particles.yaml";
  std::string interaction_source = m_source + "interaction.yaml";
  std::string output_dir = GetCurrentWorkingDir();
  std::cout << "Current directory is: " << output_dir << std::endl;

  Molecule N2("N2", true, "harmonic", particle_source);
  // Molecule N2("N2", false, true, particle_source);
  Molecule O2("O2", true, "harmonic", particle_source);
  Atom N("N", particle_source);
  Atom O("O", particle_source);

  Approximation ApproximationTest{};

  calculate(N2, N, ApproximationTest);
  calculate(O2, O, ApproximationTest);

  return 0;
}
