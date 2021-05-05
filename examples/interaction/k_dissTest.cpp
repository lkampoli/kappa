/*
   \file k_dissTest.cpp
   \brief Computing dissociation rate coefficients k_diss according to several models.
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

  if (m.name == "O2") {
    i_vals = { 0, 5, 10, 20 };
  } else {
    i_vals = { 0, 5, 10, 20, 30 };
  }

  std::vector<std::string> ALL_model_names = {"rs_thresh_cmass_vibr", 
                                              "rs_thresh_vibr",
       					      "rs_thresh_cmass_vibr",
                                              "rs_thresh_cmass",
      					      "rs_thresh",
      					      "vss_thresh_cmass_vibr",
      					      "vss_thresh_vibr",
      					      "vss_thresh_cmass",
      					      "vss_thresh",
      					      "arrh_scanlon",
      					      "arrh_park",
      					      "tm_D6k_arrh_scanlon",
      					      "tm_3T_arrh_scanlon",
      					      "tm_infty_arrh_scanlon",	
      					      "tm_D6k_arrh_park",
      					      "tm_3T_arrh_park",
      					      "tm_infty_arrh_park",
      					      "phys4entry",
      					      "ilt"};

  std::vector<kappa::models_k_diss> ALL_model_vals = {models_k_diss::model_k_diss_rs_thresh_cmass_vibr,
      						      models_k_diss::model_k_diss_rs_thresh_vibr,
      						      models_k_diss::model_k_diss_rs_thresh_cmass_vibr,
      						      models_k_diss::model_k_diss_rs_thresh_cmass,
      						      models_k_diss::model_k_diss_rs_thresh,
                                                      models_k_diss::model_k_diss_vss_thresh_cmass_vibr, 
      						      models_k_diss::model_k_diss_vss_thresh_vibr,
      						      models_k_diss::model_k_diss_vss_thresh_cmass,
      						      models_k_diss::model_k_diss_vss_thresh,
      						      models_k_diss::model_k_diss_arrh_scanlon,
      						      models_k_diss::model_k_diss_arrh_park,
      						      models_k_diss::model_k_diss_tm_D6k_arrh_scanlon,
      						      models_k_diss::model_k_diss_tm_3T_arrh_scanlon,
      						      models_k_diss::model_k_diss_tm_infty_arrh_scanlon,
      						      models_k_diss::model_k_diss_tm_D6k_arrh_park,
      						      models_k_diss::model_k_diss_tm_3T_arrh_park,
      						      models_k_diss::model_k_diss_tm_infty_arrh_park,
      						      models_k_diss::model_k_diss_phys4entry, 
      						      models_k_diss::model_k_diss_ilt};

  std::vector<std::string> non_TM_model_names = {"RS, center-of-mass line energy + vibr energy", "RS, full energy + vibr energy",
                                                 "VSS, center-of-mass line energy + vibr energy", "VSS, full energy + vibr energy",
                                                 "Phys4Entry approximation", "ILT"};

  std::vector<kappa::models_k_diss> non_TM_model_vals = {models_k_diss::model_k_diss_rs_thresh_cmass_vibr, models_k_diss::model_k_diss_rs_thresh_vibr,
                                                         models_k_diss::model_k_diss_vss_thresh_cmass_vibr, models_k_diss::model_k_diss_vss_thresh_vibr,
                                                         models_k_diss::model_k_diss_phys4entry, models_k_diss::model_k_diss_ilt};

  // std::vector<std::string> TM_model_names = {"Treanor-Marrone, U=inf, Scanlon", "Treanor-Marrone, U=D/6K, Scanlon", "Treanor-Marrone, U=3T, Scanlon"};
  // std::vector<kappa::models_k_diss> TM_model_vals = {model_k_diss_tm_infty_arrh_scanlon, model_k_diss_tm_D6k_arrh_scanlon,  model_k_diss_tm_3T_arrh_scanlon};

  // std::vector<std::string> TM_model_names = {"Treanor-Marrone, U=inf, Park", "Treanor-Marrone, U=D/6K, Park", "Treanor-Marrone, U=3T, Park"};
  // std::vector<kappa::models_k_diss> TM_model_vals = {models_k_diss::model_k_diss_tm_infty_arrh_park, models_k_diss::model_k_diss_tm_D6k_arrh_park, models_k_diss::model_k_diss_tm_3T_arrh_park};

  double res, tmp;
  outf.open( output_dir + "/DISS-REC/" + m.name + "_" + a.name + ".txt");

  outf << std::setw(20) << "T [K]";
  outf << std::setw(20) << "Vibr. l.";
  outf << std::setw(20) << "k_diss";
  outf << std::endl;

  int counter;

  for (auto T : T_vals) {
    for (auto i : i_vals) {
      for (auto diss_model: ALL_model_vals) {
        outf << std::setw(15) << T << std::setw(20);
        outf << std::setw(15) << i << std::setw(20);
        outf << std::setw(15) << appr.k_diss(T, m, inter, i, models_k_diss::model_k_diss_rs_thresh_cmass_vibr) << std::setw(20);
        outf << std::setw(15) << appr.k_diss(T, m, inter, i, models_k_diss::model_k_diss_rs_thresh_vibr) << std::setw(20);
        outf << std::setw(15) << appr.k_diss(T, m, inter, i, models_k_diss::model_k_diss_rs_thresh_cmass_vibr) << std::setw(20);
        outf << std::setw(15) << appr.k_diss(T, m, inter, i, models_k_diss::model_k_diss_rs_thresh_cmass) << std::setw(20);
        outf << std::setw(15) << appr.k_diss(T, m, inter, i, models_k_diss::model_k_diss_rs_thresh) << std::setw(20);
        outf << std::setw(15) << appr.k_diss(T, m, inter, i, models_k_diss::model_k_diss_vss_thresh_cmass_vibr) << std::setw(20);
        outf << std::setw(15) << appr.k_diss(T, m, inter, i, models_k_diss::model_k_diss_vss_thresh_vibr) << std::setw(20);
        outf << std::setw(15) << appr.k_diss(T, m, inter, i, models_k_diss::model_k_diss_vss_thresh_cmass) << std::setw(20);
        outf << std::setw(15) << appr.k_diss(T, m, inter, i, models_k_diss::model_k_diss_vss_thresh) << std::setw(20);
        outf << std::setw(15) << appr.k_diss(T, m, inter, i, models_k_diss::model_k_diss_arrh_scanlon) << std::setw(20);
        outf << std::setw(15) << appr.k_diss(T, m, inter, i, models_k_diss::model_k_diss_arrh_park) << std::setw(20);
        outf << std::setw(15) << appr.k_diss(T, m, inter, i, models_k_diss::model_k_diss_tm_D6k_arrh_scanlon) << std::setw(20);
        outf << std::setw(15) << appr.k_diss(T, m, inter, i, models_k_diss::model_k_diss_tm_3T_arrh_scanlon) << std::setw(20);
        outf << std::setw(15) << appr.k_diss(T, m, inter, i, models_k_diss::model_k_diss_tm_infty_arrh_scanlon) << std::setw(20);
        outf << std::setw(15) << appr.k_diss(T, m, inter, i, models_k_diss::model_k_diss_tm_D6k_arrh_park) << std::setw(20);
        outf << std::setw(15) << appr.k_diss(T, m, inter, i, models_k_diss::model_k_diss_tm_3T_arrh_park) << std::setw(20);
        outf << std::setw(15) << appr.k_diss(T, m, inter, i, models_k_diss::model_k_diss_tm_infty_arrh_park) << std::setw(20);
        outf << std::setw(15) << appr.k_diss(T, m, inter, i, models_k_diss::model_k_diss_phys4entry) << std::setw(20);
        outf << std::setw(15) << appr.k_diss(T, m, inter, i, models_k_diss::model_k_diss_ilt) << std::setw(20);
        outf << std::endl;
        counter++;
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

  std::cout << "Loading particles data" << std::endl;

  Molecule N2("N2", false, true, particle_source);

  Atom N("N", particle_source);

  Approximation ApproximationTest{};

  calculate(N2, N, ApproximationTest);

   return 0;
}
