
#include <iostream>
#include <fstream>

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#include "kappa.hpp"

  std::string GetCurrentWorkingDir( void ) 
  {
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
   
  Interaction inter(m, a, interaction_source);
   
  double rel_vel = 0.0;
  ofstream outf;
    
  std::cout << "Calculating sigma_diss" << endl;

  std::vector<double> t_vals = { 1000., 2000., 4000., 5000., 10000., 15000., 20000., 25000., 30000., 45000.,
                                 50000., 55000., 60000., 65000., 70000., 75000., 80000., 85000., 90000., 95000., 100000.,
                                 105000., 110000., 115000., 120000., 125000., 130000., 135000., 140000., 145000., 150000., 155000.,
                                 160000., 165000., 170000., 175000., 180000., 185000., 190000., 195000., 200000.};
  std::vector<int> i_vals;

  if (m.name == "O2") {
    i_vals = { 0, 5, 10, 20 };
  } else {
    i_vals = { 0, 5, 10, 20, 30 };
  }

  std::vector<std::string> cs_model_names = {	"RS",
						"RS, center-of-mass line energy", 
						"RS, center-of-mass line energy + vibr energy",
                                               	"VSS", 
						"VSS, center-of-mass line energy + vibr energy",
                                               	"ILT"};

  std::vector<kappa::models_cs_diss> cs_model_vals = {models_cs_diss::model_cs_diss_rs_thresh, 
                                                      models_cs_diss::model_cs_diss_rs_thresh_cmass, 
						      models_cs_diss::model_cs_diss_rs_thresh_cmass_vibr, 
						      models_cs_diss::model_cs_diss_vss_thresh, 
						      models_cs_diss::model_cs_diss_vss_thresh_cmass_vibr, 
						      models_cs_diss::model_cs_diss_ilt};
  double res, tmp;
  outf.open(m.name + "_" + a.name + ".txt");

  outf << "t;";

  for (auto i : i_vals) {
    int counter = 0;
    for (auto diss_model: cs_model_vals) {
      outf << cs_model_names[counter] << ",i=" << i << "(" << m.vibr_energy[0][i] / K_CONST_EV << " eV);"; 
      counter++;
    }
    outf << ";";
  }

  outf << endl;

  for (auto t : t_vals) {
    outf << t << ";";

   for (auto i : i_vals) {
     for (auto diss_model: cs_model_vals) {
       rel_vel = sqrt(2 * K_CONST_K * t / inter.collision_mass);
       outf << appr.crosssection_diss(rel_vel, m, inter, i, diss_model) << ";"; 
     }
     outf << ";";
   }
   outf << endl;
  }
}

int main(int argc, char** argv) {

  std::cout << "Start Test for dissociation cross sections, loading particle data" << endl;

  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
  std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
  std::string particle_source    = m_source + "particles.yaml";
  std::string interaction_source = m_source + "interaction.yaml";

  kappa::Molecule N2("N2", true,	"harmonic",	particle_source);
  kappa::Molecule O2("O2", true,	"harmonic",	particle_source);
  kappa::Atom N("N",					particle_source);
  kappa::Atom O("O",					particle_source);

  Approximation ApproximationTest{};

  calculate(N2, N, ApproximationTest);
  calculate(O2, O, ApproximationTest);

  return 0;
}
