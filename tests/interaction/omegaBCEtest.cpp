/*
  \file omegaBCEtest.cpp
*/

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

std::string GetCurrentWorkingDir( void ) {
  char buff[FILENAME_MAX];
  GetCurrentDir( buff, FILENAME_MAX );
  std::string current_working_dir(buff);
  return current_working_dir;
}

using namespace std;
using namespace kappa;

int main(int argc, char** argv) {
 
  std::cout << "Start Test for Omega integrals, loading particle data" << endl;
  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
  std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
  std::string particle_source    = m_source + "particles.yaml";
  std::string interaction_source = m_source + "interaction.yaml";

  Approximation ApproximationTest{};

  Molecule N2("N2", true, true, particle_source);
  // Molecule O2("O2", true, true, particle_source);
  // Molecule NO("NO", true, true, particle_source);

  Atom N("N", particle_source);
  // Atom O("O", particle_source);

  // std::vector<Molecule> molecules = {N2, O2, NO};
  // std::vector<Atom> atoms = {N, O};
  std::vector<Molecule> molecules = {N2};
  std::vector<Atom> atoms = {N};

  std::vector<models_omega> omega_integral_models = {	models_omega::model_omega_vss, 
							models_omega::model_omega_bornmayer, 
							models_omega::model_omega_lennardjones, 
							models_omega::model_omega_esa };
    
  std::map<models_omega, std::string> omega_model_names;

  omega_model_names.insert(std::pair<models_omega, std::string>(models_omega::model_omega_vss, "VSS"));
  omega_model_names.insert(std::pair<models_omega, std::string>(models_omega::model_omega_bornmayer, "Born-Mayer"));
  omega_model_names.insert(std::pair<models_omega, std::string>(models_omega::model_omega_lennardjones, "LJ"));
  omega_model_names.insert(std::pair<models_omega, std::string>(models_omega::model_omega_esa, "ESA"));

  double A, B, C, E, o11, o22;
  // std::vector<double> T_vals = {200., 300., 400., 500., 600., 700., 800., 900., 1000.};
  std::vector<double> T_vals = {500., 1000., 2000., 5000., 10000., 15000., 20000., 25000., 30000., 40000.};
  std::vector<double> res_vals;
  std::string curr_model;
  ofstream out_m1_m2;

  for (int i=0; i<molecules.size(); i++) {
    for (int k=i; k<molecules.size(); k++) {
      std::cout << molecules[i].name + ", " + molecules[k].name << endl;
      Interaction int_m1_m2(molecules[i], molecules[k], interaction_source);
      out_m1_m2.open(molecules[i].name + "_" + molecules[k].name + ".txt");
      out_m1_m2 << "T;";
      for (auto model : omega_integral_models) {
        curr_model = omega_model_names[model];
        out_m1_m2 << curr_model << " A*;" << curr_model << " B*;" << curr_model << " C*;" << curr_model << " E*;";
      }
      out_m1_m2 << "A* (F-K approx); B* (Wright); C* (Wright); E* (F-K approx)" << endl;
      for (auto T : T_vals) {
        out_m1_m2 << T << ";";
        for (auto model : omega_integral_models) {
          try {
            o11 = ApproximationTest.omega_integral(T, int_m1_m2, 1, 1, model, false);
            o22 = ApproximationTest.omega_integral(T, int_m1_m2, 2, 2, model, false);
            A = o22 / o11;
            B = (5 * ApproximationTest.omega_integral(T, int_m1_m2, 1, 2, model, false) - 4 * ApproximationTest.omega_integral(T, int_m1_m2, 1, 3, model, false)) / o11;
            C = ApproximationTest.omega_integral(T, int_m1_m2, 1, 2, model, false) / o11;
            E = ApproximationTest.omega_integral(T, int_m1_m2, 2, 3, model, false) / o22;
            out_m1_m2 << A << ";" << B << ";" << C << ";" << E << ";";
          } catch (const ModelParameterException &e) {
            std::cout << e.what() << endl;
            out_m1_m2 << ";";
          }
        }
	out_m1_m2 << 1.12 << ";" << 1.15 << ";" << 0.92 << ";" << 0.96 << endl;
      }
      out_m1_m2.close();
    }
  }

  for (int i=0; i<molecules.size(); i++) {
    for (int k=0; k<atoms.size(); k++) {
      Interaction int_m1_a1(molecules[i], atoms[k], interaction_source);
      std::cout << molecules[i].name + ", " + atoms[k].name << endl;
      out_m1_m2.open(molecules[i].name + "_" + atoms[k].name + ".txt");
      out_m1_m2 << "T;";
      for (auto model : omega_integral_models) {
        curr_model = omega_model_names[model];
        out_m1_m2 << curr_model << "_A;" << curr_model << "_B;" << curr_model << "_C;" << curr_model << "_E;";
      }
      out_m1_m2 << "A* (F-K approx); B* (Wright); C* (Wright); E* (F-K approx)" << endl;
      for (auto T : T_vals) {
        out_m1_m2 << T << ";";
        for (auto model : omega_integral_models) {
          try {
            o11 = ApproximationTest.omega_integral(T, int_m1_a1, 1, 1, model, false);
            o22 = ApproximationTest.omega_integral(T, int_m1_a1, 2, 2, model, false);
            A = o22 / o11;
            B = (5 * ApproximationTest.omega_integral(T, int_m1_a1, 1, 2, model, false) - 4 * ApproximationTest.omega_integral(T, int_m1_a1, 1, 3, model, false)) / o11;
            C = ApproximationTest.omega_integral(T, int_m1_a1, 1, 2, model, false) / o11;
            E = ApproximationTest.omega_integral(T, int_m1_a1, 2, 3, model, false) / o22;
            out_m1_m2 << A << ";" << B << ";" << C << ";" << E << ";";
          } catch (const ModelParameterException &e) {
            std::cout << e.what() << endl;
            out_m1_m2 << ";";
          }
        }
        out_m1_m2 << 1.12 << ";" << 1.15 << ";" << 0.92 << ";" << 0.96 << endl;
      }
      out_m1_m2.close();
    }
  }

  for (int i = 0; i<atoms.size(); i++) {
    for (int k = i; k<atoms.size(); k++) {
      Interaction int_a1_a2(atoms[i], atoms[k], interaction_source);
      std::cout << atoms[i].name + ", " + atoms[k].name << endl;
      out_m1_m2.open(atoms[i].name + "_" + atoms[k].name + ".txt");
      out_m1_m2 << "T;";
      for (auto model : omega_integral_models) {
        curr_model = omega_model_names[model];
        out_m1_m2 << curr_model << "_A;" << curr_model << "_B;" << curr_model << "_C;" << curr_model << "_E;";
      }
      out_m1_m2 << "A* (F-K approx); B* (Wright); C* (Wright); E* (F-K approx)" << endl;
      for (auto T : T_vals) {
        out_m1_m2 << T << ";";
        for (auto model : omega_integral_models) {
          try {
            o11 = ApproximationTest.omega_integral(T, int_a1_a2, 1, 1, model, false);
            o22 = ApproximationTest.omega_integral(T, int_a1_a2, 2, 2, model, false);
            A = o22 / o11;
            B = (5 * ApproximationTest.omega_integral(T, int_a1_a2, 1, 2, model, false) - 4 * ApproximationTest.omega_integral(T, int_a1_a2, 1, 3, model, false)) / o11;
            C = ApproximationTest.omega_integral(T, int_a1_a2, 1, 2, model, false) / o11;
            E = ApproximationTest.omega_integral(T, int_a1_a2, 2, 3, model, false) / o22;
            out_m1_m2 << A << ";" << B << ";" << C << ";" << E << ";";
          } catch (const ModelParameterException &e) {
            std::cout << e.what() << endl;
            out_m1_m2 << ";";
          }
        }
	out_m1_m2 << 1.12 << ";" << 1.15 << ";" << 0.92 << ";" << 0.96 << endl;
      }
      out_m1_m2.close();
    }
  }

  // just to load only one omega model
  auto model = models_omega::model_omega_vss;
  Interaction int_m1_a1(molecules[0], atoms[0], interaction_source);
  std::cout << "LC-test: " << ApproximationTest.omega_integral(1000, int_m1_a1, 1, 1, model, false) << std::endl;   

  return 0;
}
