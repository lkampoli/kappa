/*
  \file omegaTest.cpp
  \brief Computing omega integrals with different models.
*/

#include <iostream>
#include <iomanip>
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
//  Molecule O2("O2", true, true, particle_source);
//  Molecule NO("NO", true, true, particle_source);

  Atom N("N", particle_source);
//  Atom O("O", particle_source);

  //std::vector<Molecule> molecules = {N2, O2, NO};
  //std::vector<Atom> atoms = {N, O};
  std::vector<Molecule> molecules = {N2};
  std::vector<Atom> atoms = {N};

  std::vector<models_omega> omega_integral_models = {models_omega::model_omega_rs, models_omega::model_omega_vss, 
       						     models_omega::model_omega_bornmayer, models_omega::model_omega_lennardjones, 
                                                     models_omega::model_omega_esa};

  std::map<models_omega, std::string> omega_model_names;
  omega_model_names.insert(std::pair<models_omega, std::string>(models_omega::model_omega_rs, "RS"));
  omega_model_names.insert(std::pair<models_omega, std::string>(models_omega::model_omega_vss, "VSS"));
  omega_model_names.insert(std::pair<models_omega, std::string>(models_omega::model_omega_bornmayer, "Born-Mayer"));
  omega_model_names.insert(std::pair<models_omega, std::string>(models_omega::model_omega_lennardjones, "LJ"));
  omega_model_names.insert(std::pair<models_omega, std::string>(models_omega::model_omega_esa, "ESA"));

  int l = 2;
  int r = 2;
  std::cout << " Omega("<<l<<","<<r<<")" << endl;
  //std::vector<double> T_vals = {500., 1000., 2000., 5000., 10000., 15000., 20000., 25000., 30000., 40000.};
  std::vector<double> T_vals = {500.};

  // for (int i = 0; i < 159; i++) { // assume a max num. of vibr. levels a priori
  //   T_vals.push_back(500 + i * 500);
  // }

  std::vector<double> res_vals;
  std::string curr_model;

  ofstream out_m1_m2;

  // molecule + molecule
  for (int i=0; i<molecules.size(); i++) {
    for (int k=i; k<molecules.size(); k++) {

      std::cout << molecules[i].name + ", " + molecules[k].name << endl;
      Interaction int_m1_m2(molecules[i], molecules[k], interaction_source);
      std::cout << ApproximationTest.omega_integral(500., int_m1_m2, l, r, models_omega::model_omega_rs, true) << std::endl;

      out_m1_m2.open(molecules[i].name + "_" + molecules[k].name + "_" + to_string(l) + "_" + to_string(r) + ".txt");
      out_m1_m2 << std::setw(20) << "Temperature [K]"; 
      out_m1_m2 << std::setw(20) << "Collision diameter"; 
      out_m1_m2 << std::setw(20) << "Collision mass"; 

      for (auto model : omega_integral_models) {
        curr_model = omega_model_names[model];
        out_m1_m2 << std::setw(15) << curr_model << "(" << to_string(l) + "," + to_string(r) + ")" << std::setw(20);
      }
      out_m1_m2 << endl;

      for (auto T : T_vals) {
        out_m1_m2 << std::setw(20) << T << std::setw(20) << int_m1_m2["collision diameter"] << std::setw(20) << int_m1_m2["collision mass"] << std::setw(20);
        for (auto model : omega_integral_models) {

          try {
           out_m1_m2 << std::setw(20) << ApproximationTest.omega_integral(T, int_m1_m2, l, r, model, true) << std::setw(20);
          } catch (const ModelParameterException &e) {
            std::cout << e.what() << endl;
            out_m1_m2 << ";";
          }
        }
        out_m1_m2 << endl;
      }
      out_m1_m2.close();
    }
  }

  // molecule + atom
  for (int i=0; i<molecules.size(); i++) {
    for (int k=0; k<atoms.size(); k++) {

      std::cout << molecules[i].name + ", " + atoms[k].name << endl;
      Interaction int_m1_a1(molecules[i], atoms[k], interaction_source);
      std::cout << ApproximationTest.omega_integral(500., int_m1_a1, l, r, models_omega::model_omega_rs, true) << std::endl;

      out_m1_m2.open(molecules[i].name + "_" + atoms[k].name + "_" + to_string(l) + "_" + to_string(r) + ".txt");
      out_m1_m2 << std::setw(20) << "Temperature [K]";
      out_m1_m2 << std::setw(20) << "Collision diameter";
      out_m1_m2 << std::setw(20) << "Collision mass";

      for (auto model : omega_integral_models) {
        curr_model = omega_model_names[model];
        out_m1_m2 << std::setw(15) << curr_model << "(" << to_string(l) + "," + to_string(r) + ")" << std::setw(20);
      }
      out_m1_m2 << endl;
 
      for (auto T : T_vals) {
        out_m1_m2 << T << std::setw(20) << int_m1_a1["collision diameter"] << std::setw(20) << int_m1_a1["collision mass"] << std::setw(20);
        for (auto model : omega_integral_models) {

          try {
            out_m1_m2 << std::setw(20) << ApproximationTest.omega_integral(T, int_m1_a1, l, r, model, true) << std::setw(20);
          } catch (const ModelParameterException &e) {
            std::cout << e.what() << endl;
            out_m1_m2 << ";";
          }
        }
        out_m1_m2 << endl;
      }
      out_m1_m2.close();
    }
  }

  // atom + atom
  for (int i=0; i<atoms.size(); i++) {
    for (int k=i; k<atoms.size(); k++) {

       std::cout << atoms[i].name + ", " + atoms[k].name << endl;
       Interaction int_a1_a2(atoms[i], atoms[k], interaction_source);
       std::cout << ApproximationTest.omega_integral(500., int_a1_a2, l, r, models_omega::model_omega_rs, true) << std::endl;

       out_m1_m2.open(atoms[i].name + "_" + atoms[k].name + "_" + to_string(l) + "_" + to_string(r) + ".txt");
       out_m1_m2 << std::setw(20) << "Temperature [K]";
       out_m1_m2 << std::setw(20) << "Collision diameter";
       out_m1_m2 << std::setw(20) << "Collision mass";
       for (auto model : omega_integral_models) {
         curr_model = omega_model_names[model];
	 out_m1_m2 << std::setw(15) << curr_model << "(" << to_string(l) + "," + to_string(r) + ")" << std::setw(20);
       }
       out_m1_m2 << endl;
       for (auto T : T_vals) {
         out_m1_m2 << T << std::setw(20) << int_a1_a2["collision diameter"] << std::setw(20) << int_a1_a2["collision mass"] << std::setw(20);
	 for (auto model : omega_integral_models) {
	   try {
	     out_m1_m2 << std::setw(20) << ApproximationTest.omega_integral(T, int_a1_a2, l, r, model, true) << std::setw(20);
	   } catch (const ModelParameterException &e) {
	     std::cout << e.what() << endl;
	     out_m1_m2 << ";";
	   }
	 }
	 out_m1_m2 << endl;
       }
       out_m1_m2.close();
     }
   }

   

   return 0;
}
