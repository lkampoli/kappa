/*!
    \file interaction.cpp
 */

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <iomanip>

#include "exceptions.hpp"
#include "interaction.hpp"
#include "yaml.h"

using namespace std;
using namespace kappa;

Interaction::Interaction(const Particle &Particle1, const Particle &Particle2, const string &filename)
  : particle1_name(Particle1.name), particle2_name(Particle2.name) {

    vss_Tref  = 273;
    vss_dref  = 0;
    vss_omega = 0;
    vss_alpha = 0;
    vss_data  = false;
    vss_c_d   = 0;
    vss_c_cs  = 0;
    
    charge1 = Particle1.charge;
    charge2 = Particle2.charge;

    try {
      readData(particle1_name + " + " + particle2_name, filename);
    } catch (DataNotFoundException){
      try {
        readData(particle2_name + " + " + particle1_name, filename);
      } catch (DataNotFoundException) {
        vss_data = false;
      }
    }
    if ((Particle1.charge == 0) && (Particle2.charge == 0)) {
      interaction_type = interaction_types::interaction_neutral_neutral;
    } else if (((Particle1.charge != 0) && (Particle2.charge == 0)) || ((Particle1.charge == 0) && (Particle2.charge != 0))) {
      if ((Particle1.name == "e-") || (Particle2.name == "e-")) {
        interaction_type = interaction_types::interaction_neutral_electron;
      } else {
        interaction_type = interaction_types::interaction_neutral_ion;
      }
    } else {
      interaction_type = interaction_types::interaction_charged_charged;
    }

    double cd = 0.5 * (Particle1.diameter + Particle2.diameter);
    double LJe = sqrt(Particle1.LennardJones_epsilon * Particle2.LennardJones_epsilon) * pow(Particle1.diameter * Particle2.diameter, 3) / pow(cd, 6);

    data.insert(std::pair<string, double>("collision mass", Particle1.mass * Particle2.mass / (Particle1.mass + Particle2.mass)));
    collision_mass = Particle1.mass * Particle2.mass / (Particle1.mass + Particle2.mass);

    data.insert(std::pair<string, double>("collision diameter", cd));
    collision_diameter = cd;

    data.insert(std::pair<string, double>("Lennard-Jones epsilon", LJe));  // TODO check if exists
    epsilon = LJe;

    if (vss_data) {
      double gref = sqrt(2 * K_CONST_K * vss_Tref / collision_mass);
      vss_c_d = vss_dref * pow(gref, vss_omega - 0.5);
      vss_c_cs = K_CONST_PI * vss_dref * vss_dref * pow(gref, 2 * vss_omega - 1);
    }
} // Interaction

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Interaction::readData(const string &name, const string &filename) {

  YAML::Node file;
  std::vector<double> tmp; 

  try {file = YAML::LoadFile(filename);} 
  catch (const YAML::BadFile &e) {
    std::string error_string = "Could not load database file " + filename;
    throw UnopenedFileException(error_string);	
  }

  if (!file[name]) {
    try {
      std::string error_string = "No data found for " + name + " interaction in the database";}
      // throw DataNotFoundException(error_string);	}
      catch(int num) {
      // int exceptions thrown in the above try block will be handled here
      cerr << "No data found for " + name + " interaction in the database " << num << endl;}
  }

  YAML::Node interaction = file[name];
  int vss_counter = 0;

  for (YAML::const_iterator it = interaction.begin(); it != interaction.end(); ++it) {
    if (!it->second.IsNull()) {
      if (it->second.IsSequence()) {
        // vector_data.insert(std::pair<string, arma::vec>(it->first.as<std::string>(),));
        tmp = it->second.as<std::vector<double> >();
        for (int i=0; i<tmp.size(); i++) {
          data.insert(std::pair<string, double>("_" + it->first.as<std::string>() + "_" + std::to_string(i), tmp[i]));
        }
      } else {
        data.insert(std::pair<string, double>(it->first.as<std::string>(), it->second.as<double>()));
        if (it->first.as<std::string>() == "VSS, Tref") {
          vss_Tref = it->second.as<double>();
          vss_counter += 1;
        } else if (it->first.as<std::string>() == "VSS, dref") {
          vss_dref = it->second.as<double>();
          vss_counter += 1;
        } else if (it->first.as<std::string>() == "VSS, omega") {
          vss_omega = it->second.as<double>();
          vss_counter += 1;
        } else if (it->first.as<std::string>() == "VSS, alpha") {
          vss_alpha = it->second.as<double>();
          vss_counter += 1;
        }
      }
    }
  }
  if (vss_counter==4) {
    vss_data = true;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const double& kappa::Interaction::operator[](const std::string &name) const {

 try {
   return data.at(name);
 } 
 catch (const std::out_of_range &oor) {
   std::string error_string = "No " + name + " interaction parameter found for " + particle1_name + "+" + particle2_name + " interaction";
   throw kappa::DataNotFoundException(error_string.c_str());	
 }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const double& kappa::Interaction::operator[](const char* name) const {

 std::string conv_str(name);
 try {
   return data.at(conv_str);
 }
 catch (const std::out_of_range &oor) {
   std::string error_string = "No " + conv_str + " interaction parameter found for " + particle1_name + "+" + particle2_name + " interaction";
   throw kappa::DataNotFoundException(error_string.c_str());	
 }
}
  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
