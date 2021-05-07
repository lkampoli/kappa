/*! 
    \file particle.cpp
    \brief Definition of the superclass Particle.
 */

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip> 
#include <cstdlib>
#include "particle.hpp"
#include <numeric>
#include "constants.h"
#include "exceptions.hpp"
#include "yaml-cpp/yaml.h"

using namespace std;
using namespace kappa;

 //! \brief Particle constructor: implements a particle object
 //! @param name name of the particle to instantiate
 //! @param filename name of the database file (particle.yaml)
 Particle::Particle(const std::string &name, const std::string &filename) {
   readData(name, filename);
 }	

 Particle::Particle() {};

 // read particle.yaml database file
 void Particle::readData(const std::string &name, const std::string &filename) {
	
   YAML::Node file;

   try {
     file = YAML::LoadFile(filename);
   } catch (const YAML::BadFile &e) {
     std::string error_string = "Could not load database file " + filename;
     throw UnopenedFileException(error_string);
   }
    
   if (!file[name]) {
     std::string error_string = "No data found for " + name + " in the database";
     throw DataNotFoundException(error_string);
   }

   this->name = name;
   YAML::Node particle = file[name];

  if (particle	["Mass, kg"]) mass                        = particle["Mass, kg"].as<double>();              //!< particle mass
  if (particle	["Diameter, m"]) diameter                 = particle["Diameter, m"].as<double>();           //!< particle diameter
  if (particle	["Formation energy, J"]) formation_energy = particle["Formation energy, J"].as<double>();   //!< particle formation energy
  if (particle	["Charge"]) charge                        = particle["Charge"].as<int>();                   //!< particle charge
  if (particle	["Type"]) particleType                    = particle["Type"].as<string>();                  //!< particle type

  // needed to compute molecular weights & Co.
  if (particle["Stoichiometry"]) {
    YAML::Node stoi = particle["Stoichiometry"];  
    for (YAML::const_iterator it=stoi.begin();it!=stoi.end();++it) {
      stoichiometry.push_back(std::make_pair( it->first.as<std::string>(), it->second.as<int>() ));
    }
  }
	
  // not an electron
  if (name != "e-") {
  if (particle	["Parameter ε (Lennard-Jones), J"]) LennardJones_epsilon = particle["Parameter ε (Lennard-Jones), J"].as<double>();//!< depth of the LJ potential well (if particle = electron, is equal to 0)
  if (particle	["Ionization potential, J"]) ionization_potential = particle["Ionization potential, J"].as<double>(); //!< particle ionization potential (if particle = electron, is equal to 0)
  if (particle	["Electronic energy, J"]) electron_energy = particle["Electronic energy, J"].as<vector<double>>();//!< number of electronic levels of a particle (if particle = electron, is equal to 0)
  if (particle	["Statistical weight"]) statistical_weight = particle["Statistical weight"].as<vector<unsigned long>>(); //!< stat. weights of electron lev. of the particle (if particle = electron, is not initialized)
   num_electron_levels = statistical_weight.size();
  }   
} 
