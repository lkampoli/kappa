/*! 
 *  \file atom.cpp
 *  \brief The Atom class loads and stores spectroscopic data for atoms
 */

#include <fstream>
#include "atom.hpp"
#include "exceptions.hpp"
#include "yaml-cpp/yaml.h"

using namespace std;
using namespace kappa;

 Atom::Atom(const std::string &name, const std::string &filename) : Particle(name, filename) {
   readData(name, filename);
 }

 void Atom::readData(const std::string &name, const std::string &filename) {
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

   if(particle["Atomic_Number"]) atomic_number = particle["Atomic_Number"].as<int>();	//!< atomic number
   if(particle["Atomic_Weight"]) atomic_weight = particle["Atomic_Weight"].as<double>();//!< atomic weight

 } 
