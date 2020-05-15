/*!
    \file molecule.cpp
    \brief Definition of Molecule sublclass
 */

#include <fstream>
#include "molecule.hpp"
#include "exceptions.hpp"
#include "yaml.h"
#include "yaml-cpp/yaml.h"

using namespace std;
using namespace kappa;

Molecule::Molecule(const std::string &name, bool anharmonic_spectrum, bool rigid_rotator, const std::string &filename) : Particle(name, filename) {

  readData(name, filename);

    int e, i, j;
    double tmp, level_we_xe, level_we_ye, level_we_ze;

    this->rigid_rotator = rigid_rotator;
     
    // calculation of energy levels goes here
    for (e = 0; e < num_electron_levels; e++) {
      if (anharmonic_spectrum == false) {
        level_we_xe = 0.0;
        level_we_ye = 0.0;
        level_we_ze = 0.0;
        this->anharmonic_spectrum = false;
      } else {
        level_we_xe = vibr_we_xe[e];
        level_we_ye = vibr_we_ye[e];
        level_we_ze = vibr_we_ze[e];
        this->anharmonic_spectrum = true;
      }

      tmp=0.0; // vibrational energy handler
      i=0;     // vibrational levels handler

      while (tmp < diss_energy[e] - electron_energy[e]) {

        // Morse potential (the anharmonic oscillator model), eq. 1.3.
        tmp = K_CONST_H * K_CONST_C * (vibr_frequency[e] * (i + 0.5) - level_we_xe * (i + 0.5) * (i + 0.5) 
                                                                     + level_we_ye * (i + 0.5) * (i + 0.5) * (i + 0.5) 
                                                                     + level_we_ze * (i + 0.5) * (i + 0.5) * (i + 0.5) * (i + 0.5));
        if (tmp < diss_energy[e] - electron_energy[e]) {
          if (i==0) {
            // Adds the tmp value at the end of the vector level_vibr_energies
            level_vibr_energies.push_back(tmp);
            i++;	
          } else {
            if (tmp > level_vibr_energies[i - 1]) {
              level_vibr_energies.push_back(tmp);
              i++;	
            } else {
       	      break;
            }
          }
        }
      }

      // Adds the i-th vibrational level at the end of the vector num_vibr_levels
      num_vibr_levels.push_back(i);
      vibr_energy.push_back(arma::vec(level_vibr_energies));

      if (rigid_rotator==true) {

        std::vector<double> ev_level_rot_energies;
        std::vector<arma::vec> e_level_rot_energies;
        std::vector<int> e_level_rot_amt;
     
        tmp=0.0;
        j=0;

        while (tmp < diss_energy[e] - electron_energy[e]) {
          tmp = K_CONST_H * K_CONST_C * (rot_be[e] * j * (j+1));
          if (tmp < diss_energy[e] - electron_energy[e]) {
            if (j==0) {
              ev_level_rot_energies.push_back(tmp);
              j++;	
            } else {
              if (tmp > ev_level_rot_energies[j-1]) {
                ev_level_rot_energies.push_back(tmp);
                j++;	
              } else {
                break;
              }
            }
          }
       }

       for (i=0; i<num_vibr_levels[e]; i++) {
         e_level_rot_amt.push_back(j);
         e_level_rot_energies.push_back(arma::vec(ev_level_rot_energies));
       }
       rot_energy.push_back(e_level_rot_energies);
       num_rot_levels.push_back(e_level_rot_amt);

      } else { // non-rigid rotator
        std::vector<arma::vec> e_level_rot_energies;
        std::vector<int> ev_level_rot_amt;

        for (i=0; i<num_vibr_levels[e]; i++) {
          std::vector<double> ev_level_rot_energies;
          tmp=0.;
          j=0;

          while (tmp + vibr_energy[e][i] < diss_energy[e] - electron_energy[e]) {
            tmp = K_CONST_H * K_CONST_C * ((rot_be[e] - rot_ae[e] * (i + 0.5)) * j * (j + 1));
            if (tmp + vibr_energy[e][i] < diss_energy[e] - electron_energy[e]) {
              if (j==0) {
                ev_level_rot_energies.push_back(tmp);
                j++;	
              } else {
                if (tmp > ev_level_rot_energies[j - 1]) {
                  ev_level_rot_energies.push_back(tmp);
                  j++;	
                } else {
                  break;
                }
              }
            }
          }
         
          ev_level_rot_amt.push_back(j);
          e_level_rot_energies.push_back(arma::vec(ev_level_rot_energies));
        } 

        rot_energy.push_back(e_level_rot_energies);
        num_rot_levels.push_back(ev_level_rot_amt);
      }
    }
    // formula p. 10
    characteristic_vibr_temperatures = K_CONST_H * K_CONST_C * vibr_frequency / K_CONST_K;

} // Molecule

  // reads molecule particles
  void Molecule::readData(const std::string &name, const std::string &filename) {

    YAML::Node file;
    
    try {
      file = YAML::LoadFile(filename);
    } catch (const YAML::BadFile &e) {
      std::string error_string = "Could not open database file " + filename;
      throw UnopenedFileException(error_string);
      return;
    }

    if (!file[name]) {
      std::string error_string = "No data found for " + name + " in the database";
      throw DataNotFoundException(error_string);
    }

    this->name = name;
    YAML::Node particle = file[name];
    
    // TODO something if does not exist in file
    if (!particle["Dissociation energy, J"]) {
      std::string error_string = "No data found for " + name + " in the database";
      throw DataNotFoundException(error_string);
    }
    parker_const = 0;

    if (particle["Dissociation energy, J"]) diss_energy = particle["Dissociation energy, J"].as<vector<double>>();
    if (particle["Factor of symmetry"]) rot_symmetry = particle["Factor of symmetry"].as<int>();
    if (particle["Parker (zeta^infty)"]) parker_const = particle["Parker (zeta^infty)"].as<double>();
    if (particle["Reduced oscillator mass"]) reduced_osc_mass = particle["Reduced oscillator mass"].as<double>();
    if (particle["Moment of Inertia, J*s^2"]) rot_inertia = particle["Moment of Inertia, J*s^2"].as<double>();
    if (particle["Frequency of vibrations (we), m^-1"]) vibr_frequency = particle["Frequency of vibrations (we), m^-1"].as<vector<double>>();
    if (particle["wexe, m^-1"]) vibr_we_xe = particle["wexe, m^-1"].as<vector<double>>();
    if (particle["weye, m^-1"]) vibr_we_ye = particle["weye, m^-1"].as<vector<double>>();
    if (particle["weze, m^-1"]) vibr_we_ze = particle["weze, m^-1"].as<vector<double>>();
    if (particle["Be, m^-1"]) rot_be = particle["Be, m^-1"].as<vector<double>>();
    if (particle["ae, m^-1"]) rot_ae = particle["ae, m^-1"].as<vector<double>>();
    if (particle["internuclear distance, r_e , m"]) internuclear_distance = particle["internuclear distance, r_e , m"].as<double>();
    if (particle["mA/mAB"]) mA_mAB = particle["mA/mAB"].as<double>();
    if (particle["mB/mAB"]) mB_mAB = particle["mB/mAB"].as<double>();
    
  } // Molecule::readData 
