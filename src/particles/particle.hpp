/*! 
    \file particle.hpp
    \brief Declaration of the superclass Particle.
 
    The Particle class loads and stores basic data for particles; 
    the only particle that can be loaded only as a Particle object is the electron
 */

#ifndef particle_hpp
#define particle_hpp

#include <string>
#include <vector>
#include <armadillo>

namespace kappa {

class Particle {
    
 public:
 
  //! \brief Creates an Particle-type object by loading the particle data from the database
  //! @param name particle name 
  //! @filename path to the database file
  Particle(const std::string &name, const std::string &filename="particles.yaml");

  //! Creates an empty object of type Particle with all parameters initialized to zero
  Particle();

  //! Name of the particle
  std::string name;

  //! Mass of the particle
  double mass = 0.0;

  //! Diameter of the particle
  double diameter = 0.0; 

  //! The charge of a particle (expressed in elementary electric charges)
  int charge = 0;

  //! Formation energy of the particle
  double formation_energy = 0.0;
 
  //! The depth of the potential well in the Lennard-Jones potential (if the particle is an electron, is assumed to be 0)
  double LennardJones_epsilon = 0.0;

  //! The particle ionization potential (if the particle is an electron, is considered equal to 0)
  double ionization_potential = 0.0;

  //! The number of electronic levels of a particle (if the particle is an electron, is considered equal to 0)
  int num_electron_levels = 0;

  //! The energy vector of the electron levels of the particle (if the particle is an electron, is not initialized)
  arma::vec electron_energy;

  //! Vector of statistical weights of the electron levels of the particle (if the particle is an electron, is not initialized)
  arma::Col<unsigned long> statistical_weight;
 
  //! Type of particle
  std::string particleType;

  //! Stoichiometry of the species
  std::vector<std::pair<std::string,int>> stoichiometry;

  int atomic_number = 0;
  double atomic_weight = 0.0;

  //! Pair vector to store all atomic element in the database file (particles.yaml)
  std::vector<std::pair<std::string,double>> element_list;

 private:

  void readData(const std::string &name, const std::string &filename);

};
} 
# endif 
