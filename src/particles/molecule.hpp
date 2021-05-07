/*! 
    \file molecule.hpp
    \brief The Molecule class loads and stores spectroscopic data for diatomic molecules
 */

#ifndef molecule_hpp
#define molecule_hpp

#include <map>
#include <string>
#include <vector>
#include "particle.hpp"
#include "constants.h"

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>

namespace kappa {

class Molecule : public Particle {

 public:

  /*!
     \brief Creates an object of type Molecule, loading data on a diatomic molecule from the database and calculating rotational and vibrational spectra for each electronic state
     @param const std::string &name molecule name
     @param bool anharmonic_spectrum if the molecule has anharmonic spectrum (the default value is "true")
     @param bool rigid_rotator if the molecule is modelled as a rigid rotator (the default is "true")
     @param const std::string &filename path to particle.yaml database file
   */
  Molecule(const std::string &name, bool anharmonic_spectrum=true, bool rigid_rotator=true, const std::string &filename="particles.yaml");

  // if the molecule has anharmonic spectrum
  bool anharmonic_spectrum;

  // if the molecule is modelled as a rigid rotator
  bool rigid_rotator;

  //! The reduced mass of the oscillator
  //! \f$m_A \times m_B/m_A+m_B\f$ where \f$m_A\f$, \f$m_B\f$ are the mass of atom and molecule, respectively.
  double reduced_osc_mass;

  // The ratio of the mass of the first atom of the molecule to the mass of the molecule
  // \f$m_A/m_A+m_B\f$ where \f$m_A\f$, \f$m_B\f$ are the mass of atom and molecule, respectively.
  double mA_mAB;

  // The ratio of the mass of the first atom of the molecule to the mass of the molecule 
  // \f$m_B/m_A+m_B\f$ where \f$m_A\f$, \f$m_B\f$ are the mass of atom and molecule, respectively.
  double mB_mAB;

  // Rotational inertia moment of the molecule
  double rot_inertia;

  // internuclear distance
  double internuclear_distance;

  // The symmetry factor of molecules (equal to 2 for homonuclear and 1 for heteronuclear molecules)
  int rot_symmetry;

  // The vibrational frequency vector (for each electronic state)
  arma::vec vibr_frequency;

  // The vector of the anharmonicity parameter \f$\omega_ex_e\f$ (for each electronic state)
  arma::vec vibr_we_xe;

  // The vector of the anharmonicity parameter \f$\omega_ey_e\f$ (for each electronic state)
  arma::vec vibr_we_ye;

  // The vector of the anharmonicity parameter \f$\omega_ez_e\f$ (for each electronic state)
  arma::vec vibr_we_ze;

  // Vector of parameter \f$B_e\f$  
  // which describes the dependence of the rotational energy on the rotational level (for each electronic state)
  arma::vec rot_be;

  // Vector of parameter \f$\alpha_e\f$
  // which describes the dependence of the rotational energy on the rotational and vibrational levels (for each electronic state)
  arma::vec rot_ae;
  
  // Vector of vibrational characteristic temperatures (for each electronic state)
  arma::vec characteristic_vibr_temperatures;

  // Vector of dissociation energy (for each electronic state)
  arma::vec diss_energy;

  /*! 
      Array of energy vectors of rotational levels:
      1st index - number of electronic levels 
      2nd index - number of vibrational levels; 
      if the molecule is a rigid rotator, for any fixed 1st index all values are the same 
      (the rotational spectrum does not depend on the vibrational state of the molecule)
   */
  std::vector<std::vector<int> > num_rot_levels; 

  // Defines the vibrational energy levels (moved here from .cpp)
  std::vector<double> level_vibr_energies; 

  // The array of the number of vibrational levels (for each electronic state)
  std::vector<int> num_vibr_levels;

  /*! 
      Array of energy vectors of rotational levels:
      1st index - number of electronic levels 
      2nd index - number of vibrational levels; 
      if the molecule is a rigid rotator, for any fixed 1st index all values are the same 
      (the rotational spectrum does not depend on the vibrational state of the molecule)
   */
  std::vector<std::vector<arma::vec> > rot_energy;

  // The vector of energy vectors of vibrational levels (the index is the number of the electronic level)
  std::vector<arma::vec> vibr_energy;

  // The values of the constant in the Parker formula used to calculate the rotational relaxation times
  double parker_const;
  
 private:

  // overloading of Particle::readData
  void readData(const std::string &name, const std::string &filename);

}; // Molecule
} // namespace kappa
#endif /* molecule_hpp */
