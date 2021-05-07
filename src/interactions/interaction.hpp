/*!
    \file interaction.hpp
    \brief The Interaction class loads and stores data for particle-particle interactions
 */

#ifndef kappa_interaction_hpp
#define kappa_interaction_hpp

#include <string>
#include <unordered_map>

#include "particles.hpp"
#include "exceptions.hpp"

namespace kappa {

enum class interaction_types {interaction_neutral_neutral, 
			      interaction_neutral_ion, 
			      interaction_neutral_electron, 
			      interaction_charged_charged}; 

class Interaction {

        
 public:

  // Builts an object of type Interaction 
  // loading data on the interaction of two particles from the database
  // @param const Particle &Particle1   - first particle participating in the interaction
  // @param const Particle &Particle2   - second particle participating in the interaction
  // @param const std::string &filename - path to the interaction.yaml database file
  Interaction(const Particle &Particle1, const Particle &Particle2, const std::string &filename = "interaction.yaml");

  // The name of the first particle participating in the interaction
  std::string particle1_name;

  // The name of the second particle participating in the interaction
  std::string particle2_name;

  // Charge of the first particle participating in the interaction (expressed in elementary electric charges)
  int charge1;

  // Charge of the second particle 
  int charge2;

  // Reduced mass of colliding particles
  double collision_mass;

  // Diameter of collision of particles
  double collision_diameter;

  // The depth of the potential well in the Lennard-Jones potential
  double epsilon;

  // Is there data of the VSS potential for this interaction ("true" if present, "false" if there are none)
  bool vss_data;

  // The reference diameter $d_{ref}$ in the potential VSS
  double vss_dref;

  // d = d_{ref} * (g_{ref} / g)^(\omega - 0.5)
  // @param omega the parameter in the potential VSS  
  // @param d the diameter of the collision is expressed in terms of the relative velocity of the colliding particles
  // g, g_{ref} reference speed
  double vss_omega;

  double vss_alpha;

  // Reference temperature T_{ref} in the VSS potential
  double vss_Tref;

  // d = vss_{c,d} g^(0.5 - \omega)
  // Auxiliary value vss_{c,d} 
  // The collision diameter is expressed through it by the formula:
  // d_vss = d_ref * (g_ref / g)^(vss_omega - 0.5) = vss_c_d * g^(0.5 - vss_omega)
  double vss_c_d;

  // \pi * d^2 = vss_{c,cs} g^(1 - 2\omega)
  // Auxiliary value vss_{c,cs}
  // The cross section for the collision is expressed in terms of it by the formula:
  // \pi * d_vss^2 = pi d_ref^2 * (g_ref / g)^(2 * vss_omega - 1) = vss_c_cs * g^(1 - 2 * vss_omega)
  double vss_c_cs; 

   // Type of interaction, which can be:
   // @param kappa::interaction_types::interaction_neutral_neutral  - interaction of uncharged particles
   // @param kappa::interaction_types::interaction_neutral_ion      - interaction of neutral particles and ions
   // @param kappa::interaction_types::interaction_neutral_electron - interaction of neutral particles and electrons
   // @param kappa::interaction_types::interaction_charged_charged  - interaction of two charged particles
   interaction_types interaction_type;

   // arma::mat neutral_electron_elastic_coeffs;
   
   const double& operator[](const std::string &name) const;
   const double& operator[](const char* name) const;

 private:

  std::unordered_map<std::string, double> data;
  void readData(const std::string &name, const std::string &filename);
};
}
#endif /* interaction_hpp */
