/*!
    \file mixture.hpp
 */

#ifndef kappa_mixture_hpp
#define kappa_mixture_hpp

#include <vector>    
#include <map>
#include <numeric>
#include <string>

#define ARMA_DONT_PRINT_ERRORS // to avoid warning message during solve() linear systems

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo> 

#include "approximation.hpp"
#include "numeric.hpp"
#include "models.h"

namespace kappa {

class Mixture : public Approximation { 

  using Approximation::c_tr;
  using Approximation::c_rot;
  using Approximation::debye_length;

 public:
    
  Mixture(const std::vector<kappa::Molecule> &i_molecules, 
          const std::vector<kappa::Atom> &i_atoms, 
          const std::string &interactions_filename="interaction.yaml", 
          const std::string &particles_filename="particles.yaml");

  Mixture(const std::vector<kappa::Molecule> &i_molecules, 
          const std::string &interactions_filename="interaction.yaml", 
          const std::string &particles_filename="particles.yaml");

  Mixture(const std::vector<kappa::Atom> &i_atoms, 
          const std::string &interactions_filename="interaction.yaml", 
          const std::string &particles_filename="particles.yaml");

  Mixture(const kappa::Molecule &molecule, 
          const kappa::Atom &atom, 
          const std::string &interactions_filename="interaction.yaml", 
          const std::string &particles_filename="particles.yaml");

  Mixture(const kappa::Atom &atom, 
          const std::string &interactions_filename="interaction.yaml", 
          const std::string &particles_filename="particles.yaml");

  Mixture(const kappa::Molecule &molecule, 
          const std::string &interactions_filename="interaction.yaml", 
          const std::string &particles_filename="particles.yaml");

  Mixture(const std::string particle_names, 
          const std::string &interactions_filename="interaction.yaml", 
          const std::string &particles_filename="particles.yaml", 
          bool anharmonic=true, 
          bool rigid_rotators=true);

  // Returns a string with the names of particles in the mixture
  std::string get_names();

  // Returns the number of particles in the mixture
  int get_n_particles();

  // Returns the sum of the number of vibrational levels in the mixtures
  int get_n_vibr_levels();

  // Returns an array of the number of vibrational levels of molecules in the mixture
  std::vector<int> get_n_vibr_levels_array();

  // Translates an array of molar fractions into an array of mass fractions
  // @param const arma::vec &x - array of molar fractions of the mixture
  arma::vec convert_molar_frac_to_mass(const arma::vec &x);

  // Converts an array of mass fractions into an array of molar fractions
  // @param const arma::vec &y - array of mass fractions of the mixture
  arma::vec convert_mass_frac_to_molar(const arma::vec &y);

  // Returns an object corresponding to the molecule present in the mixture
  // @param const std::string &name - name of the molecule
  kappa::Molecule molecule(const std::string &name);

  // Returns an object corresponding to the atom present in the mixture
  // @param const std::string &name - name of the atom
  kappa::Atom atom(const std::string &name);
  
  kappa::Interaction interaction(const kappa::Molecule &molecule1, const kappa::Molecule &molecule2);
  kappa::Interaction interaction(const kappa::Molecule &molecule, const kappa::Atom &atom);
  kappa::Interaction interaction(const kappa::Atom &atom, const kappa::Molecule &molecule);
  kappa::Interaction interaction(const kappa::Atom &atom1, const kappa::Atom &atom2);

  double debye_length(double T, const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons);
  double debye_length(double T, const std::vector<arma::vec> &n_vl_molecule, double n_electrons);
  double debye_length(double T, const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons);
  double debye_length(double T, const arma::vec &n);

  // Calculation of the number density of the mixture
  // @param const std::vector<arma::vec> &n_vl_molecule - the population of the vibrational levels of molecules
  // @param const arma::vec &n_molecule                 - Vector of number densities of molecules
  // @param const arma::vec &n_atom                     - Vector of number densities of atoms
  // @param double n_electrons                          - the number density of electrons (the default value is 0)
  // @param const arma::vec &n                          - Vector of number densities
  double compute_n(const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons=0.0);
  double compute_n(const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons=0.0);
  double compute_n(const std::vector<arma::vec> &n_vl_molecule, double n_electrons=0.0);
  double compute_n(const arma::vec &n);

  // Calculation of the number density of molecules
  arma::vec compute_n_molecule(const std::vector<arma::vec> &n_vl_molecule);

  // Number density of the mixture
  arma::vec compute_density_array(const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons=0.0);
  arma::vec compute_density_array(const std::vector<arma::vec> &n_vl_molecule);
  double compute_density(const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons=0.0);
  double compute_density(const std::vector<arma::vec> &n_vl_molecule, double n_electrons=0.0);
  double compute_density(const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons=0.0);
  double compute_density(const arma::vec &n);

  // Compute the pressure of the mixture
  double compute_pressure(double T, const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons=0.0);

  // Calculation of the specific heat of translational degrees of freedom
  double c_tr(const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons=0.0);
  double c_tr(const std::vector<arma::vec> &n_vl_molecule, double n_electrons=0.0);
  double c_tr(const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons=0.0);
  double c_tr(const arma::vec &n);

  // Calculation of the specific heat of rotational degrees of freedom
  double c_rot(double T, const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons=0.0);
  double c_rot(double T, const std::vector<arma::vec> &n_vl_molecule, double n_electrons=0.0);
  double c_rot(double T, const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons=0.0);
  double c_rot(double T, const arma::vec &n);

  void compute_transport_coefficients(double T, 
          const std::vector<arma::vec> &n_vl_molecule, 
          const arma::vec &n_atom, 
          double n_electrons, 
          kappa::models_omega model=kappa::models_omega::model_omega_esa, 
          double perturbation=1e-9);

  void compute_transport_coefficients(double T, 
          const std::vector<arma::vec> &n_vl_molecule, 
          const arma::vec &n_atom, 
          kappa::models_omega model=kappa::models_omega::model_omega_esa, 
          double perturbation=1e-9);

  // void compute_transport_coefficients(double T, const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons, kappa::models_omega model=kappa::models_omega::model_omega_esa, double perturbation=1e-9);	
  
  // void compute_transport_coefficients(double T, const arma::vec &n_molecule, const arma::vec &n_atom, kappa::models_omega model=kappa::models_omega::model_omega_esa, double perturbation=1e-9);

  void compute_transport_coefficients(double T, 
          const std::vector<arma::vec> &n_vl_molecule, 
          double n_electrons, 
          kappa::models_omega model=kappa::models_omega::model_omega_esa, 
          double perturbation=1e-9);	

  void compute_transport_coefficients(double T,
          const std::vector<arma::vec> &n_vl_molecule, 
          kappa::models_omega model=kappa::models_omega::model_omega_esa, 
          double perturbation=1e-9);

  void compute_transport_coefficients(double T,	
          const arma::vec &n, 
          kappa::models_omega model=kappa::models_omega::model_omega_esa, 
          double perturbation=1e-9);

  double get_thermal_conductivity();
  double get_shear_viscosity();
  double get_bulk_viscosity();
  arma::vec get_thermodiffusion();

  // Returns the matrix of diffusion coefficients 
  // (for calculation it is necessary to call the function compute_transport_coefficients)
  arma::mat get_diffusion();
  arma::mat get_lite_diffusion();
  arma::vec get_binary_diffusion();
  // arma::mat binary_diffusion(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);

 protected:
        
  int inter_index(int i, int j);
  std::vector<std::string> split_string(std::string input);

  arma::vec molecule_charges;
  arma::vec atom_charges;
  arma::vec molecule_charges_sq;
  arma::vec atom_charges_sq; // squared charges

  arma::mat omega_11;
  arma::mat omega_12;
  arma::mat omega_13;
  arma::mat omega_22;
  arma::mat rot_rel_times;
  arma::vec c_rot_arr;
  arma::vec c_rot_rigid_rot_arr;
      
  arma::mat shear_viscosity_LHS;
  arma::mat thermal_conductivity_LHS;
  arma::mat thermal_conductivity_rigid_rot_LHS;

  arma::mat diffusion_LHS;
  arma::mat diffusion_rigid_rot_LHS;

  arma::mat bulk_viscosity_LHS;
  arma::mat bulk_viscosity_rigid_rot_LHS;

  arma::vec thermal_conductivity_RHS;
  arma::vec thermal_conductivity_rigid_rot_RHS;

  arma::vec diffusion_RHS;
  arma::vec diffusion_rigid_rot_RHS;

  arma::vec shear_viscosity_RHS;
  arma::vec bulk_viscosity_RHS;
  arma::vec bulk_viscosity_rigid_rot_RHS;
  arma::vec thermal_conductivity_rigid_rot_coeffs;
  arma::vec thermal_conductivity_coeffs;

  arma::vec diffusion_coeffs;
  arma::vec diffusion_rigid_rot_coeffs;

  arma::vec shear_viscosity_coeffs;
  arma::vec bulk_viscosity_coeffs;
  arma::vec bulk_viscosity_rigid_rot_coeffs;

  // if our mixture has no atoms, we pass this as a dummy parameter to functions which expect numeric density of atomic species
  arma::vec empty_n_atom; 

  std::vector<arma::vec> empty_n_vl_molecule;

  arma::vec this_n_atom;
  std::vector<arma::vec> this_n_vl_mol;
  arma::vec this_n_molecules;
  
  double this_total_n;
  double this_total_dens;
  double this_ctr;
  double this_crot;
  double this_n_electrons;

  std::vector<kappa::Molecule> molecules;
  std::vector<kappa::Atom> atoms;
  kappa::Particle electron;

  int num_molecules;
  int num_atoms;
  int n_vibr_levels_total;
  int n_particles;
  bool cache_on;
  bool all_rigid_rotators = true;
  bool is_ionized;

  // this allows to easily calculate which elements of state-to-state matrix we need to fill, given the molecule indices and the vibrational levels
  // vl_offset[0] = 0; 
  // vl_offset[1] = molecules[0].num_vibr_levels[0]; 
  // vl_offset[2] = molecules[0].num_vibr_levels[0] + molecules[1].num_vibr_levels[0], etc.
  std::vector<int> vl_offset; 
      
  std::vector<kappa::Interaction> interactions;
  std::map<std::string, int> molecule_name_map;
  std::map<std::string, int> atom_name_map;

  void init_matrices(const std::string &particles_filename);
  void add_interactions(const std::string &filename);

  void check_n_vl_molecule(const std::vector<arma::vec> &n_vl_molecule); // check sizes of arrays and test values for non-negativity
  void check_n_molecule(const arma::vec &n_molecule);                    // check sizes of arrays and test values for non-negativity
  void check_n_atom(const arma::vec &n_atom);
  void check_n(const arma::vec &n);

  void compute_omega11(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
  void compute_omega12(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
  void compute_omega13(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
  void compute_omega22(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);

  // Particle + e- interactions
  void compute_omega11(double T, double debye_length, kappa::models_omega model=kappa::models_omega::model_omega_esa);
  void compute_omega12(double T, double debye_length, kappa::models_omega model=kappa::models_omega::model_omega_esa);
  void compute_omega13(double T, double debye_length, kappa::models_omega model=kappa::models_omega::model_omega_esa);
  void compute_omega22(double T, double debye_length, kappa::models_omega model=kappa::models_omega::model_omega_esa);

  void compute_c_rot_rigid_rot(double T);
  void compute_c_rot(double T); // compute c_rot arrays
  void compute_full_crot_rigid_rot(double T);
  void compute_full_crot(double T); // compute c_rot of mixture
  void compute_rot_rel_times(double T, double n, kappa::models_omega model=kappa::models_omega::model_omega_esa);
  void inplace_compute_n_molecule(const std::vector<arma::vec> &n_vl_molecule);

  double th_cond, sh_visc, b_visc;
  
  arma::vec th_diff;
  arma::mat diff;
  arma::mat lite_diff;
  arma::vec binary_diff;

  const arma::mat &compute_bulk_viscosity_LHS(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
  const arma::vec &compute_bulk_viscosity_RHS(double T);
  const arma::vec &compute_bulk_viscosity_coeffs(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);

  const arma::mat &compute_bulk_viscosity_rigid_rot_LHS(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
  const arma::vec &compute_bulk_viscosity_rigid_rot_RHS(double T);
  const arma::vec &compute_bulk_viscosity_rigid_rot_coeffs(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);

  const arma::mat &compute_thermal_conductivity_LHS(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
  const arma::vec &compute_thermal_conductivity_RHS(double T);
  const arma::vec &compute_thermal_conductivity_coeffs(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);

  const arma::mat &compute_thermal_conductivity_rigid_rot_LHS(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
  const arma::vec &compute_thermal_conductivity_rigid_rot_RHS(double T);
  const arma::vec &compute_thermal_conductivity_rigid_rot_coeffs(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);

  // const arma::mat &compute_diffusion_LHS(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
  // const arma::vec &compute_diffusion_LHS(double T);
  // const arma::vec &compute_diffusion_RHS(double T);
  // const arma::vec &compute_diffusion_RHS(double T, int b, int n);
  // const arma::vec &compute_diffusion_coeffs(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
  const arma::mat &compute_diffusion_LHS(double T);
  const arma::mat &lite_compute_diffusion_LHS(double T);
  const arma::mat &lite_compute_diffusion_LHS(double T, int b, int n, std::string ParticleType);
  const arma::vec &compute_diffusion_RHS(double T, int b, int n);
  const arma::vec &lite_compute_diffusion_RHS(double T, int b, int n, std::string ParticleType);
  const arma::vec &compute_diffusion_coeffs(double T, int b, int n);
  const arma::vec &lite_compute_diffusion_coeffs(double T, int b, int n, std::string ParticleType);
  // const arma::mat &compute_diffusion_rigid_rot_LHS(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
  // const arma::vec &compute_diffusion_rigid_rot_RHS(double T);
  // const arma::vec &compute_diffusion_rigid_rot_coeffs(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);

  const arma::mat &compute_shear_viscosity_LHS(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
  const arma::vec &compute_shear_viscosity_RHS(double T);
  const arma::vec &compute_shear_viscosity_coeffs(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);

  double bulk_viscosity(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
  double thermal_conductivity(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);    
  double shear_viscosity(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);

  void diffusion(double T);
  void lite_diffusion(double T);
  
  void thermodiffusion(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);
  void binary_diffusion(double T, kappa::models_omega model=kappa::models_omega::model_omega_esa);

}; // class Mixture
} // namespace kappa
#endif /* kappa_mixture_hpp */
