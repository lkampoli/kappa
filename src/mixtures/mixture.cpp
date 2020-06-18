/*!
    \file mixture.cpp
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <set>
#include "mixture.hpp"
#include "kappa.hpp"

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  kappa::Mixture::Mixture(const std::vector<kappa::Molecule> &i_molecules, const std::vector<kappa::Atom> &i_atoms,
                          const std::string &interactions_filename, const std::string &particles_filename)

    : kappa::Approximation()
{

      molecules = i_molecules;
      atoms = i_atoms;
      num_molecules = molecules.size();
      num_atoms = atoms.size();
      cache_on = false;
      n_vibr_levels_total = 0;

      // TODO: read the mixture and determine if there are e-
      is_ionized = false; // avoid e- loading

      int i;

      for (i = 0; i < num_molecules; i++) {
        if (!molecules[i].rigid_rotator) {
          all_rigid_rotators = false;
        }
        molecule_name_map[molecules[i].name] = i;
        n_vibr_levels_total += molecules[i].num_vibr_levels[0];
      }

      for (i = 0; i < num_atoms; i++) {
        atom_name_map[atoms[i].name] = i;
      }

      // initialization of matrices for transport coefficients computation
      init_matrices(particles_filename);

      // according to loaded particles, creates interactions
      add_interactions(interactions_filename);
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  kappa::Mixture::Mixture(const std::vector<kappa::Molecule> &i_molecules, const std::string &interactions_filename,
                          const std::string &particles_filename)

    : kappa::Approximation()
{

      molecules = i_molecules;
      num_molecules = molecules.size();
      num_atoms = 0;
      cache_on = false;
      n_vibr_levels_total = 0;

      for (int i = 0; i < num_molecules; i++) {
        if (!molecules[i].rigid_rotator) {
          all_rigid_rotators = false;
       }
       molecule_name_map[molecules[i].name] = i;
       n_vibr_levels_total += molecules[i].num_vibr_levels[0];
      }

      init_matrices(particles_filename);
      add_interactions(interactions_filename);
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  kappa::Mixture::Mixture(const std::vector<kappa::Atom> &i_atoms, const std::string &interactions_filename, const std::string &particles_filename)

    : kappa::Approximation() {

      atoms = i_atoms;
      num_molecules = 0;
      num_atoms = atoms.size();
      cache_on = false;
      n_vibr_levels_total = 0;
      all_rigid_rotators = true;

      for (int i = 0; i < num_atoms; i++) {
        atom_name_map[atoms[i].name] = i;
      }

      init_matrices(particles_filename);
      add_interactions(interactions_filename);
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  kappa::Mixture::Mixture(const kappa::Molecule &molecule, const kappa::Atom &atom, const std::string &interactions_filename, const std::string &particles_filename)

    : kappa::Approximation() {

      molecules.push_back(molecule);
      atoms.push_back(atom);
      num_molecules = 1;
      num_atoms = 1;
      cache_on = false;
      n_vibr_levels_total = 0;

      int i;
      for (i = 0; i < num_molecules; i++) {
        if (!molecules[i].rigid_rotator) {
          all_rigid_rotators = false;
       }
       molecule_name_map[molecules[i].name] = i;
       n_vibr_levels_total += molecules[i].num_vibr_levels[0];
      }

      for (i = 0; i < num_atoms; i++) {
        atom_name_map[atoms[i].name] = i;
      }

      init_matrices(particles_filename);
      add_interactions(interactions_filename);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  kappa::Mixture::Mixture(const kappa::Molecule &molecule, const std::string &interactions_filename, const std::string &particles_filename)

    : kappa::Approximation() {

      molecules.push_back(molecule);
      num_molecules = 1;
      num_atoms = 0;
      cache_on = false;
      n_vibr_levels_total = 0;

      for (int i = 0; i < num_molecules; i++) {
        if (!molecules[i].rigid_rotator) {
          all_rigid_rotators = false;
       }
       molecule_name_map[molecules[i].name] = i;
       n_vibr_levels_total += molecules[i].num_vibr_levels[0];
      }

      init_matrices(particles_filename);
      add_interactions(interactions_filename);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  kappa::Mixture::Mixture(const kappa::Atom &atom, const std::string &interactions_filename, const std::string &particles_filename)

    : kappa::Approximation() {

      atoms.push_back(atom);
      num_molecules = 0;
      num_atoms = 1;
      cache_on = false;
      n_vibr_levels_total = 0;
      all_rigid_rotators = true;

      for (int i = 0; i < num_atoms; i++) {
        atom_name_map[atoms[i].name] = i;
      }

      init_matrices(particles_filename);
      add_interactions(interactions_filename);
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  kappa::Mixture::Mixture(const std::string particle_names, const std::string &interactions_filename, const std::string &particles_filename, bool anharmonic, bool rigid_rotators)

    : kappa::Approximation() {

      std::vector<std::string> split_names = split_string(particle_names);
      all_rigid_rotators = rigid_rotators;

      for (auto name : split_names) {
        std::cout << " Particles' name: " << name << std::endl;
        if (name != "e-") { // if there is an electron, there should be charged particles and the electron data will be loaded anyway
          try {
            kappa::Molecule mol(name, anharmonic, rigid_rotators, particles_filename);
            molecules.push_back(mol);
          } catch (const DataNotFoundException& e) {
            kappa::Atom atom(name, particles_filename);
            atoms.push_back(atom);
          }
        } else {
          is_ionized = true;
        }
      }
      num_molecules = molecules.size();
      num_atoms = atoms.size();
      cache_on = false;
      n_vibr_levels_total = 0;

      int i;
      for (i=0; i<num_molecules; i++) {
        molecule_name_map[molecules[i].name] = i;
        n_vibr_levels_total += molecules[i].num_vibr_levels[0];
      }

      for (i=0; i<num_atoms; i++) {
        atom_name_map[atoms[i].name] = i;
      }

      init_matrices(particles_filename);
      add_interactions(interactions_filename);
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // initialization of the matrix for transport coefficients
  void kappa::Mixture::init_matrices(const std::string &particles_filename) {

    int i;
    empty_n_atom.zeros(1);
    arma::vec n_vl_zero(1);
    n_vl_zero.zeros(1);
    empty_n_vl_molecule.push_back(n_vl_zero);

    this_n_molecules.zeros(num_molecules);

    n_particles = num_molecules + num_atoms;

    // initialize the offset
    if (num_molecules > 0) {
      vl_offset.push_back(0);
      for (i=0; i<num_molecules; i++) {
        vl_offset.push_back(molecules[i].num_vibr_levels[0] + vl_offset[i]);
      }

      molecule_charges.zeros(num_molecules);
      molecule_charges_sq.zeros(num_molecules);

      for (i=0; i<num_molecules; i++) {
        molecule_charges[i] = molecules[i].charge * K_CONST_ELEMENTARY_CHARGE;
        if (molecules[i].charge != 0) {
          is_ionized = true;
        }
      }
      molecule_charges_sq = molecule_charges % molecule_charges;
    } else {
      molecule_charges.zeros(1);
      molecule_charges_sq.zeros(1);
    }

    if (num_atoms > 0) {
      atom_charges.zeros(num_atoms);
      atom_charges_sq.zeros(num_atoms);
      for (i=0; i<num_atoms; i++) {
        atom_charges[i] = atoms[i].charge * K_CONST_ELEMENTARY_CHARGE;
        if (atoms[i].charge != 0) {
          is_ionized = true;
        }
     }
     atom_charges_sq = atom_charges % atom_charges;
    } else {
      atom_charges.zeros(1);
      atom_charges_sq.zeros(1);
    }

    // initialize thermal diffusion coefficient vector
    if (is_ionized) {
      electron = kappa::Particle("e-", particles_filename);
      n_particles += 1;
      th_diff.zeros(n_vibr_levels_total + num_atoms + 1);
    } else {
      th_diff.zeros(n_vibr_levels_total + num_atoms);
    }

    shear_viscosity_RHS.zeros(n_particles); // n_particles = num_molecules + num_atoms
    shear_viscosity_LHS.zeros(n_particles, n_particles);
    shear_viscosity_coeffs.zeros(n_particles);

    thermal_conductivity_RHS.zeros(3 * n_vibr_levels_total + 2 * num_atoms); // fix for ionized mixtures!
    thermal_conductivity_LHS.zeros(3 * n_vibr_levels_total + 2 * num_atoms, 3 * n_vibr_levels_total + 2 * num_atoms);
    thermal_conductivity_coeffs.zeros(3 * n_vibr_levels_total + 2 * num_atoms);

    thermal_conductivity_rigid_rot_RHS.zeros(3 * num_molecules + 2 * num_atoms);
    thermal_conductivity_rigid_rot_LHS.zeros(3 * num_molecules + 2 * num_atoms, 3 * num_molecules + 2 * num_atoms);
    thermal_conductivity_rigid_rot_coeffs.zeros(3 * num_molecules + 2 * num_atoms); // a_ci_00, a_ci_10, a_ci_01

    bulk_viscosity_rigid_rot_RHS.zeros(2 * num_molecules + num_atoms);
    bulk_viscosity_rigid_rot_LHS.zeros(2 * num_molecules + num_atoms, 2 * num_molecules + num_atoms);
    bulk_viscosity_rigid_rot_coeffs.zeros(2 * num_molecules + num_atoms);

    bulk_viscosity_RHS.zeros(2 * n_vibr_levels_total + num_atoms);
    bulk_viscosity_LHS.zeros(2 * n_vibr_levels_total + num_atoms, 2 * n_vibr_levels_total + num_atoms);
    bulk_viscosity_coeffs.zeros(2 * n_vibr_levels_total + num_atoms);

    diffusion_LHS.zeros(n_vibr_levels_total + num_atoms, n_vibr_levels_total + num_atoms);
    diffusion_RHS.zeros(n_vibr_levels_total + num_atoms);
    diffusion_coeffs.zeros(n_vibr_levels_total + num_atoms);
    diff.zeros(n_vibr_levels_total + num_atoms, n_vibr_levels_total + num_atoms);

    omega_11.zeros(n_particles, n_particles);
    omega_12.zeros(n_particles, n_particles);
    omega_13.zeros(n_particles, n_particles);
    omega_22.zeros(n_particles, n_particles);
    c_rot_arr.zeros(n_vibr_levels_total);
    c_rot_rigid_rot_arr.zeros(num_molecules);
    rot_rel_times.zeros(num_molecules, num_molecules + num_atoms);

    this_total_n = 0.0;
    this_n_electrons = 0.0;
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Defines interactions between molecules, atoms and ions
  void kappa::Mixture::add_interactions(const std::string &filename) {

    int i, j;

    for (i = 0; i < num_molecules; i++) {
      // molecule-molecule interaction
      for (j = i; j < num_molecules; j++) {
        interactions.push_back(kappa::Interaction(molecules[i], molecules[j], filename));
      }
      // molecule-atom interaction
      for (j = 0; j < num_atoms; j++) {
        interactions.push_back(kappa::Interaction(molecules[i], atoms[j], filename));
       }
       // molecule-ions interaction
       if (is_ionized) {
         interactions.push_back(kappa::Interaction(molecules[i], electron, filename));
       }
    }

    for (i = 0; i < num_atoms; i++) {
      // atom-atom interaction
      for (j = i; j < num_atoms; j++) {
        interactions.push_back(kappa::Interaction(atoms[i], atoms[j], filename));
      }
      // atom-electron interaction
      if (is_ionized) {
        interactions.push_back(kappa::Interaction(atoms[i], electron, filename));
      }
    }

    // electron-electron interaction
    if (is_ionized) {
      interactions.push_back(kappa::Interaction(electron, electron, filename));
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // returns the name of the particles in the mixture
  std::string kappa::Mixture::get_names() {

    std::string res="";
    int i;

    for (i=0; i<num_molecules; i++) {
      res += molecules[i].name + " ";
    }
    for (i=0; i<num_atoms-1; i++) {
      res += atoms[i].name + " ";
    }
    if (num_atoms > 0) {
      res += atoms[num_atoms-1].name;
    }
    if (is_ionized) {
      res += " e-";
    }
    return res;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int kappa::Mixture::get_n_particles() {
    return n_particles;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int kappa::Mixture::get_n_vibr_levels() {
    return n_vibr_levels_total;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector<int> kappa::Mixture::get_n_vibr_levels_array() {

    std::vector<int> res;
    for (int i=0; i<num_molecules; i++) {
      res.push_back(molecules[i].num_vibr_levels[0]);
    }
    return res;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // molar-mass fraction conversion
  arma::vec kappa::Mixture::convert_molar_frac_to_mass(const arma::vec &x) {

    #ifdef KAPPA_STRICT_CHECKS
    check_n(x);
    #endif

    int i;
    double rho = 0.;

    // dimensioning and initialization
    arma::vec res = arma::zeros(x.size());

    for (i=0; i<num_molecules; i++) {
      rho += x[i] * molecules[i].mass;
    }

    for (i=0; i<num_atoms; i++) {
      rho += x[num_molecules + i] * atoms[i].mass;
    }

    if (is_ionized) {
      rho += x[num_molecules + num_atoms + i] * electron.mass;
    }

    for (i=0; i<num_molecules; i++) {
      res[i] = x[i] * molecules[i].mass / rho;
    }

    for (i=0; i<num_atoms; i++) {
      res[num_molecules + i] = x[num_molecules + i] * atoms[i].mass / rho;
    }

    if (is_ionized) {
      res[num_molecules + num_atoms + i] = x[num_molecules + num_atoms + i] * electron.mass / rho;
    }

    //std::cout << res << std::endl;
    return res;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // mass-molar fraction conversion
  arma::vec kappa::Mixture::convert_mass_frac_to_molar(const arma::vec &y) {

    #ifdef KAPPA_STRICT_CHECKS
    check_n(y);
    #endif

    int i;
    double rho = 0.0;

    // dimensioning and initialization
    arma::vec res = arma::zeros(y.size());;

    for (i=0; i<num_molecules; i++) {
      rho += y[i] / molecules[i].mass;
    }

    for (i=0; i<num_atoms; i++) {
      rho += y[num_molecules + i] / atoms[i].mass;
    }

    if (is_ionized) {
      rho += y[num_molecules + num_atoms + i] / electron.mass;
    }

    for (i=0; i<num_molecules; i++) {
      res[i] = y[i] / rho / molecules[i].mass;
    }

    for (i=0; i<num_atoms; i++) {
      res[num_molecules + i] = y[num_molecules + i] / rho / atoms[i].mass;
    }

    if (is_ionized) {
      res[num_molecules + num_atoms + i] = y[num_molecules + num_atoms + i] / rho / electron.mass;
    }

    return res;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // return interaction index for particles with indices i and j
  int kappa::Mixture::inter_index(int i, int j) {

    if (i <= j) {
      return (n_particles) * i + j - i * (i + 1) / 2;
    } else {
      return inter_index(j, i);
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector<std::string> kappa::Mixture::split_string(std::string input) {

    std::vector<std::string> result;
    std::string holder;
    for (auto c : input) {
      if ((c == ',') && (holder.length() > 0)) {
        result.push_back(holder);
        holder = "";
      } else if (c != ' ') {
        holder += c;
      }
    }
    if (holder.length() > 0) {
      result.push_back(holder);
    }
    return result;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute constant volume specific heat for rigid rotator
  void kappa::Mixture::compute_c_rot_rigid_rot(double T) {

    for (int i = 0; i < num_molecules; i++) {
      c_rot_rigid_rot_arr[i] = c_rot(T, molecules[i], 0, 0);
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute constant volume specific heat
  void kappa::Mixture::compute_c_rot(double T) {

    int j = 0;
    for (int i = 0; i < num_molecules; i++) {
      for (int k = 0; k < molecules[i].num_vibr_levels[0]; k++) {
        c_rot_arr[j] = c_rot(T, molecules[i], k, 0);
        j++;
      }
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute constant volume specific heat of the mixture for rigid rotator
  void kappa::Mixture::compute_full_crot_rigid_rot(double T) {

    double res = 0;
    for (int i = 0; i < num_molecules; i++) {
      res += c_rot_rigid_rot_arr[i] * this_n_molecules[i] * molecules[i].mass;
    }
    this_crot = res / this_total_dens;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute constant volume specific heat of the mixture
  void kappa::Mixture::compute_full_crot(double T) {

    double res = 0;
    int k=0;
    for (int i=0; i<num_molecules; i++) {
      for (int j=0; j<molecules[i].num_vibr_levels[0]; j++) {
        res += c_rot_arr[k] * this_n_vl_mol[i][j] * molecules[i].mass;
        k++;
      }
    }
    this_crot = res / this_total_dens;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // rotational relaxation time (Parker's model)
  void kappa::Mixture::compute_rot_rel_times(double T, double n, kappa::models_omega model) {

    int i, j;
    for (i = 0; i < num_molecules; i++) {
      for (j = 0; j < num_molecules; j++) {
        rot_rel_times.at(i, j) = rot_relaxation_time_parker(T, n, molecules[i], interactions[inter_index(i, j)], model);
      }
      for (j = 0; j < num_atoms; j++) {
        rot_rel_times.at(i, num_molecules + j) = rot_relaxation_time_parker(T, n, molecules[i], interactions[inter_index(i, num_molecules + j)], model);
      }
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // sanity check for molecular number density vector of arma vector
  void kappa::Mixture::check_n_vl_molecule(const std::vector<arma::vec> &n_vl_molecule) {

    std::string error_string;
    if (n_vl_molecule.size() != num_molecules) {
      error_string = "Incorrect size of array of vibrational level populations of molecules, wrong number of molecules: ";
      error_string += "passed array is of size " + std::to_string(n_vl_molecule.size()) + ", should be of size " + std::to_string(num_molecules);
      throw kappa::IncorrectValueException(error_string.c_str());
    }

    for (int i = 0; i < num_molecules; i++) {
      if (molecules[i].num_vibr_levels[0] != n_vl_molecule[i].n_elem) {
        error_string = "Incorrect size of array of vibrational level populations of molecule " + molecules[i].name + ": ";
        error_string += "passed array is of size " + std::to_string(n_vl_molecule[i].n_elem) + ", should be of size " + std::to_string(molecules[i].num_vibr_levels[0]);
        throw kappa::IncorrectValueException(error_string.c_str());
      }
      for (int k = 0; k < molecules[i].num_vibr_levels[0]; k++) {
        if (n_vl_molecule[i][k] < 0) {
          error_string = "Array of vibrational level populations of molecule " + molecules[i].name + " contains at least one negative value: ";
          error_string += "n_vl_molecule[" + std::to_string(i) + "][" + std::to_string(k) + "]=" + std::to_string(n_vl_molecule[i][k]);
          throw kappa::IncorrectValueException(error_string.c_str());
        }
      }
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // sanity check for molecular number density arma vector
  void kappa::Mixture::check_n_molecule(const arma::vec &n_molecule) {

    std::string error_string;
    if (n_molecule.n_elem != num_molecules) {
      error_string = "Incorrect size of array of number densities of molecules, wrong number of molecules: ";
      error_string += "passed array is of size " + std::to_string(n_molecule.n_elem) + ", should be of size " + std::to_string(num_molecules);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    for (int i = 0; i < num_molecules; i++) {
      if (n_molecule[i] < 0) {
        error_string = "Number density of molecule " + molecules[i].name + " is negative: ";
        error_string += "n_molecule[" + std::to_string(i) + "]=" + std::to_string(n_molecule[i]);
        throw kappa::IncorrectValueException(error_string.c_str());
      }
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // sanity check for atomic number density arma vector
  void kappa::Mixture::check_n_atom(const arma::vec &n_atom) {

    std::string error_string;
    if (n_atom.n_elem != num_atoms) {
      error_string = "Incorrect size of array of number densities of atoms, wrong number of atoms: ";
      error_string += "passed array is of size " + std::to_string(n_atom.n_elem) + ", should be of size " + std::to_string(num_atoms);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    for (int i=0; i < num_atoms; i++) {
      if (n_atom[i] < 0) {
        error_string = "Number density of atom " + atoms[i].name + " is negative: ";
        error_string += "n_atom[" + std::to_string(i) + "]=" + std::to_string(n_atom[i]);
        throw kappa::IncorrectValueException(error_string.c_str());
      }
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // sanity check for a general number density arma vector
  void kappa::Mixture::check_n(const arma::vec &n) {

    std::string error_string;
    int i;
    if (n.n_elem != n_particles) {
      error_string = "Incorrect size of array of number densities of particles: ";
      error_string += "passed array is of size " + std::to_string(n.n_elem) + ", should be of size " + std::to_string(n_particles);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    for (i = 0; i < num_molecules; i++) {
      if (n[i] < 0) {
        error_string = "Number density of molecule " + molecules[i].name + " is negative: ";
        error_string += "n[" + std::to_string(i) + "]=" + std::to_string(n[i]);
        throw kappa::IncorrectValueException(error_string.c_str());
      }
    }
    for (i = 0; i < num_atoms; i++) {
      if (n[num_molecules + i] < 0) {
        error_string = "Number density of atom " + atoms[i].name + " is negative: ";
        error_string += "n[" + std::to_string(i) + "]=" + std::to_string(n[num_molecules + i]);
        throw kappa::IncorrectValueException(error_string.c_str());
      }
    }
    if (is_ionized) {
      if (n[n.n_elem - 1] < 0.0) {
        error_string = "Number density of electrons is negative: ";
        error_string += "n[" + std::to_string(n.n_elem - 1) + "]=" + std::to_string(n[n.n_elem - 1]);
        throw kappa::IncorrectValueException(error_string.c_str());
      }
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Ω_11
  void kappa::Mixture::compute_omega11(double T, kappa::models_omega model) {

    int i1, i2;
    int interindex;

    // molecule + molecule collisions
    for (i1=0; i1<num_molecules; i1++) {
      for (i2=i1; i2<num_molecules; i2++) {
        interindex = inter_index(i1, i2);
        omega_11.at(i1, i2) = omega_integral(T, interactions[interindex], 1, 1, model, true);
        omega_11.at(i2, i1) = omega_11.at(i1, i2);
      }
    }

    // molecule + atom collisions
    for (i1=0; i1<num_molecules; i1++) {
      for (i2=0; i2<num_atoms; i2++) {
        interindex = inter_index(i1, num_molecules+i2);
        omega_11.at(i1, num_molecules+i2) = omega_integral(T, interactions[interindex], 1, 1, model, true);
        omega_11.at(num_molecules+i2, i1) = omega_11.at(i1, num_molecules+i2);
      }
    }

    // atom + atom collisions
    for (i1=0; i1<num_atoms; i1++) {
      for (i2=i1; i2<num_atoms; i2++) {
        interindex = inter_index(num_molecules+i1, num_molecules+i2);
        omega_11.at(num_molecules+i1, num_molecules+i2) = omega_integral(T, interactions[interindex], 1, 1, model, true);
        omega_11.at(num_molecules+i2, num_molecules+i1) = omega_11.at(num_molecules+i1, num_molecules+i2);
      }
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void kappa::Mixture::compute_omega12(double T, kappa::models_omega model) {

    int i1, i2;
    int interindex;

    // molecule + molecule collisions
    for (i1=0; i1<num_molecules; i1++) {
      for (i2=i1; i2<num_molecules; i2++) {
        interindex = inter_index(i1, i2);
        omega_12.at(i1, i2) = omega_integral(T, interactions[interindex], 1, 2, model, true);
        omega_12.at(i2, i1) = omega_12.at(i1, i2);
      }
    }

    // molecule + atom collisions
    for (i1=0; i1<num_molecules; i1++) {
      for (i2=0; i2<num_atoms; i2++) {
        interindex = inter_index(i1, num_molecules+i2);
        omega_12.at(i1, num_molecules+i2) = omega_integral(T, interactions[interindex], 1, 2, model, true);
        omega_12.at(num_molecules+i2, i1) = omega_12.at(i1, num_molecules+i2);
      }
    }

    // atom + atom collisions
    for (i1=0; i1<num_atoms; i1++) {
      for (i2=i1; i2<num_atoms; i2++) {
        interindex = inter_index(num_molecules+i1, num_molecules+i2);
        omega_12.at(num_molecules+i1, num_molecules+i2) = omega_integral(T, interactions[interindex], 1, 2, model, true);
        omega_12.at(num_molecules+i2, num_molecules+i1) = omega_12.at(num_molecules+i1, num_molecules+i2);
      }
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void kappa::Mixture::compute_omega13(double T, kappa::models_omega model) {

    int i1, i2;
    int interindex;

    // molecule + molecule collisions
    for (i1=0; i1<num_molecules; i1++) {
      for (i2=i1; i2<num_molecules; i2++) {
        interindex = inter_index(i1, i2);
        omega_13.at(i1, i2) = omega_integral(T, interactions[interindex], 1, 3, model, true);
        omega_13.at(i2, i1) = omega_13.at(i1, i2);
      }
    }

    // molecule + atom collisions
    for (i1=0; i1<num_molecules; i1++) {
      for (i2=0; i2<num_atoms; i2++) {
        interindex = inter_index(i1, num_molecules + i2);
        omega_13.at(i1, num_molecules+i2) = omega_integral(T, interactions[interindex], 1, 3, model, true);
        omega_13.at(num_molecules+i2, i1) = omega_13.at(i1, num_molecules+i2);
      }
    }

    // atom + atom collisions
    for (i1=0; i1<num_atoms; i1++) {
      for (i2=i1; i2<num_atoms; i2++) {
        interindex = inter_index(num_molecules+i1, num_molecules+i2);
        omega_13.at(num_molecules+i1, num_molecules+i2) = omega_integral(T, interactions[interindex], 1, 3, model, true);
        omega_13.at(num_molecules+i2, num_molecules+i1) = omega_13.at(num_molecules+i1, num_molecules+i2);
      }
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void kappa::Mixture::compute_omega22(double T, kappa::models_omega model) {

    int i1, i2;
    int interindex;

    for (i1=0; i1<num_molecules; i1++) { // molecule + molecule collisions
      for (i2=i1; i2<num_molecules; i2++) {
        interindex = inter_index(i1, i2);
        omega_22.at(i1, i2) = omega_integral(T, interactions[interindex], 2, 2, model, true);
        omega_22.at(i2, i1) = omega_22.at(i1, i2);
      }
    }

    for (i1=0; i1<num_molecules; i1++) { // molecule + atom collisions
      for (i2=0; i2<num_atoms; i2++) {
        interindex = inter_index(i1, num_molecules + i2);
        omega_22.at(i1, num_molecules+i2) = omega_integral(T, interactions[interindex], 2, 2, model, true);
        omega_22.at(num_molecules+i2, i1) = omega_22.at(i1, num_molecules+i2);
      }
    }

    for (i1=0; i1<num_atoms; i1++) { // atom + atom collisions
      for (i2=i1; i2<num_atoms; i2++) {
        interindex = inter_index(num_molecules+i1, num_molecules+i2);
        omega_22.at(num_molecules+i1, num_molecules+i2) = omega_integral(T, interactions[interindex], 2, 2, model, true);
        omega_22.at(num_molecules+i2, num_molecules+i1) = omega_22.at(num_molecules+i1, num_molecules+i2);
      }
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void kappa::Mixture::compute_omega11(double T, double debye_length, kappa::models_omega model) {

    int i;

    for (i=0; i<num_molecules; i++) {
      omega_11.at(i, num_molecules+num_atoms) = omega_integral(T, interactions[inter_index(i, num_molecules+num_atoms)], 1, 1, debye_length, model, true);
      omega_11.at(num_molecules+num_atoms, i) = omega_11.at(i, num_molecules+num_atoms);
    }

    for (i=0; i<num_atoms; i++) {
      omega_11.at(num_molecules+i, num_molecules+num_atoms) = omega_integral(T, interactions[inter_index(num_molecules + i, num_molecules+num_atoms)], 1, 1, debye_length, model, true);
      omega_11.at(num_molecules+num_atoms, num_molecules+i) = omega_11.at(num_molecules+i, num_molecules+num_atoms);
    }
    omega_11.at(num_molecules+num_atoms, num_molecules+num_atoms) = omega_integral(T, interactions[inter_index(num_molecules+num_atoms, num_molecules+num_atoms)], 1, 1, debye_length, model, true);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void kappa::Mixture::compute_omega12(double T, double debye_length, kappa::models_omega model) {

    int i;

    for (i=0; i<num_molecules; i++) {
      omega_12.at(i, num_molecules+num_atoms) = omega_integral(T, interactions[inter_index(i, num_molecules+num_atoms)], 1, 2, debye_length, model, true);
      omega_12.at(num_molecules+num_atoms, i) = omega_12.at(i, num_molecules+num_atoms);
    }

    for (i=0; i<num_atoms; i++) {
      omega_12.at(num_molecules+i, num_molecules+num_atoms) = omega_integral(T, interactions[inter_index(num_molecules + i, num_molecules+num_atoms)], 1, 2, debye_length, model, true);
      omega_12.at(num_molecules+num_atoms, num_molecules+i) = omega_12.at(num_molecules+i, num_molecules+num_atoms);
    }
    omega_12.at(num_molecules+num_atoms, num_molecules+num_atoms) = omega_integral(T, interactions[inter_index(num_molecules + num_atoms, num_molecules+num_atoms)], 1, 2, debye_length, model, true);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void kappa::Mixture::compute_omega13(double T, double debye_length, kappa::models_omega model) {

    int i;

    for (i=0; i<num_molecules; i++) {
      omega_13.at(i, num_molecules+num_atoms) = omega_integral(T, interactions[inter_index(i, num_molecules+num_atoms)], 1, 3, debye_length, model, true);
      omega_13.at(num_molecules+num_atoms, i) = omega_13.at(i, num_molecules+num_atoms);
    }

    for (i=0; i<num_atoms; i++) {
      omega_13.at(num_molecules+i, num_molecules+num_atoms) = omega_integral(T, interactions[inter_index(num_molecules + i, num_molecules+num_atoms)], 1, 3, debye_length, model, true);
      omega_13.at(num_molecules+num_atoms, num_molecules+i) = omega_13.at(num_molecules+i, num_molecules+num_atoms);
    }
    omega_13.at(num_molecules+num_atoms, num_molecules+num_atoms) = omega_integral(T, interactions[inter_index(num_molecules + num_atoms, num_molecules+num_atoms)], 1, 3, debye_length, model, true);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void kappa::Mixture::compute_omega22(double T, double debye_length, kappa::models_omega model) {

    int i;
    for (i=0; i<num_molecules; i++) {
      omega_22.at(i, num_molecules+num_atoms) = omega_integral(T, interactions[inter_index(i, num_molecules+num_atoms)], 2, 2, debye_length, model, true);
      omega_22.at(num_molecules+num_atoms, i) = omega_22.at(i, num_molecules+num_atoms);
    }

    for (i=0; i<num_atoms; i++) {
      omega_22.at(num_molecules+i, num_molecules+num_atoms) = omega_integral(T, interactions[inter_index(num_molecules + i, num_molecules+num_atoms)], 2, 2, debye_length, model, true);
      omega_22.at(num_molecules+num_atoms, num_molecules+i) = omega_22.at(num_molecules+i, num_molecules+num_atoms);
    }
    omega_22.at(num_molecules+num_atoms, num_molecules+num_atoms) = omega_integral(T, interactions[inter_index(num_molecules + num_atoms, num_molecules+num_atoms)], 2, 2, debye_length, model, true);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void kappa::Mixture::inplace_compute_n_molecule(const std::vector<arma::vec> &n_vl_molecule) {

    for (int i = 0; i < num_molecules; i++) {
      this_n_molecules[i] = arma::sum(n_vl_molecule[i]);
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // molecule + molecule interaction
  kappa::Interaction kappa::Mixture::interaction(const kappa::Molecule &molecule1, const kappa::Molecule &molecule2) {
    return interactions[inter_index(molecule_name_map[molecule1.name], molecule_name_map[molecule2.name])];
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // molecule + atom interaction
  kappa::Interaction kappa::Mixture::interaction(const kappa::Molecule &molecule, const kappa::Atom &atom) {
    return interactions[inter_index(molecule_name_map[molecule.name], num_molecules + atom_name_map[atom.name])];
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // atom + molecule interaction
  kappa::Interaction kappa::Mixture::interaction(const kappa::Atom &atom, const kappa::Molecule &molecule) {
    return interaction(molecule, atom);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // atom + atom interaction
  kappa::Interaction kappa::Mixture::interaction(const kappa::Atom &atom1, const kappa::Atom &atom2) {
    return interactions[inter_index(num_molecules + atom_name_map[atom1.name], num_molecules + atom_name_map[atom2.name])];
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  kappa::Molecule kappa::Mixture::molecule(const std::string &name) {
    return molecules[molecule_name_map[name]];
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  kappa::Atom kappa::Mixture::atom(const std::string &name) {
    return atoms[atom_name_map[name]];
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // computation of Debye length
  double kappa::Mixture::debye_length(double T, const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons) {

    #ifdef KAPPA_STRICT_CHECKS
    std::string error_string;
    if (T <= 0) {
      error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    check_n_molecule(n_molecule);
    check_n_atom(n_atom);
    if (n_electrons < 0) {
      error_string = "Number density of electrons is negative: n_electrons=" + std::to_string(n_electrons);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    return sqrt(K_CONST_E0 * K_CONST_K * T / (arma::dot(n_molecule, molecule_charges_sq) + arma::dot(n_atom, atom_charges_sq) + n_electrons * K_CONST_ELEMENTARY_CHARGE * K_CONST_ELEMENTARY_CHARGE));
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Mixture::debye_length(double T, const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons) {
    return debye_length(T, compute_n_molecule(n_vl_molecule), n_atom, n_electrons);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Mixture::debye_length(double T, const std::vector<arma::vec> &n_vl_molecule, double n_electrons) {

    #ifdef KAPPA_STRICT_CHECKS
    std::string error_string;
    if (T <= 0) {
      error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    check_n_vl_molecule(n_vl_molecule);
    if (n_electrons < 0) {
      error_string = "Number density of electrons is negative: n_electrons=" + std::to_string(n_electrons);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    arma::vec n_mol = compute_n_molecule(n_vl_molecule);
    return sqrt(K_CONST_E0 * K_CONST_K * T / (arma::dot(n_mol, molecule_charges_sq) + n_electrons * K_CONST_ELEMENTARY_CHARGE * K_CONST_ELEMENTARY_CHARGE));
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Mixture::debye_length(double T, const arma::vec &n) {

    #ifdef KAPPA_STRICT_CHECKS
    if (T <= 0) {
      std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    check_n(n);
    #endif

    double dl=0.0;
    if (num_molecules > 0) {
      dl += arma::dot(n.subvec(0, num_molecules-1), molecule_charges_sq);
    }
    if (num_atoms > 0) {
      dl += arma::dot(n.subvec(num_molecules, n_particles-2), atom_charges_sq);
    }
    dl += n[n_particles - 1] * K_CONST_ELEMENTARY_CHARGE * K_CONST_ELEMENTARY_CHARGE;

    return sqrt(K_CONST_E0 * K_CONST_K * T / dl);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute total number density given the vibrational levels of each molecule and atoms (electrons are optional)
  double kappa::Mixture::compute_n(const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons) {

    #ifdef KAPPA_STRICT_CHECKS
    check_n_vl_molecule(n_vl_molecule);
    check_n_atom(n_atom);
    if (n_electrons < 0) {
      std::string error_string = "Number density of electrons is negative: n_electrons=" + std::to_string(n_electrons);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    double n = 0;
    for (int i = 0; i < num_molecules; i++) {
      n += arma::sum(n_vl_molecule[i]);
    }

    return n + arma::sum(n_atom) + n_electrons;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute total number density given the vibrational levels of each molecule
  double kappa::Mixture::compute_n(const std::vector<arma::vec> &n_vl_molecule, double n_electrons) {

    #ifdef KAPPA_STRICT_CHECKS
    check_n_vl_molecule(n_vl_molecule);
    if (n_electrons < 0) {
      std::string error_string = "Number density of electrons is negative: n_electrons=" + std::to_string(n_electrons);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    double n = 0;
    for (int i = 0; i < num_molecules; i++) {
      n += arma::sum(n_vl_molecule[i]);
    }

    return n + n_electrons;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute total number density given molecules and atoms (electrons are optional)
  double kappa::Mixture::compute_n(const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons) {

    #ifdef KAPPA_STRICT_CHECKS
    check_n_molecule(n_molecule);
    check_n_atom(n_atom);
    if (n_electrons < 0) {
      std::string error_string = "Number density of electrons is negative: n_electrons=" + std::to_string(n_electrons);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    return arma::sum(n_molecule) + arma::sum(n_atom);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute total number density given the vector of partial number densities
  double kappa::Mixture::compute_n(const arma::vec &n) {

    #ifdef KAPPA_STRICT_CHECKS
    check_n(n);
    #endif

    return arma::sum(n);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the number density for each molecule given its vibrational levels
  arma::vec kappa::Mixture::compute_n_molecule(const std::vector<arma::vec> &n_vl_molecule) {

    #ifdef KAPPA_STRICT_CHECKS
    check_n_vl_molecule(n_vl_molecule);
    #endif

    arma::vec res(num_molecules);
    for (int i = 0; i < num_molecules; i++) {
      res[i] = arma::sum(n_vl_molecule[i]);
    }

    return res;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // *** Warning: pay attention to the order in which species are stored when coupling with CFD solvers!
  arma::vec kappa::Mixture::compute_density_array(const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons) {

    #ifdef KAPPA_STRICT_CHECKS
    check_n_vl_molecule(n_vl_molecule);
    check_n_atom(n_atom);
    if (n_electrons < 0) {
      std::string error_string = "Number density of electrons is negative: n_electrons=" + std::to_string(n_electrons);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    arma::vec res(num_molecules+num_atoms); // FIXME in case of electrons will have dimensional problems
    int i=0;
    for (i = 0; i < num_molecules; i++) {
      res[i] = molecules[i].mass * arma::sum(n_vl_molecule[i]);
    }
    for (i = 0; i < num_atoms; i++) {
      res[num_molecules + i] = atoms[i].mass * n_atom[i];
    }
    if (is_ionized) {
      res[num_molecules + num_atoms + i] = electron.mass * n_electrons;
    }

    return res;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  arma::vec kappa::Mixture::compute_density_array(const std::vector<arma::vec> &n_vl_molecule) {

    #ifdef KAPPA_STRICT_CHECKS
    check_n_vl_molecule(n_vl_molecule);
    if (n_electrons < 0) {
      std::string error_string = "Number density of electrons is negative: n_electrons=" + std::to_string(n_electrons);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    arma::vec res(num_molecules+num_atoms); // FIXME in case of elelctrons will have dimension problems
    int i=0;
    for (i = 0; i < num_molecules; i++) {
      res[i] = molecules[i].mass * arma::sum(n_vl_molecule[i]);
    }

    return res;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Mixture::compute_density(const arma::vec &n) {

    #ifdef KAPPA_STRICT_CHECKS
    check_n(n);
    #endif

    int i, j=0;
    double rho=0;
    for (i = 0; i < num_molecules; i++) {
      rho += molecules[i].mass * n[i];
      j++;
    }
    for (i = 0; i < num_atoms; i++) {
      rho += atoms[i].mass * n[j];
      j++;
    }
    if (is_ionized) {
      rho += electron.mass * n[n_particles-1];
    }
    return rho;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Mixture::compute_density(const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons) {

    #ifdef KAPPA_STRICT_CHECKS
    check_n_vl_molecule(n_vl_molecule);
    check_n_atom(n_atom);
    if (n_electrons < 0) {
      std::string error_string = "Number density of electrons is negative: n_electrons=" + std::to_string(n_electrons);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    int i=0;
    double rho=0;
    for (i = 0; i < num_molecules; i++) {
      rho += molecules[i].mass * arma::sum(n_vl_molecule[i]);
    }
    for (i = 0; i < num_atoms; i++) {
      rho += atoms[i].mass * n_atom[i];
    }
    if (is_ionized) {
      rho += electron.mass * n_electrons;
    }

    return rho;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Mixture::compute_density(const std::vector<arma::vec> &n_vl_molecule, double n_electrons) {

    #ifdef KAPPA_STRICT_CHECKS
    check_n_vl_molecule(n_vl_molecule);
    if (n_electrons < 0) {
      std::string error_string = "Number density of electrons is negative: n_electrons=" + std::to_string(n_electrons);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    double rho=0;
    for (int i=0; i<num_molecules; i++) {
      rho += molecules[i].mass * arma::sum(n_vl_molecule[i]);
    }
    if (is_ionized) {
      rho += electron.mass * n_electrons;
    }

    return rho;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Mixture::compute_density(const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons) {

    #ifdef KAPPA_STRICT_CHECKS
    check_n_molecule(n_molecule);
    check_n_atom(n_atom);
    if (n_electrons < 0) {
      std::string error_string = "Number density of electrons is negative: n_electrons=" + std::to_string(n_electrons);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    int i = 0;
    double rho = 0;
    for (i = 0; i < num_molecules; i++) {
      rho += molecules[i].mass * n_molecule[i];
    }
    for (i = 0; i < num_atoms; i++) {
      rho += atoms[i].mass * n_atom[i];
    }
    if (is_ionized) {
      rho += electron.mass * n_electrons;
    }

    return rho;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Mixture::compute_pressure(double T, const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons) {

    #ifdef KAPPA_STRICT_CHECKS
    check_n_vl_molecule(n_vl_molecule);
    check_n_atom(n_atom);
    if (n_electrons < 0) {
      std::string error_string = "Number density of electrons is negative: n_electrons=" + std::to_string(n_electrons);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    int i=0;
    double n=0;
    for (i = 0; i < num_molecules; i++) {
      n += arma::sum(n_vl_molecule[i]);
    }
    for (i = 0; i < num_atoms; i++) {
      n += n_atom[i];
    }
    if (is_ionized) {
      n += n_electrons;
    }

    // pressure
    return n * K_CONST_K * T;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute translational constant volume specific heat
  double kappa::Mixture::c_tr(const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons) {

    #ifdef KAPPA_STRICT_CHECKS
    check_n_vl_molecule(n_vl_molecule);
    check_n_atom(n_atom);
    if (n_electrons < 0) {
      std::string error_string = "Number density of electrons is negative: n_electrons=" + std::to_string(n_electrons);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    double n = 0;
    for (int i = 0; i < num_molecules; i++) {
      n += arma::sum(n_vl_molecule[i]);
    }

    // eq. 1.23
    return 1.5 * K_CONST_K * (n + arma::sum(n_atom) + n_electrons) / compute_density(n_vl_molecule, n_atom, n_electrons);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Mixture::c_tr(const std::vector<arma::vec> &n_vl_molecule, double n_electrons) {

    #ifdef KAPPA_STRICT_CHECKS
    check_n_vl_molecule(n_vl_molecule);
    if (n_electrons < 0) {
      std::string error_string = "Number density of electrons is negative: n_electrons=" + std::to_string(n_electrons);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    double n = 0;

    for (int i = 0; i < num_molecules; i++) {
      n += arma::sum(n_vl_molecule[i]);
    }

    return 1.5 * K_CONST_K * (n + n_electrons) / compute_density(n_vl_molecule, n_electrons);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Mixture::c_tr(const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons) {

    #ifdef KAPPA_STRICT_CHECKS
    check_n_molecule(n_molecule);
    check_n_atom(n_atom);
    if (n_electrons < 0) {
      std::string error_string = "Number density of electrons is negative: n_electrons=" + std::to_string(n_electrons);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    return 1.5 * K_CONST_K * (arma::sum(n_molecule) + arma::sum(n_atom) + n_electrons) / compute_density(n_molecule, n_atom, n_electrons);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Mixture::c_tr(const arma::vec &n) {

    #ifdef KAPPA_STRICT_CHECKS
    // if (T <= 0) {
    //   std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
    //   throw kappa::IncorrectValueException(error_string.c_str());
    // }
    check_n(n);
    #endif

    return 1.5 * K_CONST_K * arma::sum(n) / compute_density(n);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute rotational constant volume specific heat
  double kappa::Mixture::c_rot(	double T, const std::vector<arma::vec> &n_vl_molecule, const arma::vec &n_atom, double n_electrons) {

    #ifdef KAPPA_STRICT_CHECKS
    std::string error_string;
    if (T <= 0) {
      error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    check_n_vl_molecule(n_vl_molecule);
    check_n_atom(n_atom);
    if (n_electrons < 0) {
      error_string = "Number density of electrons is negative: n_electrons=" +
      std::to_string(n_electrons);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    // eq. 1.22
    double res = 0;
    for (int i = 0; i < num_molecules; i++) {
      for (int j = 0; j < molecules[i].num_vibr_levels[0]; j++) {
        res += c_rot(T, molecules[i], j, 0) * n_vl_molecule[i][j] * molecules[i].mass;
      }
    }
    return res / compute_density(n_vl_molecule, n_atom, n_electrons);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Mixture::c_rot(double T, const std::vector<arma::vec> &n_vl_molecule, double n_electrons) {

    #ifdef KAPPA_STRICT_CHECKS
    std::string error_string;
    if (T <= 0) {
      error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    check_n_vl_molecule(n_vl_molecule);
    if (n_electrons < 0) {
      error_string = "Number density of electrons is negative: n_electrons=" +
      std::to_string(n_electrons);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    double res = 0;
    // TODO: check if cache_on!
    for (int i = 0; i < num_molecules; i++) {
      for (int j = 0; j < molecules[i].num_vibr_levels[0]; j++) {
        res += c_rot(T, molecules[i], j, 0) * n_vl_molecule[i][j] * molecules[i].mass;
      }
    }
    return res / compute_density(n_vl_molecule, n_electrons);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Mixture::c_rot(	double T, const arma::vec &n_molecule, const arma::vec &n_atom, double n_electrons) {

    #ifdef KAPPA_STRICT_CHECKS
    std::string error_string;
    if (T <= 0) {
      error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    check_n_molecule(n_molecule);
    check_n_atom(n_atom);
    if (n_electrons < 0) {
      error_string = "Number density of electrons is negative: n_electrons=" +
      std::to_string(n_electrons);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (!all_rigid_rotators) {
      error_string = "Mixture contains non-rigid rotator molecules, cannot compute c_rot without vibrational level populations";
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    double res = 0;
    for (int i = 0; i < num_molecules; i++) {
      res += c_rot(T, molecules[i], 0, 0) * n_molecule[i] * molecules[i].mass;
    }

    return res / compute_density(n_molecule, n_atom, n_electrons);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Mixture::c_rot(	double T, const arma::vec &n) {

    #ifdef KAPPA_STRICT_CHECKS
    std::string error_string;
    if (T <= 0) {
      error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    check_n(n);
    if (!all_rigid_rotators) {
      error_string = "Mixture contains non-rigid rotator molecules, cannot compute c_rot without vibrational level populations";
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif
    if (num_molecules == 0) {
      return 0.0;
    }  else {
      double res = 0;
      for (int i = 0; i < num_molecules; i++) {
        res += c_rot(T, molecules[i], 0, 0) * n[i] * molecules[i].mass;
      }

      return res / compute_density(n);
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the left hand side of the algebraic system for bulk viscosity coefficients
  const arma::mat &kappa::Mixture::compute_bulk_viscosity_LHS(double T, kappa::models_omega model) {

    double n = this_total_n;
    double A_cd, eta_cd, coll_mass, n_ij, rm, kT=K_CONST_K * T, tmp; // xi_rot = rm * tau_rot
    int i, j, k, l, p1, p2, o1, o2;

    // molecule + molecule collisions for molecules of different species or at different vibrational levels
    for (i=0; i<num_molecules; i++) {
      // std::cout << " I = " << i << std::endl;
      for (j=i; j<num_molecules; j++) {
      // std::cout << " J = " << j << std::endl;

        coll_mass = interactions[inter_index(i, j)].collision_mass;

        A_cd   = 0.50  *      omega_22.at(i, j) / omega_11.at(i, j);	// eq. 5.62
        eta_cd = 0.625 * kT / omega_22.at(i, j); 			// 0.625 = 5/8
        rm     = 32    * n  * omega_22.at(i, j) / (5 * K_CONST_PI);

        for (k=0; k<molecules[i].num_vibr_levels[0]; k++) {

          p1 = 2 * (vl_offset[i] + k);
          o1 = vl_offset[i] + k;

          for (l=0; l<molecules[j].num_vibr_levels[0]; l++) {

            n_ij = this_n_vl_mol[i][k] * this_n_vl_mol[j][l];
            p2 = 2 * (vl_offset[j] + l);
            o2 = vl_offset[j] + l;

            if ((j!=i) || (l>k)) {

               bulk_viscosity_LHS.at(p1    , p2) = - 5 * kT * n_ij * coll_mass / (A_cd * eta_cd * (molecules[i].mass + molecules[j].mass))
                                                   + 4 *  T * n_ij * coll_mass * (molecules[i].mass * c_rot_arr[o1] / (rm * rot_rel_times.at(i, j))
                                                                               +  molecules[j].mass * c_rot_arr[o2] / (rm * rot_rel_times.at(j, i)))
                                                                               / (K_CONST_PI * eta_cd * (molecules[i].mass + molecules[j].mass)); // 1100
              //std::cout << "bulk_viscosity_LHS.at("<<p1<<","<<p2<<") = " << bulk_viscosity_LHS.at(p1, p2) / (n * n) << std::endl;

              bulk_viscosity_LHS.at(p1 + 1, p2) = - 4 * T * n_ij * molecules[j].mass * molecules[j].mass * c_rot_arr[o2] /
                                                  ((molecules[i].mass + molecules[j].mass) * K_CONST_PI * eta_cd * rm * rot_rel_times.at(j,i)); // 0110
              //std::cout << "bulk_viscosity_LHS.at("<<p1+1<<","<<p2<<") = " << bulk_viscosity_LHS.at(p1+1, p2) / (n * n) << std::endl;

              bulk_viscosity_LHS.at(p1    , p2 + 1) = - 4 * T * n_ij * molecules[i].mass * molecules[i].mass * c_rot_arr[o1] /
                                                      ((molecules[i].mass + molecules[j].mass) * K_CONST_PI * eta_cd * rm * rot_rel_times.at(i,j)); // 1001
              //std::cout << "bulk_viscosity_LHS.at("<<p1<<","<<p2+1<<") = " << bulk_viscosity_LHS.at(p1, p2+1) / (n * n) << std::endl;

              bulk_viscosity_LHS.at(p1 + 1, p2 + 1) = 0; // 0011
              //std::cout << "bulk_viscosity_LHS.at("<<p1+1<<","<<p2+1<<") = " << bulk_viscosity_LHS.at(p1+1, p2+1) / (n * n) << std::endl;

              // symmetrization
              bulk_viscosity_LHS.at(p2    , p1    ) = bulk_viscosity_LHS.at(p1    , p2    );
              bulk_viscosity_LHS.at(p2    , p1 + 1) = bulk_viscosity_LHS.at(p1 + 1, p2    );
              bulk_viscosity_LHS.at(p2 + 1, p1    ) = bulk_viscosity_LHS.at(p1    , p2 + 1);
              bulk_viscosity_LHS.at(p2 + 1, p1 + 1) = bulk_viscosity_LHS.at(p1 + 1, p2 + 1);
              //std::cout << "bulk_viscosity_LHS.at("<<p2<<","<<p1<<") = " << bulk_viscosity_LHS.at(p2, p1) / (n * n) << std::endl;
              //std::cout << "bulk_viscosity_LHS.at("<<p2<<","<<p1+1<<") = " << bulk_viscosity_LHS.at(p2, p1+1) / (n * n) << std::endl;
              //std::cout << "bulk_viscosity_LHS.at("<<p2+1<<","<<p1<<") = " << bulk_viscosity_LHS.at(p2+1, p1) / (n * n) << std::endl;
              //std::cout << "bulk_viscosity_LHS.at("<<p2+1<<","<<p1+1<<") = " << bulk_viscosity_LHS.at(p2+1, p1+1) / (n * n) << std::endl;

            }
          }
        }
      }
    }

    //std::cout << "bulk_viscosity_LHS MOL * MOL" << std::endl;

    // molecule + molecule collisions for identical molecular species at identical vibrational levels
    for (i=0; i<num_molecules; i++) {

       coll_mass = interactions[inter_index(i, i)].collision_mass;

       A_cd   = 0.5   *      omega_22.at(i, i) / omega_11.at(i, i);	// eqn. 5.62
       eta_cd = 0.625 * kT / omega_22.at(i, i); 			// 0.625 = 5/8
       rm     = 32    * n  * omega_22.at(i, i) / (5 * K_CONST_PI);

       for (k=0; k<molecules[i].num_vibr_levels[0]; k++) {

          p1 = 2 * (vl_offset[i] + k);
          o1 = vl_offset[i] + k;

          n_ij = this_n_vl_mol[i][k] * this_n_vl_mol[i][k];

          tmp = 4 * T * n_ij * molecules[i].mass * c_rot_arr[o1] / (K_CONST_PI * eta_cd * rm * rot_rel_times.at(i,i));

          bulk_viscosity_LHS.at(p1    , p1   )  =  tmp; // 1100
          bulk_viscosity_LHS.at(p1 + 1, p1   )  = -tmp; // 0110
          bulk_viscosity_LHS.at(p1    , p1 + 1) = -tmp; // 1001
          bulk_viscosity_LHS.at(p1 + 1, p1 + 1) =  tmp; // 0011
          //std::cout << "PRE MOLMOL" << std::endl;
          //std::cout << T << " " << n_ij << " " << c_rot_arr[o1] << " " << eta_cd << " " << rm << " " << rot_rel_times.at(i,i) << std::endl;
          //std::cout << "bulk_viscosity_LHS.at("<<p1<<","<<p1<<") = " << bulk_viscosity_LHS.at(p1, p1) / (n * n) << std::endl;
          //std::cout << "bulk_viscosity_LHS.at("<<p1+1<<","<<p1<<") = " << bulk_viscosity_LHS.at(p1+1, p1) / (n * n) << std::endl;
          //std::cout << "bulk_viscosity_LHS.at("<<p1<<","<<p1+1<<") = " << bulk_viscosity_LHS.at(p1, p1+1) / (n * n) << std::endl;
          //std::cout << "bulk_viscosity_LHS.at("<<p1+1<<","<<p1+1<<") = " << bulk_viscosity_LHS.at(p1+1, p1+1) / (n * n) << std::endl;
       }

       for (j=0; j<num_molecules; j++) {

         coll_mass = interactions[inter_index(i, j)].collision_mass;

         A_cd   = 0.5   *      omega_22.at(i, j) / omega_11.at(i, j);
         eta_cd = 0.625 * kT / omega_22.at(i, j);
         rm     = 32    * n  * omega_22.at(i, j) / (5 * K_CONST_PI);

         for (k = 0; k < molecules[i].num_vibr_levels[0]; k++) {

           p1 = 2 * (vl_offset[i] + k);
           o1 = vl_offset[i] + k;

           for (l=0; l<molecules[j].num_vibr_levels[0]; l++) {
             if ((k != l) || (i != j)) {

               n_ij = this_n_vl_mol[i][k] * this_n_vl_mol[j][l];
               o2 = vl_offset[j] + l;

               // eq. 5.78
               bulk_viscosity_LHS.at(p1    , p1    ) += 5 * kT * n_ij * coll_mass / ((molecules[i].mass + molecules[j].mass) * A_cd * eta_cd) +  4 * T  * n_ij * molecules[j].mass * molecules[j].mass * (molecules[i].mass * c_rot_arr[o1] / (rm * rot_rel_times.at(i, j)) +  molecules[j].mass * c_rot_arr[o2] / (rm * rot_rel_times.at(j, i))) / (K_CONST_PI * eta_cd * (molecules[i].mass + molecules[j].mass) * (molecules[i].mass + molecules[j].mass)); // 1100

               //std::cout << "bulk_viscosity_LHS.at("<<p1<<","<<p1<<") = " << bulk_viscosity_LHS.at(p1, p1) / (n * n) << std::endl;

               // eq. 5.80
               bulk_viscosity_LHS.at(p1 + 1, p1    ) -= 4 * T * n_ij * molecules[j].mass * molecules[i].mass * c_rot_arr[o1] / (rm * rot_rel_times.at(i, j) * K_CONST_PI * eta_cd * (molecules[i].mass + molecules[j].mass)); // 0110

               //std::cout << "bulk_viscosity_LHS.at("<<p1+1<<","<<p1<<") = " << bulk_viscosity_LHS.at(p1+1, p1) / (n * n) << std::endl;

               bulk_viscosity_LHS.at(p1    , p1 + 1) -= 4 * T * n_ij * molecules[j].mass * molecules[i].mass * c_rot_arr[o1] / (rm * rot_rel_times.at(i, j) * K_CONST_PI * eta_cd * (molecules[i].mass + molecules[j].mass)); // 1001

               //std::cout << "bulk_viscosity_LHS.at("<<p1<<","<<p1+1<<") = " << bulk_viscosity_LHS.at(p1, p1+1) / (n * n) << std::endl;

               // eq. 5.82
               bulk_viscosity_LHS.at(p1 + 1, p1 + 1) += 4 * T * n_ij * molecules[i].mass * c_rot_arr[o1] / (K_CONST_PI * eta_cd * rm * rot_rel_times.at(i, j)); // 0011

               //std::cout << "bulk_viscosity_LHS.at("<<p1+1<<","<<p1+1<<") = " << bulk_viscosity_LHS.at(p1+1, p1+1) / (n * n) << std::endl;
             }
           }
         }
       }

       for (j=0; j<num_atoms; j++) {

         coll_mass = interactions[inter_index(i, num_molecules + j)].collision_mass;

         A_cd   = 0.50  *      omega_22.at(i, num_molecules + j) / omega_11.at(i, num_molecules + j);
         eta_cd = 0.625 * kT / omega_22.at(i, num_molecules + j);
         rm     = 32    * n  * omega_22.at(i, num_molecules + j) / (5 * K_CONST_PI);

         for (k=0; k<molecules[i].num_vibr_levels[0]; k++) {

           p1 = 2 * (vl_offset[i] + k);
           o1 = vl_offset[i] + k;

           n_ij = this_n_vl_mol[i][k] * this_n_atom[j];

           bulk_viscosity_LHS.at(p1    , p1    ) += 5 * kT * n_ij * coll_mass / ((molecules[i].mass + atoms[j].mass) * A_cd * eta_cd)
                                                              +  4 * T * n_ij * atoms[j].mass * atoms[j].mass
                                                              * molecules[i].mass * c_rot_arr[o1] / (rm * rot_rel_times.at(i, num_molecules + j))
                                                              / (K_CONST_PI * eta_cd * (molecules[i].mass+atoms[j].mass) * (molecules[i].mass+atoms[j].mass)); // 1100


           bulk_viscosity_LHS.at(p1 + 1, p1    ) -= 4 * T * n_ij * atoms[j].mass * molecules[i].mass * c_rot_arr[o1] / (rm * rot_rel_times.at(i, num_molecules + j) * K_CONST_PI * eta_cd * (molecules[i].mass + atoms[j].mass)); // 0110

           bulk_viscosity_LHS.at(p1    , p1 + 1) -= 4 * T * n_ij * atoms[j].mass * molecules[i].mass * c_rot_arr[o1] / (rm * rot_rel_times.at(i, num_molecules + j) * K_CONST_PI * eta_cd * (molecules[i].mass + atoms[j].mass)); // 1001

           bulk_viscosity_LHS.at(p1 + 1, p1 + 1) += 4 * T * n_ij * molecules[i].mass * c_rot_arr[o1] / (K_CONST_PI * eta_cd * rm * rot_rel_times.at(i, num_molecules + j)); // 0011
        }
      }
    }

    // molecule + atom collisions
    for (i=0; i<num_molecules; i++) {
      for (j=0; j<num_atoms; j++) {

        coll_mass = interactions[inter_index(i, num_molecules + j)].collision_mass;

        A_cd   = 0.50  *      omega_22.at(i, num_molecules + j) / omega_11.at(i, num_molecules + j);
        eta_cd = 0.625 * kT / omega_22.at(i, num_molecules + j);
        rm     = 32    * n  * omega_22.at(i, num_molecules + j) / (5 * K_CONST_PI);

        for (k = 0; k < molecules[i].num_vibr_levels[0]; k++) {

          n_ij = this_n_vl_mol[i][k] * this_n_atom[j];

          p1 = 2 * (vl_offset[i] + k);
          o1 =  vl_offset[i] + k;

          bulk_viscosity_LHS.at(p1    , n_vibr_levels_total * 2 + j) = - 5 * kT * n_ij * coll_mass / (A_cd * eta_cd * (molecules[i].mass + atoms[j].mass))
                                                                       + 4 *  T * n_ij * coll_mass * (molecules[i].mass * c_rot_arr[o1] / (rm * rot_rel_times.at(i, num_molecules + j)))
                                                                        / (K_CONST_PI * eta_cd * (molecules[i].mass + atoms[j].mass)); // 1100
          //std::cout << "bulk_viscosity_LHS.at("<<p1<<","<<n_vibr_levels_total * 2 + j<<") = " << bulk_viscosity_LHS.at(p1, n_vibr_levels_total * 2 + j) / (n * n) << std::endl;

          bulk_viscosity_LHS.at(p1 + 1, n_vibr_levels_total * 2 + j) = -4 * T * n_ij * molecules[i].mass * molecules[i].mass * c_rot_arr[o1]
                                                                                    / ((molecules[i].mass + atoms[j].mass)
                                                                                    * K_CONST_PI * eta_cd * rm * rot_rel_times.at(i, num_molecules + j)); // 0110
          //std::cout << "bulk_viscosity_LHS.at("<<p1+1<<","<<n_vibr_levels_total * 2 + j<<") = " << bulk_viscosity_LHS.at(p1+1, n_vibr_levels_total * 2 + j) / (n * n) << std::endl;

          // symmetrization
          bulk_viscosity_LHS.at(n_vibr_levels_total * 2 + j, p1    ) = bulk_viscosity_LHS.at(p1    , n_vibr_levels_total * 2 + j);
          bulk_viscosity_LHS.at(n_vibr_levels_total * 2 + j, p1 + 1) = bulk_viscosity_LHS.at(p1 + 1, n_vibr_levels_total * 2 + j);
          //std::cout << "bulk_viscosity_LHS.at("<<n_vibr_levels_total * 2 + j<<","<<p1<<") = " << bulk_viscosity_LHS.at(n_vibr_levels_total * 2 + j,p1) / (n * n) << std::endl;
          //std::cout << "bulk_viscosity_LHS.at("<<n_vibr_levels_total * 2 + j<<","<<p1+1<<") = " << bulk_viscosity_LHS.at(n_vibr_levels_total * 2 + j,p1+1) / (n * n) << std::endl;
        }
      }
    }

    // atoms + atom collisions for different atomic species
    for (i=0; i<num_atoms-1; i++) {
      for (j=i+1; j<num_atoms; j++) {

        coll_mass = interactions[inter_index(num_molecules + i, num_molecules + j)].collision_mass;

        n_ij = this_n_atom[i] * this_n_atom[j];

        A_cd = 0.5 * omega_22.at(num_molecules + i, num_molecules + j) / omega_11.at(num_molecules + i, num_molecules + j);
        eta_cd = 0.625 * kT / omega_22.at(num_molecules + i, num_molecules + j);

        bulk_viscosity_LHS.at(n_vibr_levels_total * 2 + i, n_vibr_levels_total * 2 + j) = -5 * kT * n_ij * coll_mass / (A_cd * eta_cd * (atoms[i].mass + atoms[j].mass)) ; // 1100
        //std::cout << "bulk_viscosity_LHS.at("<<n_vibr_levels_total * 2 + i<<","<<n_vibr_levels_total * 2 + j<<") = " << bulk_viscosity_LHS.at(n_vibr_levels_total * 2 + i, n_vibr_levels_total * 2 + j) / (n * n) << std::endl;

        // symmetrization
        bulk_viscosity_LHS.at(n_vibr_levels_total * 2 + j, n_vibr_levels_total * 2 + i) = bulk_viscosity_LHS.at(n_vibr_levels_total * 2 + i, n_vibr_levels_total * 2 + j);
        //std::cout << "bulk_viscosity_LHS.at("<<n_vibr_levels_total * 2 + j<<","<<n_vibr_levels_total * 2 + i<<") = " << bulk_viscosity_LHS.at(n_vibr_levels_total * 2 + j, n_vibr_levels_total * 2 + i) / (n * n) << std::endl;

      }
    }

    // atom + atom collisions for identical atomic species
    for (i=0; i<num_atoms; i++) {

      bulk_viscosity_LHS.at(n_vibr_levels_total * 2 + i, n_vibr_levels_total * 2 + i) = 0; // 1100

      for (j=0; j<num_molecules; j++) {

        coll_mass = interactions[inter_index(num_molecules + i, j)].collision_mass;

        A_cd   = 0.50  *      omega_22.at(num_molecules + i, j) / omega_11.at(num_molecules + i, j);
        eta_cd = 0.625 * kT / omega_22.at(num_molecules + i, j);
        rm     = 32    * n  * omega_22.at(num_molecules + i, j) / (5 * K_CONST_PI);

        for (k=0; k<molecules[j].num_vibr_levels[0]; k++) {

          n_ij = this_n_atom[i] * this_n_vl_mol[j][k];
          o1 =  vl_offset[i] + k;

          bulk_viscosity_LHS.at(n_vibr_levels_total * 2 + i, n_vibr_levels_total * 2 + i) += 5 * kT * n_ij * coll_mass /
                          ((atoms[i].mass + molecules[j].mass) * A_cd * eta_cd)
 											    +  4 *  T * n_ij * molecules[j].mass * molecules[j].mass
                          * molecules[j].mass * c_rot_arr[o1] / (rm * rot_rel_times.at(j, num_molecules + i))
                          / (K_CONST_PI*eta_cd*(atoms[i].mass+molecules[j].mass)*(atoms[i].mass+molecules[j].mass)); // 1100;
        }
      }

      for (j=0; j<num_atoms; j++) {
        if (j!=i) {

          coll_mass = interactions[inter_index(num_molecules + i, num_molecules + j)].collision_mass;

          n_ij = this_n_atom[i] * this_n_atom[j];

          A_cd   = 0.5   *      omega_22.at(num_molecules + i, num_molecules + j) / omega_11.at(num_molecules + i, num_molecules + j);
          eta_cd = 0.625 * kT / omega_22.at(num_molecules + i, num_molecules + j);

          bulk_viscosity_LHS.at(n_vibr_levels_total * 2 + i, n_vibr_levels_total * 2 + i) += 5 * kT * n_ij * coll_mass / ((atoms[i].mass + atoms[j].mass) * A_cd * eta_cd); // 1100;
        }
      }
    }

    bulk_viscosity_LHS /= (n * n);

    // std::cout << "NONRIGID, BEFORE: " << arma::det(bulk_viscosity_LHS * 1e15) << std::endl;

    // condition so that the system is not singular
    // we divide by 10^20 so that the order of magnitude of the terms is roughly the same as the brackets
    for (i=0; i<num_molecules; i++) {
      for (k=0; k<molecules[i].num_vibr_levels[0]; k++) {

        p1 = 2 * (vl_offset[i] + k);

        bulk_viscosity_LHS.at(0, p1) = this_n_vl_mol[i][k] * this_ctr / (n * 1e20);
        bulk_viscosity_LHS.at(0, p1 + 1) = this_n_vl_mol[i][k] * c_rot_arr[vl_offset[i] + k] / (n * 1e20);

      }
    }

    for (i=0; i<num_atoms; i++) {
      bulk_viscosity_LHS.at(0,  n_vibr_levels_total * 2 + i) = this_n_atom[i] * this_ctr / (n * 1e20);
    }

    // std::cout << "NONRIGID, AFTER: " << arma::det(bulk_viscosity_LHS * 1e15) << std::endl;

//    /////////
//    for (i=0; i<num_molecules; i++) {
//      for (j=i; j<num_molecules; j++) {
//        for (k=0; k<molecules[i].num_vibr_levels[0]; k++) {
//          p1 = 2 * (vl_offset[i] + k);
//          o1 = vl_offset[i] + k;
//          for (l=0; l<molecules[j].num_vibr_levels[0]; l++) {
//            p2 = 2 * (vl_offset[j]);
//            o2 = vl_offset[j] + l;
//            std::cout << "bulk_viscosity_LHS.at("<<p1<<","<<p2<<") = " << bulk_viscosity_LHS.at(p1, p2) << std::endl;
//            std::cout << "bulk_viscosity_LHS.at("<<p1+1<<","<<p2<<") = " << bulk_viscosity_LHS.at(p1+1, p2) << std::endl;
//            std::cout << "bulk_viscosity_LHS.at("<<p1<<","<<p2+1<<") = " << bulk_viscosity_LHS.at(p1, p2+1) << std::endl;
//            std::cout << "bulk_viscosity_LHS.at("<<p1+1<<","<<p2+1<<") = " << bulk_viscosity_LHS.at(p1+1, p2+1) << std::endl;
//          }
//        }
//      }
//    }

    /////////

    //std::cout << "T = " << T << std::endl;
    //std::cout << "bulk_viscosity_LHS" << bulk_viscosity_LHS << std::endl;
    return bulk_viscosity_LHS;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the right hand side of the algebraic system for bulk viscosity coefficients
  const arma::vec &kappa::Mixture::compute_bulk_viscosity_RHS(double T) {

    int i, j=0, k;

    //std::cout << "this_ctr" << this_ctr << std::endl;
    //std::cout << "this_crot" << this_crot << std::endl;
    //std::cout << "this_total_n" << this_total_n << std::endl;
    //std::cout << "c_rot_arr" << c_rot_arr << std::endl;

    //std::cout << "PRINT" << std::endl;
    for (i=0; i<num_molecules; i++) {
      for (k=0; k<molecules[i].num_vibr_levels[0]; k++) {

        //std::cout << molecules[i].mass << std::endl;
        //std::cout << i << " " << k << " " << " " << this_n_vl_mol[i][k] << " "
        //                                         << this_n_vl_mol[i][k] / this_total_n << " "
        //                                         << this_n_vl_mol[i][k] * molecules[i].mass / this_total_dens << std::endl;

        bulk_viscosity_RHS[2 * j] = -this_n_vl_mol[i][k] * this_crot / this_total_n; // even
        bulk_viscosity_RHS[2 * j + 1] = this_n_vl_mol[i][k] * molecules[i].mass * c_rot_arr[j] / this_total_dens; // odd
        j++;
      }
    }

    for (i=0; i<num_atoms; i++) {

      //std::cout << atoms[i].mass << std::endl;

      //bulk_viscosity_rigid_rot_RHS[num_molecules * 2 + i] = -this_n_atom[i] * this_crot / this_total_n;
      //bulk_viscosity_RHS[num_molecules * 2 + i] = -this_n_atom[i] * this_crot / this_total_n;
      //lk
      //std::cout << i << " " << -this_n_atom[i] / this_total_n << std::endl;
      bulk_viscosity_RHS[n_vibr_levels_total * 2 + i] = -this_n_atom[i] * this_crot / this_total_n;
    }

    bulk_viscosity_RHS[0] = 0; // normalization condition as to make the system linearly independent
    bulk_viscosity_RHS /= (this_ctr + this_crot);

    //std::cout << "T = " << T << std::endl;
    //std::cout << "bulk_viscosity_RHS" << bulk_viscosity_RHS << std::endl;
    //std::cout << "size bulk_viscosity_RHS" << bulk_viscosity_RHS.size() << std::endl;
    return bulk_viscosity_RHS;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the bulk viscosity coefficients
  const arma::vec &kappa::Mixture::compute_bulk_viscosity_coeffs(double T, kappa::models_omega model) {

    compute_bulk_viscosity_LHS(T, model);
    compute_bulk_viscosity_RHS(T);

    std::cout << "compute_bulk_viscosity_coeffs" << std::endl;
    std::cout << "bulk_viscosity_RHS" << std::endl;
    std::cout << bulk_viscosity_RHS << std::endl;
    std::cout << "bulk_viscosity_LHS" << std::endl;
    std::cout << bulk_viscosity_LHS << std::endl;
    bulk_viscosity_coeffs = arma::solve(bulk_viscosity_LHS * 1e20, bulk_viscosity_RHS) * 1e20;

    return bulk_viscosity_coeffs;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the left hand side of the algebraic system for bulk viscosity coefficients in the rigid rotator approximation
  const arma::mat &kappa::Mixture::compute_bulk_viscosity_rigid_rot_LHS(double T, kappa::models_omega model) {

    double n = this_total_n;

    double A_cd, eta_cd, coll_mass, n_ij, rm, kT=K_CONST_K * T, tmp; // xi_rot = rm * tau_rot
    int i, j;

    // molecule + molecule collisions for different molecular species
    for (i=0; i<num_molecules - 1; i++) {
      for (j=i+1; j<num_molecules; j++) {

        coll_mass = interactions[inter_index(i, j)].collision_mass;

        n_ij = this_n_molecules[i] * this_n_molecules[j];

        A_cd = 0.5 * omega_22.at(i, j) / omega_11.at(i, j);
        eta_cd = 0.625 * kT / omega_22.at(i, j);
        rm = 32 * n * omega_22.at(i, j) / (5 * K_CONST_PI);

        bulk_viscosity_rigid_rot_LHS.at(2*i,2*j) = -5 * kT * n_ij * coll_mass /
                                                   (A_cd * eta_cd * (molecules[i].mass + molecules[j].mass))
                                                    + 4 * T * n_ij * coll_mass *
                                                   (molecules[i].mass * c_rot_rigid_rot_arr[i] / (rm * rot_rel_times.at(i, j))
                                                  + molecules[j].mass * c_rot_rigid_rot_arr[j] / (rm * rot_rel_times.at(j, i)))
                                                  / (K_CONST_PI * eta_cd * (molecules[i].mass + molecules[j].mass)); // 1100

        bulk_viscosity_rigid_rot_LHS.at(2*i+1,2*j) = -4*T * n_ij * molecules[j].mass * molecules[j].mass * c_rot_rigid_rot_arr[j]
                                                      / ((molecules[i].mass + molecules[j].mass) * K_CONST_PI *
                                                          eta_cd * rm * rot_rel_times.at(j,i)); // 1001

        bulk_viscosity_rigid_rot_LHS.at(2*i,2*j+1) = -4*T * n_ij * molecules[i].mass * molecules[i].mass * c_rot_rigid_rot_arr[i]
                                                               / ((molecules[i].mass + molecules[j].mass) *
                                                                   K_CONST_PI * eta_cd * rm * rot_rel_times.at(i,j)); // 0110

        bulk_viscosity_rigid_rot_LHS.at(2*i+1,2*j+1) = 0; // 0011

        bulk_viscosity_rigid_rot_LHS.at(2 * j, 2 * i) = bulk_viscosity_rigid_rot_LHS.at(2 * i, 2 * j);
        bulk_viscosity_rigid_rot_LHS.at(2 * j, 2 * i + 1) = bulk_viscosity_rigid_rot_LHS.at(2 * i + 1, 2 * j);

        bulk_viscosity_rigid_rot_LHS.at(2 * j + 1, 2 * i) = bulk_viscosity_rigid_rot_LHS.at(2 * i, 2 * j + 1);
        bulk_viscosity_rigid_rot_LHS.at(2 * j + 1, 2 * i + 1) = bulk_viscosity_rigid_rot_LHS.at(2 * i + 1, 2 * j + 1);
      }
    }

    // molecule + molecule collisions for identical molecular species
    for (i=0; i<num_molecules; i++) {

      coll_mass = interactions[inter_index(i, i)].collision_mass;
      A_cd = 0.5 * omega_22.at(i, i) / omega_11.at(i, i);
      eta_cd = 0.625 * kT / omega_22.at(i, i);
      rm = 32 * n * omega_22.at(i, i) / (5 * K_CONST_PI);
      n_ij = this_n_molecules[i] * this_n_molecules[i];

      tmp = 4 * T * n_ij * molecules[i].mass * c_rot_rigid_rot_arr[i] / (K_CONST_PI * eta_cd * rm * rot_rel_times.at(i,i));

      bulk_viscosity_rigid_rot_LHS.at(2 * i, 2 * i) = tmp; 		// 1100
      bulk_viscosity_rigid_rot_LHS.at(2 * i + 1, 2 * i) = -tmp; 	// 1001

      bulk_viscosity_rigid_rot_LHS.at(2 * i, 2 * i + 1) = -tmp; 	// 0110
      bulk_viscosity_rigid_rot_LHS.at(2 * i + 1, 2 * i + 1) = tmp; 	// 0011

      for (j=0; j<num_molecules; j++) {
        if (j != i) {

          coll_mass = interactions[inter_index(i, j)].collision_mass;
          n_ij = this_n_molecules[i] * this_n_molecules[j];
          A_cd = 0.5 * omega_22.at(i, j) / omega_11.at(i, j);
          rm = 32 * n * omega_22.at(i, j) / (5 * K_CONST_PI);
          eta_cd = 0.625 * kT / omega_22.at(i, j);

           bulk_viscosity_rigid_rot_LHS.at(2*i,2*i) += 5*kT * n_ij * coll_mass /
                                                       ((molecules[i].mass + molecules[j].mass) * A_cd * eta_cd) +
                                                       4 * T * n_ij * molecules[j].mass * molecules[j].mass
                                                       * (molecules[i].mass * c_rot_rigid_rot_arr[i] /
                                                         (rm * rot_rel_times.at(i, j)) + molecules[j].mass *
                                                         c_rot_rigid_rot_arr[j] / (rm * rot_rel_times.at(j, i)))
                                                         / (K_CONST_PI * eta_cd * (molecules[i].mass + molecules[j].mass) *
                                                         (molecules[i].mass + molecules[j].mass)); // 1100

           bulk_viscosity_rigid_rot_LHS.at(2*i+1,2*i) -= 4*T * n_ij * molecules[j].mass * molecules[i].mass *
                                                         c_rot_rigid_rot_arr[i] / (rm * rot_rel_times.at(i, j) *
                                                         K_CONST_PI * eta_cd * (molecules[i].mass + molecules[j].mass)); // 1001

           bulk_viscosity_rigid_rot_LHS.at(2*i,2*i+1) -= 4*T * n_ij * molecules[j].mass * molecules[i].mass
                                                         * c_rot_rigid_rot_arr[i] / (rm * rot_rel_times.at(i, j)
                                                         * K_CONST_PI * eta_cd * (molecules[i].mass + molecules[j].mass)); // 0110

            bulk_viscosity_rigid_rot_LHS.at(2*i+1,2*i+1) += 4*T * n_ij * molecules[i].mass * c_rot_rigid_rot_arr[i] /
                                                            (K_CONST_PI * eta_cd * rm * rot_rel_times.at(i, j)); // 0011

        }
      }

      for (j=0; j<num_atoms; j++) {

        coll_mass = interactions[inter_index(i, num_molecules + j)].collision_mass;
        n_ij = this_n_molecules[i] * this_n_atom[j];
        A_cd = 0.5 * omega_22.at(i, num_molecules + j) / omega_11.at(i, num_molecules + j);
        eta_cd = 0.625 * kT / omega_22.at(i, num_molecules + j);
        rm = 32 * n * omega_22.at(i, num_molecules + j) / (5 * K_CONST_PI);

        bulk_viscosity_rigid_rot_LHS.at(2 * i, 2 * i) += 5 * kT * n_ij * coll_mass / ((molecules[i].mass + atoms[j].mass) * A_cd * eta_cd) +  4 * T * n_ij * atoms[j].mass * atoms[j].mass
                                                               * molecules[i].mass * c_rot_rigid_rot_arr[i] / (rm * rot_rel_times.at(i, num_molecules + j))
                                                               / (K_CONST_PI * eta_cd * (molecules[i].mass + atoms[j].mass) * (molecules[i].mass + atoms[j].mass)); // 1100

        bulk_viscosity_rigid_rot_LHS.at(2 * i + 1, 2 * i) -= 4 * T * n_ij * atoms[j].mass * molecules[i].mass * c_rot_rigid_rot_arr[i]
                                                                   / (rm * rot_rel_times.at(i, num_molecules + j) * K_CONST_PI * eta_cd * (molecules[i].mass + atoms[j].mass)); // 1001

        bulk_viscosity_rigid_rot_LHS.at(2 * i, 2 * i + 1) -= 4 * T * n_ij * atoms[j].mass * molecules[i].mass * c_rot_rigid_rot_arr[i]
                                                                   / (rm * rot_rel_times.at(i, num_molecules + j) * K_CONST_PI * eta_cd * (molecules[i].mass + atoms[j].mass)); // 0110
        bulk_viscosity_rigid_rot_LHS.at(2 * i + 1, 2 * i + 1) += 4 * T * n_ij * molecules[i].mass * c_rot_rigid_rot_arr[i] / (K_CONST_PI * eta_cd * rm * rot_rel_times.at(i, num_molecules + j)); // 0011
      }
    }

    // molecule + atom collisions
    for (i=0; i<num_molecules; i++) {
      for (j=0; j<num_atoms; j++) {

        coll_mass = interactions[inter_index(i, num_molecules + j)].collision_mass;
        n_ij = this_n_molecules[i] * this_n_atom[j];
        A_cd = 0.5 * omega_22.at(i, num_molecules + j) / omega_11.at(i, num_molecules + j);
        eta_cd = 0.625 * kT / omega_22.at(i, num_molecules + j);
        rm = 32 * n * omega_22.at(i, num_molecules + j) / (5 * K_CONST_PI);

        bulk_viscosity_rigid_rot_LHS.at(2 * i, num_molecules * 2 + j) = -5 * kT * n_ij * coll_mass / (A_cd * eta_cd * (molecules[i].mass + atoms[j].mass))
                                                                          + 4 * T * n_ij * coll_mass * (molecules[i].mass * c_rot_rigid_rot_arr[i] / (rm * rot_rel_times.at(i, num_molecules + j)))
                                                                              / (K_CONST_PI * eta_cd * (molecules[i].mass + atoms[j].mass));  // 1100
        bulk_viscosity_rigid_rot_LHS.at(2 * i + 1, num_molecules * 2 + j) = -4 * T * n_ij * molecules[i].mass * molecules[i].mass * c_rot_rigid_rot_arr[i]
                                                                                 / ((molecules[i].mass + atoms[j].mass) * K_CONST_PI * eta_cd * rm * rot_rel_times.at(i, num_molecules + j)); // 0110

        // symmetrization
        bulk_viscosity_rigid_rot_LHS.at(num_molecules * 2 + j, 2 * i) = bulk_viscosity_rigid_rot_LHS.at(2 * i, num_molecules * 2 + j);
        bulk_viscosity_rigid_rot_LHS.at(num_molecules * 2 + j, 2 * i + 1) = bulk_viscosity_rigid_rot_LHS.at(2 * i + 1, num_molecules * 2 + j);
      }
    }

    for (i=0; i<num_atoms-1; i++) { // atoms + atom collisions for different atomic species
      for (j=i+1; j<num_atoms; j++) {

        coll_mass = interactions[inter_index(num_molecules + i, num_molecules + j)].collision_mass;
        n_ij = this_n_atom[i] * this_n_atom[j];
        A_cd = 0.5 * omega_22.at(num_molecules + i, num_molecules + j) / omega_11.at(num_molecules + i, num_molecules + j);
        eta_cd = 0.625 * kT / omega_22.at(num_molecules + i, num_molecules + j); // 0.625 = 5/8

        bulk_viscosity_rigid_rot_LHS.at(num_molecules * 2 + i, num_molecules * 2 + j) = -5 * kT * n_ij * coll_mass / (A_cd * eta_cd * (atoms[i].mass + atoms[j].mass)) ; // 1100

        // symmetrization
        bulk_viscosity_rigid_rot_LHS.at(num_molecules * 2 + j, num_molecules * 2 + i) = bulk_viscosity_rigid_rot_LHS.at(num_molecules * 2 + i, num_molecules * 2 + j);
      }
    }

    for (i=0; i<num_atoms; i++) { // atom + atom collisions for identical atomic species

      bulk_viscosity_rigid_rot_LHS.at(num_molecules * 2 + i, num_molecules * 2 + i) = 0; // 0000

      for (j=0; j<num_molecules; j++) {

        coll_mass = interactions[inter_index(num_molecules + i, j)].collision_mass;
        n_ij = this_n_atom[i] * this_n_molecules[j];
        A_cd = 0.5 * omega_22.at(num_molecules + i, j) / omega_11.at(num_molecules + i, j);
        rm = 32 * n * omega_22.at(num_molecules + i, j) / (5 * K_CONST_PI);
        eta_cd = 0.625 * kT / omega_22.at(num_molecules + i, j); // 0.625 = 5/8

        bulk_viscosity_rigid_rot_LHS.at(num_molecules * 2 + i, num_molecules * 2 + i) += 5 * kT * n_ij * coll_mass / ((atoms[i].mass + molecules[j].mass) * A_cd * eta_cd) +  4 * T * n_ij * molecules[j].mass * molecules[j].mass
                                                                                           * molecules[j].mass * c_rot_rigid_rot_arr[j] / (rm * rot_rel_times.at(j, num_molecules + i))
                                                                                           / (K_CONST_PI * eta_cd * (atoms[i].mass + molecules[j].mass) * (atoms[i].mass + molecules[j].mass)); // 1100;
      }

      for (j=0; j<num_atoms; j++) {
        if (j!=i) {

          coll_mass = interactions[inter_index(num_molecules + i, num_molecules + j)].collision_mass;
          n_ij = this_n_atom[i] * this_n_atom[j];
          A_cd = 0.5 * omega_22.at(num_molecules + i, num_molecules + j) / omega_11.at(num_molecules + i, num_molecules + j);
          eta_cd = 0.625 * kT / omega_22.at(num_molecules + i, num_molecules + j); // 0.625 = 5/8

          bulk_viscosity_rigid_rot_LHS.at(num_molecules * 2 + i, num_molecules * 2 + i) += 5 * kT * n_ij * coll_mass / ((atoms[i].mass + atoms[j].mass) * A_cd * eta_cd); // 1100;
        }
      }
    }

    bulk_viscosity_rigid_rot_LHS /= (n * n);

    // double ctr = c_tr(n_vl_molecule, n_atom);
    for (i=0; i<num_molecules; i++) {

      bulk_viscosity_rigid_rot_LHS.at(0, 2 * i) = this_n_molecules[i] * this_ctr / (n * 1e20);
      bulk_viscosity_rigid_rot_LHS.at(0, 2 * i + 1) = this_n_molecules[i] * c_rot_rigid_rot_arr[i] / (n * 1e20);
    }
    for (i=0; i<num_atoms; i++) {
      bulk_viscosity_rigid_rot_LHS.at(0, num_molecules * 2 + i) = this_n_atom[i] * this_ctr / (n * 1e20);
    }

    return bulk_viscosity_rigid_rot_LHS;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the right hand side of the algebraic system for bulk viscosity coefficients in the rigid rotator approximation
  const arma::vec &kappa::Mixture::compute_bulk_viscosity_rigid_rot_RHS(double T) {

    int i;

    for (i = 0; i < num_molecules; i++) {
      bulk_viscosity_rigid_rot_RHS[2 * i    ] = -this_n_molecules[i] * this_crot / this_total_n;
      bulk_viscosity_rigid_rot_RHS[2 * i + 1] =  this_n_molecules[i] * molecules[i].mass * c_rot_rigid_rot_arr[i] / this_total_dens;
    }

    for (i = 0; i < num_atoms; i++) {
      bulk_viscosity_rigid_rot_RHS[num_molecules * 2 + i] = -this_n_atom[i] * this_crot / this_total_n;
    }

    bulk_viscosity_rigid_rot_RHS[0] = 0; // normalization condition
    bulk_viscosity_rigid_rot_RHS /= (this_ctr + this_crot);

    //std::cout << this_ctr << " " << this_crot << " " <<  " " << c_rot_rigid_rot_arr[0] << std::endl;
    //std::cout << bulk_viscosity_rigid_rot_RHS << std::endl;
    return bulk_viscosity_rigid_rot_RHS;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the bulk viscosity coefficients in the rigid rotator approximation
  const arma::vec &kappa::Mixture::compute_bulk_viscosity_rigid_rot_coeffs(double T, kappa::models_omega model) {

    compute_bulk_viscosity_rigid_rot_LHS(T, model);
    compute_bulk_viscosity_rigid_rot_RHS(T);

    std::cout << "compute_bulk_viscosity_coeffs_RR" << std::endl;
    std::cout << "bulk_viscosity_RHS_RR" << std::endl;
    std::cout << bulk_viscosity_rigid_rot_RHS << std::endl;
    std::cout << "bulk_viscosity_LHS_RR" << std::endl;
    std::cout << bulk_viscosity_rigid_rot_LHS << std::endl;

    bulk_viscosity_rigid_rot_coeffs = arma::solve(bulk_viscosity_rigid_rot_LHS * 1e17, bulk_viscosity_rigid_rot_RHS) * 1e17;

    return bulk_viscosity_rigid_rot_coeffs;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // bulk viscosity, eq. 5.11
  double kappa::Mixture::bulk_viscosity(double T, kappa::models_omega model) {

    int i;
    double res=0;

    if (all_rigid_rotators) {

      compute_bulk_viscosity_rigid_rot_coeffs(T, model);

      for (i=0; i<num_molecules; i++) {
        res += this_n_molecules[i] * bulk_viscosity_rigid_rot_coeffs[2 * i];
      }

      for (i=0; i<num_atoms; i++) {
        res += this_n_atom[i] * bulk_viscosity_rigid_rot_coeffs[num_molecules * 2 + i];
      }
    } else {
      int k, j=0;
      compute_bulk_viscosity_coeffs(T, model);
      for (i=0; i<num_molecules; i++) {
        for (k=0; k<molecules[i].num_vibr_levels[0]; k++) {
          res += this_n_vl_mol[i][k] * bulk_viscosity_coeffs[2 * j];
          j++;
        }
      }

      for (i=0; i<num_atoms; i++) {
        res += this_n_atom[i] * bulk_viscosity_coeffs[2 * n_vibr_levels_total + i];
      }
    }

    std::cout << "T = " << T << std::endl;
    std::cout << "bulk_viscosity = " << -res * K_CONST_K * T / this_total_n << std::endl;
    return -res * K_CONST_K * T / this_total_n;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the left hand side of the algebraic system for thermal conductivity coefficients
  const arma::mat &kappa::Mixture::compute_thermal_conductivity_LHS(double T, kappa::models_omega model) {

    double n = this_total_n;
    double A_cd, B_cd, C_cd, D_cd, coll_mass, D_cd_rot, n_ij, rm, kT=K_CONST_K * T; // xi_rot = rm * tau_rot
    int i, j, k, l, p1, p2, o1, o2;

    // molecule + molecule collisions for molecules of different species or at different vibrational levels
    for (i=0; i<num_molecules; i++) {
      for (j=i; j<num_molecules; j++) {

        coll_mass = interactions[inter_index(i, j)].collision_mass;

        A_cd = ( 0.5 ) *      omega_22.at(i, j) / omega_11.at(i, j); // eq. 5.60
        B_cd = (1./3.) * (5 * omega_12.at(i, j) - omega_13.at(i, j)) / omega_11.at(i, j);
        C_cd = (1./3.) *      omega_12.at(i, j) / omega_11.at(i, j);
        D_cd = (3./16) *      kT/(n * coll_mass * omega_11.at(i, j)); // eq. 5.61
        rm   = ( 32  ) *  n * omega_22.at(i, j) / (5 * K_CONST_PI);

        std::cout << rm << std::endl;
        std::cout << "coll_mass = " << coll_mass << std::endl;
        std::cout << "Omega11 = " << omega_11.at(i, j) << std::endl;
        std::cout << "Omega12 = " << omega_12.at(i, j) << std::endl;
        std::cout << "Omega13 = " << omega_13.at(i, j) << std::endl;
        std::cout << "Omega22 = " << omega_22.at(i, j) << std::endl;
        std::cout << "Acd = " << A_cd << std::endl;
        std::cout << "Bcd = " << B_cd << std::endl;
        std::cout << "Ccd = " << C_cd << std::endl;
        std::cout << "Dcd = " << D_cd << std::endl;
        std::cout << "n = " << n << std::endl;

        // loop over vibrational levels of each molecule i
        for (k=0; k<molecules[i].num_vibr_levels[0]; k++) {

          p1 = 3 * (vl_offset[i] + k); // cols
          o1 =      vl_offset[i] + k ; // rows

          // loop over vibrational levels of each molecule j
          for (l=0; l<molecules[j].num_vibr_levels[0]; l++) {

            // std::cout << this_n_vl_mol[i][k] << std::endl;
            // std::cout << this_n_vl_mol[j][l] << std::endl;
            n_ij = this_n_vl_mol[i][k] * this_n_vl_mol[j][l];

            //std::cout << this_n_vl_mol[i][k] << " " << this_n_vl_mol[j][l] << std::endl;
            //std::cout << " bo, where am i" << std::endl;

            p2 = 3 * (vl_offset[j] + l);
            o2 =      vl_offset[j] + l ;

            if ((j!=i) || (l>k)) { // d,k ∈ S_ci, i.e. δ_ik * δ_cd = 0

              //std::cout << "rot_rel_times.at(i, j) << " " << rot_rel_times.at(j, i)" << std::endl;
              //std::cout << rot_rel_times.at(i, j) << " " << rot_rel_times.at(j, i) << std::endl;

              thermal_conductivity_LHS.at(p1    , p2    ) = -1.5 * kT * n_ij / (n * D_cd); // 0000
              thermal_conductivity_LHS.at(p1 + 1, p2    ) = 0.75 * kT * n_ij * (6 * C_cd - 5) * molecules[j].mass / ((molecules[i].mass + molecules[j].mass) * n * D_cd); // 1000
              thermal_conductivity_LHS.at(p1 + 2, p2    ) = 0; // 0010
              thermal_conductivity_LHS.at(p1    , p2 + 1) = 0.75 * kT * n_ij * (6 * C_cd - 5) * molecules[i].mass / ((molecules[i].mass + molecules[j].mass) * n * D_cd); // 0100 = 1000
              thermal_conductivity_LHS.at(p1 + 1, p2 + 1) = -1.5 * kT * n_ij * coll_mass / ((molecules[i].mass + molecules[j].mass) * n * D_cd) * (13.75 - 3 * B_cd - 4 * A_cd - (20./3.) * A_cd / (K_CONST_K * K_CONST_PI)
                                                                             * (molecules[i].mass * c_rot_arr[o1] / (rm * rot_rel_times.at(i, j))
                                                                             +  molecules[j].mass * c_rot_arr[o2] / (rm * rot_rel_times.at(j, i)))); // 1100 // 13.75 = 55/4
              thermal_conductivity_LHS.at(p1 + 2, p2 + 1) = -6 * T / (n * D_cd * K_CONST_PI) * A_cd * n_ij * molecules[i].mass * molecules[i].mass * c_rot_arr[o1] / (rm * rot_rel_times.at(i, j) * (molecules[i].mass + molecules[j].mass)); // 0110 = 1001
              thermal_conductivity_LHS.at(p1    , p2 + 2) = 0; // 0001 = 0010
              thermal_conductivity_LHS.at(p1 + 1, p2 + 2) = -6 * T / (n * D_cd * K_CONST_PI) * A_cd * n_ij * molecules[j].mass * molecules[j].mass * c_rot_arr[o2] / (rm * rot_rel_times.at(j, i) * (molecules[i].mass + molecules[j].mass)); // 1001 = 0110
              thermal_conductivity_LHS.at(p1 + 2, p2 + 2) = 0; // 0011

              // symmetrization
              thermal_conductivity_LHS.at(p2    , p1    ) = thermal_conductivity_LHS.at(p1    , p2    );
              thermal_conductivity_LHS.at(p2    , p1 + 1) = thermal_conductivity_LHS.at(p1 + 1, p2    );
              thermal_conductivity_LHS.at(p2    , p1 + 2) = thermal_conductivity_LHS.at(p1 + 2, p2    );
              thermal_conductivity_LHS.at(p2 + 1, p1    ) = thermal_conductivity_LHS.at(p1    , p2 + 1);
              thermal_conductivity_LHS.at(p2 + 1, p1 + 1) = thermal_conductivity_LHS.at(p1 + 1, p2 + 1);
              thermal_conductivity_LHS.at(p2 + 1, p1 + 2) = thermal_conductivity_LHS.at(p1 + 2, p2 + 1);
              thermal_conductivity_LHS.at(p2 + 2, p1    ) = thermal_conductivity_LHS.at(p1    , p2 + 2);
              thermal_conductivity_LHS.at(p2 + 2, p1 + 1) = thermal_conductivity_LHS.at(p1 + 1, p2 + 2);
              thermal_conductivity_LHS.at(p2 + 2, p1 + 2) = thermal_conductivity_LHS.at(p1 + 2, p2 + 2);
            }
          }
        }
      }
    }

    // molecule + molecule collisions for identical molecular species at identical vibrational levels
    for (i=0; i<num_molecules; i++) {

      coll_mass = interactions[inter_index(i, i)].collision_mass;

      A_cd = ( 0.5 ) *      omega_22.at(i, i) / omega_11.at(i, i); // A_cc = 1
      B_cd = (1./3.) * (5 * omega_12.at(i, i) - omega_13.at(i, i)) / omega_11.at(i, i);
      C_cd = (1./3.) *      omega_12.at(i, i) / omega_11.at(i, i);
      D_cd = (3./8.) * kT/(n*molecules[i].mass* omega_11.at(i, i));
      rm   = ( 32. ) * n *  omega_22.at(i, i) / (5 * K_CONST_PI);
      D_cd_rot = D_cd;

      for (k=0; k<molecules[i].num_vibr_levels[0]; k++) {

        p1 = 3 * (vl_offset[i] + k);
        o1 = vl_offset[i] + k;

        n_ij = this_n_vl_mol[i][k] * this_n_vl_mol[i][k];

        thermal_conductivity_LHS.at(p1    , p1    ) = 0;	// 0000
        thermal_conductivity_LHS.at(p1 + 1, p1    ) = 0; 	// 1000
        thermal_conductivity_LHS.at(p1 + 2, p1    ) = 0; 	// 0010
        thermal_conductivity_LHS.at(p1    , p1 + 1) = 0; 	// 0100

        thermal_conductivity_LHS.at(p1 + 1, p1 + 1) = 1.5 * kT * n_ij * A_cd * (2 + (20./3) * molecules[i].mass * c_rot_arr[o1] / (K_CONST_K * K_CONST_PI * rm * rot_rel_times.at(i, i))) / (n * D_cd); // 1100

        thermal_conductivity_LHS.at(p1 + 2, p1 + 1) = -6 * T * A_cd * n_ij * molecules[i].mass * c_rot_arr[o1] / (K_CONST_PI * rm * rot_rel_times.at(i, i) * n * D_cd); // 0110

        thermal_conductivity_LHS.at(p1    , p1 + 2) = 0; // 0001

        thermal_conductivity_LHS.at(p1 + 1, p1 + 2) = -6 * T * A_cd * n_ij * molecules[i].mass * c_rot_arr[o1] / (K_CONST_PI * rm * rot_rel_times.at(i, i) * n * D_cd); // 1001

        thermal_conductivity_LHS.at(p1 + 2, p1 + 2) = T * n_ij * (molecules[i].mass * c_rot_arr[o1] / n) * (1.5 / D_cd_rot + 3.6 * A_cd / (K_CONST_PI * D_cd * rm * rot_rel_times.at(i, i))); // 0011 // 3.6=18/5
     }

     for (j=0; j<num_molecules; j++) {

       coll_mass = interactions[inter_index(i, j)].collision_mass;

       A_cd = 0.5 * omega_22.at(i, j) / omega_11.at(i, j);
       B_cd = (1./3.) * (5 * omega_12.at(i, j) - omega_13.at(i, j)) / omega_11.at(i, j);
       C_cd = (1./3.) * omega_12.at(i, j) / omega_11.at(i, j);
       D_cd = (3./16) * kT / (n * coll_mass * omega_11.at(i, j));
       D_cd_rot = D_cd;
       rm = 32 * n * omega_22.at(i, j) / (5 * K_CONST_PI);

       for (k=0; k<molecules[i].num_vibr_levels[0]; k++) {

         p1 = 3 * (vl_offset[i] + k);
         o1 =      vl_offset[i] + k ;

         for (l=0; l<molecules[j].num_vibr_levels[0]; l++) {
           if ((k != l) || (i != j)) {

             n_ij = this_n_vl_mol[i][k] * this_n_vl_mol[j][l];
             o2 = vl_offset[j] + l;

             // eq. 5.65
             thermal_conductivity_LHS.at(p1    , p1) += 1.5 * kT * n_ij / (D_cd * n); // 0000

             thermal_conductivity_LHS.at(p1 + 1, p1) -= 0.75 * kT * n_ij * molecules[j].mass * (6 * C_cd - 5) / ((molecules[i].mass + molecules[j].mass) * n * D_cd); // 1000

             thermal_conductivity_LHS.at(p1    , p1 + 1) -= 0.75 * kT * n_ij * molecules[j].mass * (6 * C_cd - 5) / ((molecules[i].mass + molecules[j].mass) * n * D_cd); // 0100

             thermal_conductivity_LHS.at(p1 + 1, p1 + 1) += 1.5 * kT * coll_mass * n_ij * (7.5 * molecules[i].mass / molecules[j].mass
                                                         + 6.25 * molecules[j].mass / molecules[i].mass - 3 * molecules[j].mass * B_cd / molecules[i].mass
                                                         + 4 * A_cd + (20./3) * A_cd * (molecules[i].mass * c_rot_arr[o1] / (rm * rot_rel_times.at(i, j))
                                                         + molecules[j].mass * c_rot_arr[o2] / (rm * rot_rel_times.at(j, i)))
                                                         / (K_CONST_PI * K_CONST_K)) / (n * D_cd * (molecules[i].mass + molecules[j].mass)); // 1100

             thermal_conductivity_LHS.at(p1 + 2, p1 + 1) -= 6. * T * A_cd * n_ij * molecules[i].mass * c_rot_arr[o1] / (rm * rot_rel_times.at(i, j)) * molecules[i].mass
                                                          / (n * D_cd * K_CONST_PI * (molecules[i].mass + molecules[j].mass)); // 0110

             thermal_conductivity_LHS.at(p1 + 1, p1 + 2) -= 6. * T * A_cd * n_ij * molecules[i].mass * c_rot_arr[o1] / (rm * rot_rel_times.at(i, j)) * molecules[i].mass
                                                          / (n * D_cd * K_CONST_PI * (molecules[i].mass + molecules[j].mass)); // 1001

             thermal_conductivity_LHS.at(p1 + 2, p1 + 2) += T * n_ij * molecules[i].mass * c_rot_arr[o1] * (1.5 / D_cd_rot + 3.6 * A_cd * molecules[i].mass
                                                          / (K_CONST_PI * D_cd * molecules[j].mass * rm * rot_rel_times.at(i, j))) / n; // 0011
           }
         }
       }
     }

     for (j=0; j<num_atoms; j++) {

       coll_mass = interactions[inter_index(i, num_molecules + j)].collision_mass;
       A_cd = 0.5 * omega_22.at(i, num_molecules + j) / omega_11.at(i, num_molecules + j);
       B_cd = (1./3.) * (5 * omega_12.at(i, num_molecules + j) - omega_13.at(i, num_molecules + j)) / omega_11.at(i, num_molecules + j);
       C_cd = (1./3.) * omega_12.at(i, num_molecules + j) / omega_11.at(i, num_molecules + j);
       D_cd = (3./16) * kT / (n * coll_mass * omega_11.at(i, num_molecules + j));
       D_cd_rot = D_cd;
       rm = 32 * n * omega_22.at(i, num_molecules + j) / (5 * K_CONST_PI);

       for (k=0; k<molecules[i].num_vibr_levels[0]; k++) {

         p1 = 3 * (vl_offset[i] + k);
         o1 = vl_offset[i] + k;
         n_ij = this_n_vl_mol[i][k] * this_n_atom[j];

         //std::cout << "rot_rel_times.at(MOL, AT)" << std::endl;
         //std::cout << rot_rel_times.at(i, num_molecules + j)  << std::endl;

         thermal_conductivity_LHS.at(p1    , p1)     += 1.5 * kT * n_ij / (D_cd * n); // 0000

         thermal_conductivity_LHS.at(p1 + 1, p1)     -= 0.75 * kT * n_ij * atoms[j].mass * (6 * C_cd - 5) / ((molecules[i].mass + atoms[j].mass) * n * D_cd); 	// 1000

         thermal_conductivity_LHS.at(p1    , p1 + 1) -= 0.75 * kT * n_ij * atoms[j].mass * (6 * C_cd - 5) / ((molecules[i].mass + atoms[j].mass) * n * D_cd); 	// 0100

         thermal_conductivity_LHS.at(p1 + 1, p1 + 1) += 1.5 * kT * coll_mass * n_ij *
   						                                          (7.5 * molecules[i].mass / atoms[j].mass +
 							                                          6.25 * atoms[j].mass / molecules[i].mass - 3 * atoms[j].mass * B_cd / molecules[i].mass
                                                        + 4 * A_cd + (20./3) * A_cd * (molecules[i].mass * c_rot_arr[o1] / (rm * rot_rel_times.at(i, num_molecules + j))) /
							                                          (K_CONST_PI * K_CONST_K)) / (n * D_cd * (molecules[i].mass + atoms[j].mass)); // 1100

         thermal_conductivity_LHS.at(p1 + 2, p1 + 1) -= 6. * T * A_cd * n_ij * molecules[i].mass * c_rot_arr[o1] / (rm * rot_rel_times.at(i, num_molecules + j))
                                                      * molecules[i].mass / (n * D_cd * K_CONST_PI * (molecules[i].mass + atoms[j].mass)); // 0110

         thermal_conductivity_LHS.at(p1 + 1, p1 + 2) -= 6. * T * A_cd * n_ij * molecules[i].mass * c_rot_arr[o1] / (rm * rot_rel_times.at(i, num_molecules + j)) * molecules[i].mass
                                                      / (n * D_cd * K_CONST_PI * (molecules[i].mass + atoms[j].mass)); // 0110

         thermal_conductivity_LHS.at(p1 + 2, p1 + 2) += T * n_ij * molecules[i].mass * c_rot_arr[o1] * (1.5 / D_cd_rot + 3.6 * A_cd * molecules[i].mass /
					                                    		    (K_CONST_PI * D_cd * atoms[j].mass * rm * rot_rel_times.at(i, num_molecules + j))) / n; // 0011
       }
      }
    }

    // molecule + atom collisions
    for (i=0; i<num_molecules; i++) {
      for (j=0; j<num_atoms; j++) {

        coll_mass = interactions[inter_index(i, num_molecules + j)].collision_mass;

        A_cd = 0.5 * omega_22.at(i, num_molecules + j) / omega_11.at(i, num_molecules + j);
        B_cd = (1./3.) * (5 * omega_12.at(i, num_molecules + j) - omega_13.at(i, num_molecules + j)) / omega_11.at(i, num_molecules + j);
        C_cd = (1./3.) * omega_12.at(i, num_molecules + j) / omega_11.at(i, num_molecules + j);
        D_cd = (3./16) * kT / (n * coll_mass * omega_11.at(i, num_molecules + j));
        rm = 32 * n * omega_22.at(i, num_molecules + j) / (5 * K_CONST_PI);

        for (k=0; k<molecules[i].num_vibr_levels[0]; k++) {

          n_ij = this_n_vl_mol[i][k] * this_n_atom[j];
          p1 = 3 * (vl_offset[i] + k);
          o1 =  vl_offset[i] + k;

          thermal_conductivity_LHS.at(p1, n_vibr_levels_total * 3 + 2 * j) = -1.5 * kT * n_ij / (n * D_cd); // 0000
          thermal_conductivity_LHS.at(p1 + 1, n_vibr_levels_total * 3 + 2 * j) =  0.75 * kT * n_ij * (6 * C_cd - 5) * atoms[j].mass / ((molecules[i].mass + atoms[j].mass) * n * D_cd); // 1000
          thermal_conductivity_LHS.at(p1 + 2, n_vibr_levels_total * 3 + 2 * j) = 0; // 0010
          thermal_conductivity_LHS.at(p1, n_vibr_levels_total*3 + 2*j + 1) = 0.75 * kT * n_ij * (6 * C_cd - 5) * molecules[i].mass / ((molecules[i].mass + atoms[j].mass) * n * D_cd); // 0100
          thermal_conductivity_LHS.at(p1 + 1, n_vibr_levels_total*3 + 2*j + 1) = -1.5 * kT * n_ij * coll_mass / (n * D_cd * (molecules[i].mass + atoms[j].mass))
                                                                                  * (13.75 - 3 * B_cd - 4 * A_cd - (20./3.) * A_cd / (K_CONST_K * K_CONST_PI)
                                                                                  * (molecules[i].mass * c_rot_arr[o1] / (rm * rot_rel_times.at(i, num_molecules + j)))); // 1100, 13.75 = 55/4
        /*
          thermal_conductivity_LHS.at(p1 + 2, n_vibr_levels_total*3 + 2*j + 1) = -6 * T / (n * D_cd * K_CONST_PI) * A_cd * n_ij * molecules[i].mass * molecules[i].mass
                                                                                  * c_rot_arr(o1) / (rm * rot_rel_times.at(i, num_molecules + j) * (molecules[i].mass + molecules[j].mass)); // 0110
         */
          // BugFix by Qizhen Hong
          thermal_conductivity_LHS.at(p1 + 2, n_vibr_levels_total*3 + 2*j + 1) = -6 * T / (n * D_cd * K_CONST_PI) * A_cd * n_ij * molecules[i].mass * molecules[i].mass
                                                                                  * c_rot_arr(o1) / (rm * rot_rel_times.at(i, num_molecules + j) * (molecules[i].mass + atoms[j].mass)); // 0110
          // symmetrization
          thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * j, p1    ) = thermal_conductivity_LHS.at(p1    , n_vibr_levels_total * 3 + 2 * j);
          thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * j, p1 + 1) = thermal_conductivity_LHS.at(p1 + 1, n_vibr_levels_total * 3 + 2 * j);
          thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * j, p1 + 2) = thermal_conductivity_LHS.at(p1 + 2, n_vibr_levels_total * 3 + 2 * j);
          thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * j + 1, p1    ) = thermal_conductivity_LHS.at(p1    , n_vibr_levels_total * 3 + 2 * j + 1);
          thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * j + 1, p1 + 1) = thermal_conductivity_LHS.at(p1 + 1, n_vibr_levels_total * 3 + 2 * j + 1);
          thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * j + 1, p1 + 2) = thermal_conductivity_LHS.at(p1 + 2, n_vibr_levels_total * 3 + 2 * j + 1);
        }
      }
    }

    // atoms + atom collisions for different atomic species
    for (i=0; i<num_atoms-1; i++) {
      for (j=i+1; j<num_atoms; j++) {

        coll_mass = interactions[inter_index(num_molecules + i, num_molecules + j)].collision_mass;

        n_ij = this_n_atom[i] * this_n_atom[j];

        A_cd = 0.5 * omega_22.at(num_molecules + i, num_molecules + j) / omega_11.at(num_molecules + i, num_molecules + j);
        B_cd = (1./3.) * (5 * omega_12.at(num_molecules + i, num_molecules + j) - omega_13.at(num_molecules + i, num_molecules + j)) / omega_11.at(num_molecules + i, num_molecules + j);
        C_cd = (1./3.) * omega_12.at(num_molecules + i, num_molecules + j) / omega_11.at(num_molecules + i, num_molecules + j);
        D_cd = (3./16) * kT / (n * coll_mass * omega_11.at(num_molecules + i, num_molecules + j));

        thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * i, n_vibr_levels_total * 3 + 2 * j) = -1.5 * kT * n_ij / (n * D_cd); // 0000
        thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * i + 1, n_vibr_levels_total * 3 + 2 * j) = 0.75 * kT * n_ij * (6 * C_cd - 5) * atoms[j].mass /
														((atoms[i].mass + atoms[j].mass) * n * D_cd); // 1000

        thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * i, n_vibr_levels_total * 3 + 2 * j + 1) =  0.75 * kT * n_ij * (6 * C_cd - 5) * atoms[i].mass / 																		((atoms[i].mass + atoms[j].mass) * n * D_cd);

        thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * i + 1, n_vibr_levels_total * 3 + 2 * j + 1) = -1.5 * kT * n_ij * coll_mass / (n * D_cd * (atoms[i].mass + atoms[j].mass))
                                                                                                                   * (13.75 - 3 * B_cd - 4 * A_cd); // 1100, 13.75 = 55/4

        // symmetrization
        thermal_conductivity_LHS.at(n_vibr_levels_total*3+2*j, n_vibr_levels_total*3+2*i) = thermal_conductivity_LHS.at(n_vibr_levels_total*3+2*i, n_vibr_levels_total*3+2*j);
        thermal_conductivity_LHS.at(n_vibr_levels_total*3+2*j, n_vibr_levels_total*3+2*i+1) = thermal_conductivity_LHS.at(n_vibr_levels_total*3+2*i+1, n_vibr_levels_total*3+2*j);
        thermal_conductivity_LHS.at(n_vibr_levels_total*3+2*j+1, n_vibr_levels_total*3+2*i) = thermal_conductivity_LHS.at(n_vibr_levels_total*3+2*i, n_vibr_levels_total*3+2*j+1);
        thermal_conductivity_LHS.at(n_vibr_levels_total*3+2*j+1, n_vibr_levels_total*3+2*i+1) = thermal_conductivity_LHS.at(n_vibr_levels_total*3+2*i+1, n_vibr_levels_total*3+2*j +1);
     }
    }

    // atom + atom collisions for identical atomic species
    for (i=0; i<num_atoms; i++) {

      coll_mass = interactions[inter_index(num_molecules + i, num_molecules + i)].collision_mass;
      n_ij = this_n_atom[i] * this_n_atom[i];

      A_cd = ( 0.5 ) *      omega_22.at(num_molecules + i, num_molecules + i) / omega_11.at(num_molecules + i, num_molecules + i);
      B_cd = (1./3.) * (5 * omega_12.at(num_molecules + i, num_molecules + i) - omega_13.at(num_molecules + i, num_molecules + i)) / omega_11.at(num_molecules + i, num_molecules + i);
      C_cd = (1./3.) *      omega_12.at(num_molecules + i, num_molecules + i) / omega_11.at(num_molecules + i, num_molecules + i);
      D_cd = (3./8.) *                                kT / (n * atoms[i].mass * omega_11.at(num_molecules + i, num_molecules + i));

      thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * i    , n_vibr_levels_total * 3 + 2 * i    ) = 0; 					// 0000
      thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * i + 1, n_vibr_levels_total * 3 + 2 * i    ) = 0; 					// 1000
      thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * i    , n_vibr_levels_total * 3 + 2 * i + 1) = 0; 					// 0100
      thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * i + 1, n_vibr_levels_total * 3 + 2 * i + 1) = 1.5 * kT * n_ij * 2 * A_cd / (n * D_cd);	// 1100

      for (j=0; j<num_molecules; j++) {

        coll_mass = interactions[inter_index(num_molecules + i, j)].collision_mass;

        A_cd = ( 0.5 ) *      omega_22.at(num_molecules + i, j) / omega_11.at(num_molecules + i, j);
        B_cd = (1./3.) * (5 * omega_12.at(num_molecules + i, j) - omega_13.at(num_molecules + i, j)) / omega_11.at(num_molecules + i, j);
        C_cd = (1./3.) *      omega_12.at(num_molecules + i, j) / omega_11.at(num_molecules + i, j);
        D_cd = (3./16) *                    kT / (n * coll_mass * omega_11.at(num_molecules + i, j));

        rm = 32 * n * omega_22.at(num_molecules + i, j) / (5 * K_CONST_PI);

        for (k=0; k<molecules[j].num_vibr_levels[0]; k++) {

          n_ij = this_n_atom[i] * this_n_vl_mol[j][k];
          // BugFix by Qizhen Hong
          //o1 =  vl_offset[i] + k;
          o1 =  vl_offset[j] + k;

          thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * i    , n_vibr_levels_total * 3 + 2 * i) +=  1.5 * kT * n_ij / (D_cd * n); // 0000
          thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * i + 1, n_vibr_levels_total * 3 + 2 * i) -= 	0.75 * kT * n_ij * molecules[j].mass * (6 * C_cd - 5) / ((atoms[i].mass + molecules[j].mass) * n * D_cd); // 1000
          thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * i    , n_vibr_levels_total * 3 + 2 * i + 1) -=  0.75 * kT * n_ij * molecules[j].mass * (6 * C_cd - 5) / ((atoms[i].mass + molecules[j].mass) * n * D_cd); // 0100
/*
          thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * i + 1, n_vibr_levels_total * 3 + 2 * i + 1) +=
		1.5 * kT * coll_mass * n_ij * (7.5 * atoms[i].mass / molecules[j].mass + 6.25 * molecules[j].mass / atoms[i].mass - 3 * molecules[j].mass * B_cd / atoms[i].mass
 		+ 4 * A_cd + (20./3) * A_cd * (molecules[j].mass * c_rot_arr[o1] / (rm * rot_rel_times.at(num_molecules + j, i))) / (K_CONST_PI * K_CONST_K)) /
               												(n * D_cd * (atoms[i].mass + molecules[j].mass)); // 1100
*/
          // BugFix by Qizhen Hong
          thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * i + 1, n_vibr_levels_total * 3 + 2 * i + 1) +=
		1.5 * kT * coll_mass * n_ij * (7.5 * atoms[i].mass / molecules[j].mass + 6.25 * molecules[j].mass / atoms[i].mass - 3 * molecules[j].mass * B_cd / atoms[i].mass
 		+ 4 * A_cd + (20./3) * A_cd * (molecules[j].mass * c_rot_arr[o1] / (rm * rot_rel_times.at(j,num_molecules + i))) / (K_CONST_PI * K_CONST_K)) /
               												(n * D_cd * (atoms[i].mass + molecules[j].mass)); // 1100
        }
      }

      for (j=0; j<num_atoms; j++) {
        if (j!=i) {

          coll_mass = interactions[inter_index(num_molecules + i, num_molecules + j)].collision_mass;
          n_ij = this_n_atom[i] * this_n_atom[j];
          A_cd = ( 0.5 ) *      omega_22.at(num_molecules + i, num_molecules + j) / omega_11.at(num_molecules + i, num_molecules + j);
          B_cd = (1./3.) * (5 * omega_12.at(num_molecules + i, num_molecules + j) - omega_13.at(num_molecules + i, num_molecules + j)) / omega_11.at(num_molecules + i, num_molecules + j);
          C_cd = (1./3.) *      omega_12.at(num_molecules + i, num_molecules + j) / omega_11.at(num_molecules + i, num_molecules + j);
          D_cd = (3./16) *      	 		             kT / (n * coll_mass * omega_11.at(num_molecules + i, num_molecules + j));

          thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * i    , n_vibr_levels_total * 3 + 2 * i) +=  1.5 * kT * n_ij / (D_cd * n); // 0000

          thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * i + 1, n_vibr_levels_total * 3 + 2 * i) -=  0.75 * kT * n_ij * atoms[j].mass * (6 * C_cd - 5) / ((atoms[i].mass + atoms[j].mass) * n * D_cd); // 1000

          thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * i    , n_vibr_levels_total * 3 + 2 * i + 1) -=  0.75 * kT * n_ij * atoms[j].mass * (6 * C_cd - 5) / ((atoms[i].mass + atoms[j].mass) * n * D_cd); // 0100

          thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + 2 * i + 1, n_vibr_levels_total * 3 + 2 * i + 1) += 1.5 * kT * coll_mass * n_ij *
				(7.5 * atoms[i].mass / atoms[j].mass + 6.25 * atoms[j].mass / atoms[i].mass - 3 * atoms[j].mass * B_cd / atoms[i].mass + 4 * A_cd) / ((atoms[i].mass + atoms[j].mass) * n * D_cd); // 1100
        }
      }
    }

    // division by n*n comes from the definition of mass molar fraction, see pp. 146
    thermal_conductivity_LHS /= (n * n);

    // condition so that the system is not singular
    // we divide by 10^40 so that the order of magnitude of the terms is roughly the same as the brackets
    for (i=0; i<num_molecules; i++) {
      for (k=0; k<molecules[i].num_vibr_levels[0]; k++) {

        p1 = 3 * (vl_offset[i] + k);

        thermal_conductivity_LHS.at(0, p1    ) = molecules[i].mass * this_n_vl_mol[i][k] / (n * 1e24);
        thermal_conductivity_LHS.at(0, p1 + 1) = 0;
        thermal_conductivity_LHS.at(0, p1 + 2) = 0;

      }
    }

    for (i=0; i<num_atoms; i++) {

      thermal_conductivity_LHS.at(0, n_vibr_levels_total * 3 + 2 * i    ) = atoms[i].mass * this_n_atom[i] / (n * 1e24);
      thermal_conductivity_LHS.at(0, n_vibr_levels_total * 3 + 2 * i + 1) = 0;
    }

    //std::cout << "T = " << T << std::endl;
    //std::cout << "thermal_conductivity_LHS" << thermal_conductivity_LHS << std::endl;
    return thermal_conductivity_LHS;

  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the right hand side of the algebraic system for thermal conductivity coefficients
  const arma::vec &kappa::Mixture::compute_thermal_conductivity_RHS(double T) {

    int i, k, j=0;

    for (i=0; i<num_molecules; i++) {
      for (k=0; k<molecules[i].num_vibr_levels[0]; k++) {

        // eqns. 5.16
        thermal_conductivity_RHS[3 * j    ] = 0;
        thermal_conductivity_RHS[3 * j + 1] = 7.5 * K_CONST_K * this_n_vl_mol[i][k];
        thermal_conductivity_RHS[3 * j + 2] = 3 * molecules[i].mass * this_n_vl_mol[i][k] * c_rot_arr[j];
        j++;
      }
    }

    for (i=0; i<num_atoms; i++) {
      thermal_conductivity_RHS[3 * n_vibr_levels_total + 2 * i    ] = 0;
      thermal_conductivity_RHS[3 * n_vibr_levels_total + 2 * i + 1] = 7.5 * K_CONST_K * this_n_atom[i];
    }

    thermal_conductivity_RHS *= T / this_total_n;

    //std::cout << "T = " << T << std::endl;
    //std::cout << "thermal_conductivity_RHS" << thermal_conductivity_RHS << std::endl;
    return thermal_conductivity_RHS;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // solve the algebraic system, computing the thermal conductivity coefficients in the STS approach
  const arma::vec &kappa::Mixture::compute_thermal_conductivity_coeffs(double T, kappa::models_omega model) {

    compute_thermal_conductivity_LHS(T, model);

    compute_thermal_conductivity_RHS(T);

    thermal_conductivity_coeffs = arma::solve(thermal_conductivity_LHS * 1e40, thermal_conductivity_RHS * 1e20) * 1e20;

    //std::cout << "thermal_conductivity_coeffs" << std::endl;
    //std::cout << thermal_conductivity_coeffs << std::endl;
    return thermal_conductivity_coeffs;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the left hand side of the algebraic system for thermal conductivity coefficients in the rigid rotator approximation
  const arma::mat &kappa::Mixture::compute_thermal_conductivity_rigid_rot_LHS(double T, kappa::models_omega model) {

    double n = this_total_n;
    double A_cd, B_cd, C_cd, D_cd, coll_mass, D_cd_rot, n_ij, rm, kT=K_CONST_K * T; // xi_rot = rm * tau_rot
    int i, j;

    // molecule + molecule collisions for different molecular species
    for (i=0; i<num_molecules - 1; i++) {
      for (j=i+1; j<num_molecules; j++) {

        coll_mass = interactions[inter_index(i, j)].collision_mass;
        n_ij = this_n_molecules[i] * this_n_molecules[j];
        A_cd = 0.5 * omega_22.at(i, j) / omega_11.at(i, j);
        B_cd = (1./3.) * (5 * omega_12.at(i, j) - omega_13.at(i, j)) / omega_11.at(i, j);
        C_cd = (1./3.) * omega_12.at(i, j) / omega_11.at(i, j);
        D_cd = (3./16) * kT / (n * coll_mass * omega_11.at(i, j));
        rm = 32 * n * omega_22.at(i, j) / (5 * K_CONST_PI);

        thermal_conductivity_rigid_rot_LHS.at(3*i,  3*j) = -1.5 * kT * n_ij / (n * D_cd); // 0000

        thermal_conductivity_rigid_rot_LHS.at(3*i+1,3*j) =  0.75 * kT * n_ij * (6 * C_cd - 5) * molecules[j].mass /
                                                            ((molecules[i].mass + molecules[j].mass) * n * D_cd); // 1000

        thermal_conductivity_rigid_rot_LHS.at(3*i+2,3*j) =  0; // 0010

        thermal_conductivity_rigid_rot_LHS.at(3*i,  3*j+1) = 0.75 * kT * n_ij * (6 * C_cd - 5) * molecules[i].mass /
                                                             ((molecules[i].mass + molecules[j].mass) * n * D_cd);

        thermal_conductivity_rigid_rot_LHS.at(3*i+1,3*j+1) = -1.5 * kT * n_ij * coll_mass /
                                                              (n * D_cd * (molecules[i].mass + molecules[j].mass))
                                                               * (13.75 - 3 * B_cd - 4 * A_cd - (20./3.) * A_cd
                                                               / (K_CONST_K * K_CONST_PI)
                                                               * (molecules[i].mass * c_rot_rigid_rot_arr[i]
                                                               / (rm * rot_rel_times.at(i, j))
                                                               + molecules[j].mass * c_rot_rigid_rot_arr[j] /
                                                               (rm * rot_rel_times.at(j, i)))); // 1100, 13.75 = 55/4

        thermal_conductivity_rigid_rot_LHS.at(3*i+2,3*j+1) = -6 * T / (n * D_cd * K_CONST_PI) * A_cd * n_ij
                                                              * molecules[i].mass * molecules[i].mass
                                                              * c_rot_rigid_rot_arr[i] / (rm * rot_rel_times.at(i, j)
                                                              * (molecules[i].mass + molecules[j].mass)); // 0110

        thermal_conductivity_rigid_rot_LHS.at(3*i,3*j+2) = 0; // 0001

        thermal_conductivity_rigid_rot_LHS.at(3*i+1,3*j+2) = -6 * T / (n * D_cd * K_CONST_PI) * A_cd * n_ij
                                                              * molecules[j].mass * molecules[j].mass
                                                              * c_rot_rigid_rot_arr[j] / (rm * rot_rel_times.at(j, i)
                                                              * (molecules[i].mass + molecules[j].mass)); // 1001

        thermal_conductivity_rigid_rot_LHS.at(3*i+2,3*j+2) = 0; // 0011

        // symmetrization
        thermal_conductivity_rigid_rot_LHS.at(3*j,3*i)     = thermal_conductivity_rigid_rot_LHS.at(3*i,3*j);
        thermal_conductivity_rigid_rot_LHS.at(3*j,3*i+1)   = thermal_conductivity_rigid_rot_LHS.at(3*i+1,3*j);
        thermal_conductivity_rigid_rot_LHS.at(3*j,3*i+2)   = thermal_conductivity_rigid_rot_LHS.at(3*i+2,3*j);
        thermal_conductivity_rigid_rot_LHS.at(3*j+1,3*i)   = thermal_conductivity_rigid_rot_LHS.at(3*i,3*j+1);
        thermal_conductivity_rigid_rot_LHS.at(3*j+1,3*i+1) = thermal_conductivity_rigid_rot_LHS.at(3*i+1,3*j+1);
        thermal_conductivity_rigid_rot_LHS.at(3*j+1,3*i+2) = thermal_conductivity_rigid_rot_LHS.at(3*i+2,3*j+1);
        thermal_conductivity_rigid_rot_LHS.at(3*j+2,3*i)   = thermal_conductivity_rigid_rot_LHS.at(3*i,3*j+2);
        thermal_conductivity_rigid_rot_LHS.at(3*j+2,3*i+1) = thermal_conductivity_rigid_rot_LHS.at(3*i+1,3*j+2);
        thermal_conductivity_rigid_rot_LHS.at(3*j+2,3*i+2) = thermal_conductivity_rigid_rot_LHS.at(3*i+2,3*j+2);
      }
    }

    for (i=0; i<num_molecules; i++) { // molecule + molecule collisions for identical molecular species

      coll_mass = interactions[inter_index(i, i)].collision_mass;
      A_cd = 0.5 * omega_22.at(i, i) / omega_11.at(i, i);
      B_cd = (1./3.) * (5 * omega_12.at(i, i) - omega_13.at(i, i)) / omega_11.at(i, i);
      C_cd = (1./3.) * omega_12.at(i, i) / omega_11.at(i, i);
      D_cd = (3./8) * kT / (n * molecules[i].mass * omega_11.at(i, i));
      D_cd_rot = D_cd;
      rm = 32 * n * omega_22.at(i, i) / (5 * K_CONST_PI);
      n_ij = this_n_molecules[i] * this_n_molecules[i];

      thermal_conductivity_rigid_rot_LHS.at(3 * i, 3 * i) = 0; // 0000
      thermal_conductivity_rigid_rot_LHS.at(3 * i + 1, 3 * i) =  0; // 1000
      thermal_conductivity_rigid_rot_LHS.at(3 * i + 2, 3 * i) = 0; // 0010
      thermal_conductivity_rigid_rot_LHS.at(3 * i, 3 * i + 1) = 0; // 0100

      thermal_conductivity_rigid_rot_LHS.at(3 * i + 1, 3 * i + 1) = 1.5 * kT * n_ij * A_cd * (2 + (20./3)
                                                                    * molecules[i].mass * c_rot_rigid_rot_arr[i]
                                                                    / (K_CONST_K * K_CONST_PI * rm * rot_rel_times.at(i, i)))
                                                                    / (n * D_cd); // 1100

      thermal_conductivity_rigid_rot_LHS.at(3 * i + 2, 3 * i + 1) = -6 * T * A_cd * n_ij * molecules[i].mass
                                                                    * c_rot_rigid_rot_arr[i]
                                                                    / (K_CONST_PI * rm * rot_rel_times.at(i, i) * n * D_cd); //0110

      thermal_conductivity_rigid_rot_LHS.at(3 * i, 3 * i + 2) = 0; // 0001
      thermal_conductivity_rigid_rot_LHS.at(3 * i + 1, 3 * i + 2) = -6 * T * A_cd * n_ij * molecules[i].mass
                                                                    * c_rot_rigid_rot_arr[i]
                                                                    / (K_CONST_PI * rm * rot_rel_times.at(i, i) * n * D_cd); //1001

      thermal_conductivity_rigid_rot_LHS.at(3*i+2,3*i+2) = T * n_ij * (molecules[i].mass * c_rot_rigid_rot_arr[i] / n)
                                                           * (1.5 / D_cd_rot + 3.6 * A_cd
                                                           / (K_CONST_PI * D_cd * rm * rot_rel_times.at(i, i))); // 0011, 3.6=18/5

      for (j=0; j<num_molecules; j++) {
        if (j != i) {

          coll_mass = interactions[inter_index(i, j)].collision_mass;
          n_ij = this_n_molecules[i] * this_n_molecules[j];
          A_cd = 0.5 * omega_22.at(i, j) / omega_11.at(i, j);
          B_cd = (1./3.) * (5 * omega_12.at(i, j) - omega_13.at(i, j)) / omega_11.at(i, j);
          C_cd = (1./3.) * omega_12.at(i, j) / omega_11.at(i, j);
          D_cd = (3./16) * kT / (n * coll_mass * omega_11.at(i, j));
          D_cd_rot = D_cd;
          rm = 32 * n * omega_22.at(i, j) / (5 * K_CONST_PI);

          thermal_conductivity_rigid_rot_LHS.at(3 * i, 3 * i) += 1.5 * kT * n_ij / (D_cd * n); // 0000
          thermal_conductivity_rigid_rot_LHS.at(3 * i + 1, 3 * i) -= 0.75 * kT * n_ij * molecules[j].mass * (6 * C_cd - 5) / ((molecules[i].mass + molecules[j].mass) * n * D_cd); // 1000
          // thermal_conductivity_rigid_rot_LHS.at(3 * i + 2, 3 * i) += 0.; // 0010

          thermal_conductivity_rigid_rot_LHS.at(3 * i, 3 * i + 1) -= 0.75 * kT * n_ij * molecules[j].mass * (6 * C_cd - 5) / ((molecules[i].mass + molecules[j].mass) * n * D_cd); // 0100
          thermal_conductivity_rigid_rot_LHS.at(3 * i + 1, 3 * i + 1) += 1.5 * kT * coll_mass * n_ij * (7.5 * molecules[i].mass / molecules[j].mass
                                                                                                        + 6.25 * molecules[j].mass / molecules[i].mass - 3 * molecules[j].mass * B_cd / molecules[i].mass
                                                                                                        + 4 * A_cd + (20./3) * A_cd * (molecules[i].mass * c_rot_rigid_rot_arr[i] / (rm * rot_rel_times.at(i, j))
                                                                                                                                       + molecules[j].mass * c_rot_rigid_rot_arr[j] / (rm * rot_rel_times.at(j, i)))
                                                                                                                             / (K_CONST_PI * K_CONST_K)) / (n * D_cd * (molecules[i].mass + molecules[j].mass)); // 1100
          thermal_conductivity_rigid_rot_LHS.at(3 * i + 2, 3 * i + 1) -= 6. * T * A_cd * n_ij * molecules[i].mass * c_rot_rigid_rot_arr[i] / (rm * rot_rel_times.at(i, j)) * molecules[i].mass
                                                                    / (n * D_cd * K_CONST_PI * (molecules[i].mass + molecules[j].mass)); // 0110

          // thermal_conductivity_rigid_rot_LHS.at(3 * i, 3 * i + 2) += 0.; // 0001
          thermal_conductivity_rigid_rot_LHS.at(3 * i + 1, 3 * i + 2) -= 6. * T * A_cd * n_ij * molecules[i].mass * c_rot_rigid_rot_arr[i] / (rm * rot_rel_times.at(i, j)) * molecules[i].mass
                                                                    / (n * D_cd * K_CONST_PI * (molecules[i].mass + molecules[j].mass)); // 1001
          thermal_conductivity_rigid_rot_LHS.at(3 * i + 2, 3 * i + 2) += T * n_ij * molecules[i].mass * c_rot_rigid_rot_arr[i] * (1.5 / D_cd_rot
                                                                                                                                  + 3.6 * A_cd * molecules[i].mass
                                                                                                                                     / (K_CONST_PI * D_cd * molecules[j].mass * rm * rot_rel_times.at(i, j))) / n; // 0011

        }
      }

      for (j=0; j<num_atoms; j++) {

        coll_mass = interactions[inter_index(i, num_molecules + j)].collision_mass;
        n_ij = this_n_molecules[i] * this_n_atom[j];
        A_cd = 0.5 * omega_22.at(i, num_molecules + j) / omega_11.at(i, num_molecules + j);
        B_cd = (1./3.) * (5 * omega_12.at(i, num_molecules + j) - omega_13.at(i, num_molecules + j)) / omega_11.at(i, num_molecules + j);
        C_cd = (1./3.) * omega_12.at(i, num_molecules + j) / omega_11.at(i, num_molecules + j);
        D_cd = (3./16) * kT / (n * coll_mass * omega_11.at(i, num_molecules + j));
        D_cd_rot = D_cd;
        rm = 32 * n * omega_22.at(i, num_molecules + j) / (5 * K_CONST_PI);

        thermal_conductivity_rigid_rot_LHS.at(3 * i, 3 * i) += 1.5 * kT * n_ij / (D_cd * n); // 0000
        thermal_conductivity_rigid_rot_LHS.at(3 * i + 1, 3 * i) -= 0.75 * kT * n_ij * atoms[j].mass * (6 * C_cd - 5) / ((molecules[i].mass + atoms[j].mass) * n * D_cd); // 1000
        // thermal_conductivity_rigid_rot_LHS.at(3 * i + 2, 3 * i) += 0.; // 0010

        thermal_conductivity_rigid_rot_LHS.at(3 * i, 3 * i + 1) -= 0.75 * kT * n_ij * atoms[j].mass * (6 * C_cd - 5) / ((molecules[i].mass + atoms[j].mass) * n * D_cd); // 0100
        thermal_conductivity_rigid_rot_LHS.at(3 * i + 1, 3 * i + 1) += 1.5 * kT * coll_mass * n_ij * (7.5 * molecules[i].mass / atoms[j].mass
                                                                                                   + 6.25 * atoms[j].mass / molecules[i].mass - 3 * atoms[j].mass * B_cd / molecules[i].mass
                                                                                                   + 4 * A_cd + (20./3) * A_cd * (molecules[i].mass * c_rot_rigid_rot_arr[i] / (rm * rot_rel_times.at(i, num_molecules + j)))
                                                                                                                        / (K_CONST_PI * K_CONST_K)) / (n * D_cd * (molecules[i].mass + atoms[j].mass)); // 1100

        thermal_conductivity_rigid_rot_LHS.at(3 * i + 2, 3 * i + 1) -= 6. * T * A_cd * n_ij * molecules[i].mass * c_rot_rigid_rot_arr[i] / (rm * rot_rel_times.at(i, num_molecules + j)) * molecules[i].mass
                                                                       / (n * D_cd * K_CONST_PI * (molecules[i].mass + atoms[j].mass)); // 0110

        // thermal_conductivity_rigid_rot_LHS.at(3 * i, 3 * i + 2) += 0.; // 0001
        thermal_conductivity_rigid_rot_LHS.at(3 * i + 1, 3 * i + 2) -= 6. * T * A_cd * n_ij * molecules[i].mass * c_rot_rigid_rot_arr[i] / (rm * rot_rel_times.at(i, num_molecules + j)) * molecules[i].mass
                                                                       / (n * D_cd * K_CONST_PI * (molecules[i].mass + atoms[j].mass)); // 0110
        thermal_conductivity_rigid_rot_LHS.at(3 * i + 2, 3 * i + 2) += T * n_ij * molecules[i].mass * c_rot_rigid_rot_arr[i] * (1.5 / D_cd_rot + 3.6 * A_cd * molecules[i].mass
                                                                                                                                                     / (K_CONST_PI * D_cd * atoms[j].mass
                                                                                                                                                     * rm * rot_rel_times.at(i, num_molecules + j))) / n; // 0011
      }
    }

    for (i=0; i<num_molecules; i++) { // molecule + atom collisions
      for (j=0; j<num_atoms; j++) {

        coll_mass = interactions[inter_index(i, num_molecules + j)].collision_mass;
        n_ij = this_n_molecules[i] * this_n_atom[j];
        A_cd = 0.5 * omega_22.at(i, num_molecules + j) / omega_11.at(i, num_molecules + j);
        B_cd = (1./3.) * (5 * omega_12.at(i, num_molecules + j) - omega_13.at(i, num_molecules + j)) / omega_11.at(i, num_molecules + j);
        C_cd = (1./3.) * omega_12.at(i, num_molecules + j) / omega_11.at(i, num_molecules + j);
        D_cd = (3./16) * kT / (n * coll_mass * omega_11.at(i, num_molecules + j));
        rm = 32 * n * omega_22.at(i, num_molecules + j) / (5 * K_CONST_PI);

        thermal_conductivity_rigid_rot_LHS.at(3 * i, num_molecules * 3 + 2 * j) = -1.5 * kT * n_ij / (n * D_cd); // 0000
        thermal_conductivity_rigid_rot_LHS.at(3 * i + 1, num_molecules * 3 + 2 * j) =  0.75 * kT * n_ij * (6 * C_cd - 5) * atoms[j].mass / ((molecules[i].mass + atoms[j].mass) * n * D_cd); // 1000
        thermal_conductivity_rigid_rot_LHS.at(3 * i + 2, num_molecules * 3 + 2 * j) = 0; // 0010


        thermal_conductivity_rigid_rot_LHS.at(3 * i, num_molecules * 3 + 2 * j + 1) = 0.75 * kT * n_ij * (6 * C_cd - 5) * molecules[i].mass / ((molecules[i].mass + atoms[j].mass) * n * D_cd); // 0100
        thermal_conductivity_rigid_rot_LHS.at(3 * i + 1, num_molecules * 3 + 2 * j + 1) = -1.5 * kT * n_ij * coll_mass / (n * D_cd * (molecules[i].mass + atoms[j].mass))
                                                                                             * (13.75 - 3 * B_cd - 4 * A_cd - (20./3.) * A_cd / (K_CONST_K * K_CONST_PI)
                                                                                                * (molecules[i].mass * c_rot_rigid_rot_arr[i] / (rm * rot_rel_times.at(i, num_molecules + j)))); // 1100, 13.75 = 55/4
        thermal_conductivity_rigid_rot_LHS.at(3 * i + 2, num_molecules * 3 + 2 * j + 1) = -6 * T / (n * D_cd * K_CONST_PI) * A_cd * n_ij * molecules[i].mass * molecules[i].mass
                                                                                           * c_rot_rigid_rot_arr[i] / (rm * rot_rel_times.at(i, num_molecules + j) * (molecules[i].mass + molecules[j].mass)); // 0110

        // symmetrization
        thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * j, 3 * i) = thermal_conductivity_rigid_rot_LHS.at(3 * i, num_molecules * 3 + 2 * j);
        thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * j, 3 * i + 1) = thermal_conductivity_rigid_rot_LHS.at(3 * i + 1, num_molecules * 3 + 2 * j);
        thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * j, 3 * i + 2) = thermal_conductivity_rigid_rot_LHS.at(3 * i + 2, num_molecules * 3 + 2 * j);

        thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * j + 1, 3 * i) = thermal_conductivity_rigid_rot_LHS.at(3 * i, num_molecules * 3 + 2 * j + 1);
        thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * j + 1, 3 * i + 1) = thermal_conductivity_rigid_rot_LHS.at(3 * i + 1, num_molecules * 3 + 2 * j + 1);
        thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * j + 1, 3 * i + 2) = thermal_conductivity_rigid_rot_LHS.at(3 * i + 2, num_molecules * 3 + 2 * j + 1);
      }
    }

    for (i=0; i<num_atoms-1; i++) { // atoms + atom collisions for different atomic species
      for (j=i+1; j<num_atoms; j++) {

        coll_mass = interactions[inter_index(num_molecules + i, num_molecules + j)].collision_mass;
        n_ij = this_n_atom[i] * this_n_atom[j];
        A_cd = 0.5 * omega_22.at(num_molecules + i, num_molecules + j) / omega_11.at(num_molecules + i, num_molecules + j);
        B_cd = (1./3.) * (5 * omega_12.at(num_molecules + i, num_molecules + j) - omega_13.at(num_molecules + i, num_molecules + j)) / omega_11.at(num_molecules + i, num_molecules + j);
        C_cd = (1./3.) * omega_12.at(num_molecules + i, num_molecules + j) / omega_11.at(num_molecules + i, num_molecules + j);
        D_cd = (3./16) * kT / (n * coll_mass * omega_11.at(num_molecules + i, num_molecules + j));

        thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * i, num_molecules * 3 + 2 * j) = -1.5 * kT * n_ij / (n * D_cd); // 0000
        thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * i + 1, num_molecules * 3 + 2 * j) = 0.75 * kT * n_ij * (6 * C_cd - 5) * atoms[j].mass / ((atoms[i].mass + atoms[j].mass) * n * D_cd); // 1000

        thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * i, num_molecules * 3 + 2 * j + 1) =  0.75 * kT * n_ij * (6 * C_cd - 5) * atoms[i].mass / ((atoms[i].mass + atoms[j].mass) * n * D_cd);
        thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * i + 1, num_molecules * 3 + 2 * j + 1) = -1.5 * kT * n_ij * coll_mass / (n * D_cd * (atoms[i].mass + atoms[j].mass))
                                                                                                               * (13.75 - 3 * B_cd - 4 * A_cd); // 1100, 13.75 = 55/4

        // symmetrization
        thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * j, num_molecules * 3 + 2 * i) = thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * i, num_molecules * 3 + 2 * j);
        thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * j, num_molecules * 3 + 2 * i + 1) = thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * i + 1, num_molecules * 3 + 2 * j);

        thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * j + 1, num_molecules * 3 + 2 * i) = thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * i, num_molecules * 3 + 2 * j + 1);
        thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * j + 1, num_molecules * 3 + 2 * i + 1) = thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * i + 1, num_molecules * 3 + 2 * j + 1);
      }
    }

    for (i=0; i<num_atoms; i++) { // atom + atom collisions for identical atomic species

      coll_mass = interactions[inter_index(num_molecules + i, num_molecules + i)].collision_mass;
      A_cd = 0.5 * omega_22.at(num_molecules + i, num_molecules + i) / omega_11.at(num_molecules + i, num_molecules + i);
      B_cd = (1./3.) * (5 * omega_12.at(num_molecules + i, num_molecules + i) - omega_13.at(num_molecules + i, num_molecules + i)) / omega_11.at(num_molecules + i, num_molecules + i);
      C_cd = (1./3.) * omega_12.at(num_molecules + i, num_molecules + i) / omega_11.at(num_molecules + i, num_molecules + i);
      D_cd = (3./8) * kT / (n * atoms[i].mass * omega_11.at(num_molecules + i, num_molecules + i));
      n_ij = this_n_atom[i] * this_n_atom[i];

      thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * i, num_molecules * 3 + 2 * i) = 0; // 0000
      thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * i + 1, num_molecules * 3 + 2 * i) =  0; // 1000

      thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * i, num_molecules * 3 + 2 * i + 1) = 0; // 0100
      thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * i + 1, num_molecules * 3 + 2 * i + 1) = 1.5 * kT * n_ij * 2 * A_cd / (n * D_cd); // 1100,


      for (j=0; j<num_molecules; j++) {

        coll_mass = interactions[inter_index(num_molecules + i, j)].collision_mass;
        n_ij = this_n_atom[i] * this_n_molecules[j];
        A_cd = 0.5 * omega_22.at(num_molecules + i, j) / omega_11.at(num_molecules + i, j);
        B_cd = (1./3.) * (5 * omega_12.at(num_molecules + i, j) - omega_13.at(num_molecules + i, j)) / omega_11.at(num_molecules + i, j);
        C_cd = (1./3.) * omega_12.at(num_molecules + i, j) / omega_11.at(num_molecules + i, j);
        D_cd = (3./16) * kT / (n * coll_mass * omega_11.at(num_molecules + i, j));
        rm = 32 * n * omega_22.at(num_molecules + i, j) / (5 * K_CONST_PI);

        thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * i, num_molecules * 3 + 2 * i) += 1.5 * kT * n_ij / (D_cd * n); // 0000
        thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * i + 1, num_molecules * 3 + 2 * i) -= 0.75 * kT * n_ij * molecules[j].mass * (6 * C_cd - 5) / ((atoms[i].mass + molecules[j].mass) * n * D_cd); // 1000

        thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * i, num_molecules * 3 + 2 * i + 1) -= 0.75 * kT * n_ij * molecules[j].mass * (6 * C_cd - 5) / ((atoms[i].mass + molecules[j].mass) * n * D_cd); // 0100
        thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * i + 1, num_molecules * 3 + 2 * i + 1) += 1.5 * kT * coll_mass * n_ij * (7.5 * atoms[i].mass / molecules[j].mass
                                                                                                   + 6.25 * molecules[j].mass / atoms[i].mass - 3 * molecules[j].mass * B_cd / atoms[i].mass
                                                                                                   + 4 * A_cd + (20./3) * A_cd * (molecules[j].mass * c_rot_rigid_rot_arr[j] / (rm * rot_rel_times.at(num_molecules + j, i)))
                                                                                                                        / (K_CONST_PI * K_CONST_K)) / (n * D_cd * (atoms[i].mass + molecules[j].mass)); // 1100
      }

      for (j=0; j<num_atoms; j++) {
        if (j!=i) {

          coll_mass = interactions[inter_index(num_molecules + i, num_molecules + j)].collision_mass;
          n_ij = this_n_atom[i] * this_n_atom[j];
          A_cd = 0.5 * omega_22.at(num_molecules + i, num_molecules + j) / omega_11.at(num_molecules + i, num_molecules + j);
          B_cd = (1./3.) * (5 * omega_12.at(num_molecules + i, num_molecules + j) - omega_13.at(num_molecules + i, num_molecules + j)) / omega_11.at(num_molecules + i, num_molecules + j);
          C_cd = (1./3.) * omega_12.at(num_molecules + i, num_molecules + j) / omega_11.at(num_molecules + i, num_molecules + j);
          D_cd = (3./16) * kT / (n * coll_mass * omega_11.at(num_molecules + i, num_molecules + j));

          thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * i, num_molecules * 3 + 2 * i) += 1.5 * kT * n_ij / (D_cd * n); // 0000
          thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * i + 1, num_molecules * 3 + 2 * i) -= 0.75 * kT * n_ij * atoms[j].mass * (6 * C_cd - 5) / ((atoms[i].mass + atoms[j].mass) * n * D_cd); // 1000

          thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * i, num_molecules * 3 + 2 * i + 1) -= 0.75 * kT * n_ij * atoms[j].mass * (6 * C_cd - 5) / ((atoms[i].mass + atoms[j].mass) * n * D_cd); // 0100
          thermal_conductivity_rigid_rot_LHS.at(num_molecules * 3 + 2 * i + 1, num_molecules * 3 + 2 * i + 1) += 1.5 * kT * coll_mass * n_ij * (7.5 * atoms[i].mass / atoms[j].mass
                                                                                                                                 + 6.25 * atoms[j].mass / atoms[i].mass - 3 * atoms[j].mass * B_cd / atoms[i].mass
                                                                                                                                 + 4 * A_cd) / (n * D_cd * (atoms[i].mass + atoms[j].mass)); // 1100
        }
      }
    }

    thermal_conductivity_rigid_rot_LHS /= (n * n);

    for (i=0; i<num_molecules; i++) { // condition so that the system is not singular; we divide by 10^40 so that the order of magnitude of the terms is roughly the same as the brackets
      thermal_conductivity_rigid_rot_LHS.at(0, 3 * i) = molecules[i].mass * this_n_molecules[i] / (this_total_n * 1e24);
      thermal_conductivity_rigid_rot_LHS.at(0, 3 * i + 1) = 0;
      thermal_conductivity_rigid_rot_LHS.at(0, 3 * i + 2) = 0;
    }
    for (i=0; i<num_atoms; i++) {
      thermal_conductivity_rigid_rot_LHS.at(0, num_molecules * 3 + 2 * i) = atoms[i].mass * this_n_atom[i] / (this_total_n * 1e24);
      thermal_conductivity_rigid_rot_LHS.at(0, num_molecules * 3 + 2 * i + 1) = 0;
    }

    return thermal_conductivity_rigid_rot_LHS;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the right hand side of the algebraic system for thermal conductivity coefficients in the rigid rotator approximation
  const arma::vec &kappa::Mixture::compute_thermal_conductivity_rigid_rot_RHS(double T) {

    int i;
    for (i=0; i<num_molecules; i++) {
      thermal_conductivity_rigid_rot_RHS[3 * i] = 0;
      thermal_conductivity_rigid_rot_RHS[3 * i + 1] = 7.5 * K_CONST_K * this_n_molecules[i];
      thermal_conductivity_rigid_rot_RHS[3 * i + 2] = 3 * molecules[i].mass * this_n_molecules[i] * c_rot_rigid_rot_arr[i];
    }

    for (i=0; i<num_atoms; i++) {
      thermal_conductivity_rigid_rot_RHS[num_molecules * 3 + 2 * i] = 0;
      thermal_conductivity_rigid_rot_RHS[num_molecules * 3 + 2 * i + 1] = 7.5 * K_CONST_K * this_n_atom[i];
    }

    thermal_conductivity_rigid_rot_RHS *= T / this_total_n;

    return thermal_conductivity_rigid_rot_RHS;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the thermal conductivity coefficients in the rigid rotator approximation
  const arma::vec &kappa::Mixture::compute_thermal_conductivity_rigid_rot_coeffs(double T, kappa::models_omega model) {

    compute_thermal_conductivity_rigid_rot_LHS(T, model);
    compute_thermal_conductivity_rigid_rot_RHS(T);

    thermal_conductivity_rigid_rot_coeffs = arma::solve(thermal_conductivity_rigid_rot_LHS * 1e40, thermal_conductivity_rigid_rot_RHS * 1e20) * 1e20;

    return thermal_conductivity_rigid_rot_coeffs;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // thermal conductivity, eq. 5.7
  double kappa::Mixture::thermal_conductivity(double T, kappa::models_omega model) {

    int i;
    double res=0;

    if (all_rigid_rotators) { // rigid rotator model
      compute_thermal_conductivity_rigid_rot_coeffs(T, model);
      for (i=0; i<num_molecules; i++) {
        res += this_n_molecules[i] * (1.25 * K_CONST_K * thermal_conductivity_rigid_rot_coeffs[3 * i + 1] + 0.50 * thermal_conductivity_rigid_rot_coeffs[3 * i + 2] * c_rot_rigid_rot_arr[i] * molecules[i].mass);
      }

      for (i=0; i<num_atoms; i++) {
        res += this_n_atom[i] * 1.25 * K_CONST_K * thermal_conductivity_rigid_rot_coeffs[num_molecules * 3 + 2 * i + 1];
      }
    } else { // non-rigid rotator model
        int k, j=0;
        compute_thermal_conductivity_coeffs(T, model);

        //#pragma omp parallel private(i, j, k) collapse(2)
        for (i=0; i<num_molecules; i++) {
          for (k=0; k<molecules[i].num_vibr_levels[0]; k++) {
            res += this_n_vl_mol[i][k] * (1.25 * K_CONST_K * thermal_conductivity_coeffs[3 * j + 1] + 0.50 * thermal_conductivity_coeffs[3 * j + 2] * c_rot_arr[j] * molecules[i].mass);
            j++;
          }
        }

        for (i=0; i<num_atoms; i++) {
          res += this_n_atom[i] * 1.25 * K_CONST_K * thermal_conductivity_coeffs[3 * n_vibr_levels_total + 2 * i + 1];
        }
    }

    return res / this_total_n;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // thermal diffusion, eq. 5.9
  void kappa::Mixture::thermodiffusion(double T, kappa::models_omega model) {

    int i, k, j=0;

    if (all_rigid_rotators) { // rigid rotator model
      for (i=0; i<num_molecules; i++) {
        for (k=0; k<molecules[i].num_vibr_levels[0]; k++) {
          th_diff[j] = thermal_conductivity_rigid_rot_coeffs[3 * i];
          j++;
        }
      }

      for (i=0; i<num_atoms; i++) {
        th_diff[j] = thermal_conductivity_rigid_rot_coeffs[num_molecules * 3 + 2 * i];
        j++;
      }
    } else { // non-rigid rotator model
      for (i=0; i<num_molecules; i++) {
        for (k=0; k<molecules[i].num_vibr_levels[0]; k++) {
          //th_diff[j] = thermal_conductivity_coeffs[3 * i];
          // BugFix by Qizhen Hong
          th_diff[j] = thermal_conductivity_coeffs[3 * j];
          j++;
        }
      }

      for (i=0; i<num_atoms; i++) {
        th_diff[j] = thermal_conductivity_coeffs[3 * n_vibr_levels_total + 2 * i];
        j++;
      }
    }
    th_diff = - 0.5 * th_diff / this_total_n;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the left hand side of the algebraic system for shear viscosity coefficients
  const arma::mat &kappa::Mixture::compute_shear_viscosity_LHS(double T, kappa::models_omega model) {

    int i1, i2;
    double tmp, kT = K_CONST_K * T;
    double n = this_total_n;

    // molecule + molecule collisions for different molecular species
    for (i1=0; i1<num_molecules-1; i1++) {
      for (i2=i1+1; i2<num_molecules; i2++) {

        // eq. 5.75 and eq. 5.62
        tmp = 3.2 * (this_n_molecules[i1] / n) * (this_n_molecules[i2] / n) *
		           (molecules[i1].mass * molecules[i2].mass) * omega_22.at(i1, i2) /
                     (kT * (molecules[i1].mass + molecules[i2].mass) * (molecules[i1].mass + molecules[i2].mass));

        tmp *= 1 - 10 * omega_11.at(i1, i2) / (3. * omega_22.at(i1, i2));

        shear_viscosity_LHS.at(i1, i2) = tmp;
        shear_viscosity_LHS.at(i2, i1) = tmp;
      }
    }
/*
    std::cout << " mass = " << molecules[0].mass  << std::endl;
    std::cout << "omega22 = " << omega_22 << std::endl;
    std::cout << "T = " << T << std::endl;
    std::cout << "shear_viscosity_LHS" << std::endl;
    std::cout << shear_viscosity_LHS << std::endl;
*/
    // molecule + molecule collisions for identical species
    for (i1=0; i1<num_molecules; i1++) {

      tmp = 1.6 * (this_n_molecules[i1] / n) * (this_n_molecules[i1] / n) * omega_22.at(i1, i1) / kT;  // 1.6 = 8/5

      for (i2=0; i2<num_molecules; i2++) {
        if (i2 != i1) {

          tmp += (3.2 * (this_n_molecules[i1] / n) * (this_n_molecules[i2] / n) * omega_22.at(i1, i2) / kT) *
		 (10. * omega_11.at(i1, i2) * molecules[i2].mass * molecules[i1].mass /
                 ((molecules[i2].mass + molecules[i1].mass) * (molecules[i2].mass + molecules[i1].mass) * 3. * omega_22.at(i1, i2)) +
                 molecules[i2].mass * molecules[i2].mass / ((molecules[i2].mass + molecules[i1].mass) * (molecules[i2].mass + molecules[i1].mass)));
        }
      }

      // loop over all atomic species
      for (i2=0; i2<num_atoms; i2++) {

        tmp += (3.2 * (this_n_molecules[i1] / n) * (this_n_atom[i2] / n) * omega_22.at(i1, num_molecules + i2) / kT)
               *  (10. * omega_11.at(i1, num_molecules + i2) / (3. * omega_22.at(i1, num_molecules + i2))
               * atoms[i2].mass * molecules[i1].mass / ((atoms[i2].mass + molecules[i1].mass) * (atoms[i2].mass + molecules[i1].mass))
               + atoms[i2].mass * atoms[i2].mass / ((atoms[i2].mass + molecules[i1].mass) * (atoms[i2].mass + molecules[i1].mass)));
      }

      shear_viscosity_LHS.at(i1, i1) = tmp;
    }

    // molecule + atom collisions
    for (i1=0; i1<num_molecules; i1++) {
      // loop over all atoms
      for (i2=0; i2<num_atoms; i2++) {

        tmp = 3.2 * (this_n_molecules[i1] / n) * (this_n_atom[i2] / n) * (molecules[i1].mass * atoms[i2].mass)
                  * omega_22.at(i1, num_molecules + i2) / (kT * (molecules[i1].mass + atoms[i2].mass) * (molecules[i1].mass + atoms[i2].mass)); // 3.2 = 16/5

        tmp *= 1. - 10. * omega_11.at(i1, num_molecules + i2) / (3. * omega_22.at(i1, num_molecules + i2));

        shear_viscosity_LHS.at(i1, num_molecules + i2) = tmp;
        shear_viscosity_LHS.at(num_molecules + i2, i1) = tmp;
      }
    }

    // atom + atom collisions for different atomic species
    for (i1=0; i1<num_atoms-1; i1++) {
      for (i2=i1+1; i2<num_atoms; i2++) {

        tmp = 3.2 * (this_n_atom[i1] / n) * (this_n_atom[i2] / n) * (atoms[i1].mass * atoms[i2].mass / (atoms[i1].mass + atoms[i2].mass))
                                         * omega_22.at(num_molecules + i1, num_molecules + i2) / (kT * (atoms[i1].mass + atoms[i2].mass)); // 3.2 = 16/5

        tmp *= 1 - 10 * omega_11.at(num_molecules + i1, num_molecules + i2) / (3. * omega_22.at(num_molecules + i1, num_molecules + i2));

        shear_viscosity_LHS.at(num_molecules + i1, num_molecules + i2) = tmp;
        shear_viscosity_LHS.at(num_molecules + i2, num_molecules + i1) = tmp;
     }
    }

    // atom + atom collisions for identical species
    for (i1=0; i1<num_atoms; i1++) {

      tmp = 1.6 * (this_n_atom[i1] / n) * (this_n_atom[i1] / n) * omega_22.at(num_molecules + i1, num_molecules + i1) / kT;  // 1.6 = 8/5

      // loop over all molecular species
      for (i2=0; i2<num_molecules; i2++) {

        tmp += (3.2 * (this_n_atom[i1] / n) * (this_n_molecules[i2] / n) * omega_22.at(i2, num_molecules + i1) / kT)
                * (10. * omega_11.at(i2, num_molecules + i1) / (3. * omega_22.at(i2, num_molecules + i1))
                * molecules[i2].mass * atoms[i1].mass / ((molecules[i2].mass + atoms[i1].mass) * (molecules[i2].mass + atoms[i1].mass))
                + molecules[i2].mass * molecules[i2].mass / ((molecules[i2].mass + atoms[i1].mass) * (molecules[i2].mass + atoms[i1].mass)));
      }

      // loop over all atomic species
      for (i2=0; i2<num_atoms; i2++) {
        if (i2 != i1) {
          tmp += (3.2 * (this_n_atom[i1] / n) * (this_n_atom[i2] / n) * omega_22.at(num_molecules + i2, num_molecules + i1) / kT)
              *  (10. * omega_11.at(num_molecules + i2, num_molecules + i1) / (3. * omega_22.at(num_molecules + i2, num_molecules + i1))
              * atoms[i2].mass * atoms[i1].mass / ((atoms[i2].mass + atoms[i1].mass) * (atoms[i2].mass + atoms[i1].mass))
              + atoms[i2].mass * atoms[i2].mass / ((atoms[i2].mass + atoms[i1].mass) * (atoms[i2].mass + atoms[i1].mass)));
        }
      }

      shear_viscosity_LHS.at(num_molecules + i1, num_molecules + i1) = tmp;
    }

    //std::cout << "T = " << T << std::endl;
    //std::cout << "shear_viscosity_LHS" << shear_viscosity_LHS << std::endl;
    return shear_viscosity_LHS;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the right hand side of the algebraic system for shear viscosity coefficients
  const arma::vec &kappa::Mixture::compute_shear_viscosity_RHS(double T) {

    shear_viscosity_RHS.zeros();

    int i;
    for (i = 0; i < num_molecules; i++) {
      shear_viscosity_RHS[i] = this_n_molecules[i];
    }

    for (i = 0; i < num_atoms; i++) {
      shear_viscosity_RHS[num_molecules + i] = this_n_atom[i];
    }

    shear_viscosity_RHS *= 2 / (K_CONST_K * T * this_total_n);
/*
    std::cout << "this_n_molecules[i]/n = " << this_n_molecules[0]/this_total_n << std::endl;
    std::cout << "this_n_atom[i]/n = " << this_n_atom[0]/this_total_n << std::endl;
    std::cout << "K_CONST_K = " << K_CONST_K << std::endl;
    std::cout << "T = " << T << std::endl;
    std::cout << "2kT = " << 2 / (K_CONST_K * T) << std::endl;
    std::cout << "shear_viscosity_RHS" << std::endl;
    std::cout << shear_viscosity_RHS << std::endl;
    std::cout << shear_viscosity_RHS[0] + shear_viscosity_RHS[1] << std::endl;
*/
    //std::cout << "T = " << T << std::endl;
    //std::cout << "shear_viscosity_RHS" << shear_viscosity_RHS << std::endl;
    return shear_viscosity_RHS;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the shear viscosity coefficients
  const arma::vec &kappa::Mixture::compute_shear_viscosity_coeffs(double T, kappa::models_omega model) {

    compute_shear_viscosity_LHS(T, model);
    compute_shear_viscosity_RHS(T);

    shear_viscosity_coeffs = arma::solve(shear_viscosity_LHS * 1e20, shear_viscosity_RHS) * 1e20;

    return shear_viscosity_coeffs;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // shear viscosity, eq. 5.10
  double kappa::Mixture::shear_viscosity(double T, kappa::models_omega model) {

    compute_shear_viscosity_coeffs(T, model);

    int i;
    double res=0;

    for (i = 0; i < num_molecules; i++) {
      res += this_n_molecules[i] * shear_viscosity_coeffs[i];
    }

    for (i = 0; i < num_atoms; i++) {
      res += this_n_atom[i] * shear_viscosity_coeffs[num_molecules + i];
    }

    return res * K_CONST_K * T / (2 * this_total_n);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the left hand side of the algebraic system for diffusion coefficients
  const arma::mat &kappa::Mixture::compute_diffusion_LHS(double T) {

    int i, k, offset=0;

    // we just need to copy data from the thermal conductivity LHS matrix
    for (i=0; i<n_vibr_levels_total; i++) {
      for (k=0; k<n_vibr_levels_total; k++) {
        diffusion_LHS(i, k) = thermal_conductivity_LHS.at(i * 3, k * 3);
      }

      for (k=0; k<num_atoms; k++) {
        diffusion_LHS(i, n_vibr_levels_total + k) = thermal_conductivity_LHS.at(i * 3, n_vibr_levels_total * 3 + k * 2);
        diffusion_LHS(n_vibr_levels_total + k, i) = thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + k * 2, i * 3);
      }
    }

    for (i=0; i<num_atoms; i++) {
      for (k=0; k<num_atoms; k++) {
        diffusion_LHS(n_vibr_levels_total + i, n_vibr_levels_total + k) =
                                    thermal_conductivity_LHS.at(n_vibr_levels_total * 3 + i * 2, n_vibr_levels_total * 3 + k * 2);
      }
    }

    // now, auxiliary condition - first row of matrix is replaced
    for (i=0; i<num_molecules; i++) {
      for (k=0; k<molecules[i].num_vibr_levels[0]; k++) {
        diffusion_LHS.at(0, offset) = this_n_vl_mol[i][k] * molecules[i].mass / (this_total_dens * 1e41);
        //diffusion_LHS.at(0, offset) = this_n_vl_mol[i][k] * molecules[i].mass / (this_total_dens);
        offset += 1;
      }
    }

    for (i=0; i<num_atoms; i++) {
      diffusion_LHS.at(0, offset) = this_n_atom[i] * atoms[i].mass / (this_total_dens * 1e41);
      //diffusion_LHS.at(0, offset) = this_n_atom[i] * atoms[i].mass / (this_total_dens);
      offset += 1;
    }

    //std::cout << "diffusion_LHS" << std::endl;
    //std::cout << diffusion_LHS << std::endl;

    return diffusion_LHS;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the left hand side of the simplified algebraic system for diffusion coefficients
  // const arma::mat &kappa::Mixture::lite_compute_diffusion_LHS(double T, int b, int n, std::string ParticleType) {}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the right hand side of the algebraic system for diffusion coefficients
  const arma::vec &kappa::Mixture::compute_diffusion_RHS(double T, int b, int n) {

    /*
    Computes RHS for diffusion
    from eq. 5.17, for fixed (c,i,b,n) we have
    rhs(c,i,b,n) = 3kT * (delta(c,b) * delta(i,n) - rho_ci / rho)

    We also need an auxiliary condition for each system:
    sum over (c,i) of rho_ci * d^bn_ci,0 / rho = 0

    So this computes the RHS for a fixed b (# of particle) and a fixed n (# of vibrational level)
    */

    int i, k, offset=0;

    for (i=0; i<num_molecules; i++) { // molecules
      for (k=0; k<molecules[i].num_vibr_levels[0]; k++) { // vibr. levels
	      // BugFix by Anna Becina:
        // this seems to fix the asymmetricity of transport coeffs. matrix
        //if (i==b || k==n) {
	      if (i==b && k==n) {
          diffusion_RHS.at(offset) += 1;
        }
        diffusion_RHS.at(offset) -= this_n_vl_mol[i][k] * molecules[i].mass / this_total_dens;
        offset += 1;
      }
    }

    for (i=0; i<num_atoms; i++) { // over atoms
      if ((num_molecules + i)==b) {
        diffusion_RHS.at(offset) += 1;
      }
      diffusion_RHS.at(offset) -= this_n_atom[i] * atoms[i].mass / this_total_dens;
      offset += 1;
    }

    diffusion_RHS *= 3 * K_CONST_K * T;
    diffusion_RHS[0] = 0; // auxiliary condition

    //std::cout << "diffusion_RHS" << std::endl;
    //std::cout << diffusion_RHS << std::endl;

    return diffusion_RHS;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the right hand side of the simplified algebraic system for diffusion coefficients
  // const arma::vec &kappa::Mixture::lite_compute_diffusion_RHS(double T, int b, int n, std::string ParticleType) {}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the diffusion coefficients
  const arma::vec &kappa::Mixture::compute_diffusion_coeffs(double T, int b, int n) {

    // compute_diffusion_LHS(T); - the diffusion matrix does not change, so we compute it before this function!
    compute_diffusion_RHS(T, b, n);

    // diffusion_coeffs = arma::solve(diffusion_LHS * 1e40, diffusion_RHS * 1e20) * 1e20;
    // diffusion_coeffs = arma::solve(diffusion_LHS * 1e40, diffusion_RHS * 1e20, arma::solve_opts::no_approx);
    // diffusion_coeffs = arma::solve(diffusion_LHS * 1e40, diffusion_RHS * 1e20, arma::solve_opts::fast);
    // diffusion_coeffs = arma::solve(diffusion_LHS * 1e40, diffusion_RHS * 1e20, arma::solve_opts::equilibrate);
    diffusion_coeffs = arma::solve(diffusion_LHS, diffusion_RHS);

    return diffusion_coeffs;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the diffusion coefficients
  // const arma::vec &kappa::Mixture::lite_compute_diffusion_coeffs(double T, int b, int n, std::string ParticleType) {}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void kappa::Mixture::diffusion(double T) {

    compute_diffusion_LHS(T);

    int i, k, offset=0;
    for (i=0; i<num_molecules; i++) { // molecules, b in eq. 5.19
      for (k=0; k<molecules[i].num_vibr_levels[0]; k++) { // vibr. levels, n
        diff.col(offset) = compute_diffusion_coeffs(T, i, k);
        //std::cout << "diffusion_coeffs at offset = " << offset << std::endl;
        //std::cout << compute_diffusion_coeffs(T, i, k) << std::endl;
        offset += 1;
      }
    }
    for (i=0; i<num_atoms; i++) {
      diff.col(offset) = compute_diffusion_coeffs(T, num_molecules + i, 0);
      offset += 1;
    }

    diff /= (2 * this_total_n); // D_cidk of eq. 5.8
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // void kappa::Mixture::lite_diffusion(double T) {}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Diffusion coefficients for binary and diatomic mixtures
  // WATCH OUT! Coefficients are stored in the following order: DA2iA2i | DA2iA2k | DA2A | DAA
  void kappa::Mixture::binary_diffusion(double T, kappa::models_omega model) {

    int i, j, c, k, d;
    double n      = this_total_n;
    double rho    = compute_density(this_n_vl_mol, this_n_atom);
    double rho_A2 = compute_density(this_n_vl_mol);
    double rho_A  = this_n_atom[0] * atoms[0].mass;
    double Y_A2   = rho_A2/rho;
    double Y_A    = rho_A/rho;
    double coll_mass_A_A2 = interactions[inter_index(0, num_molecules)].collision_mass;
    double coll_mass_A2_A2 = interactions[inter_index(0, 0)].collision_mass;
    double m_A2   = molecules[0].mass;
    double m_A    = atoms[0].mass;
    double D_AA2  = (3./16) * (K_CONST_K * T)/(n * coll_mass_A_A2 * omega_11.at(0, num_molecules));
    double D_A2   = (3./8.) * (K_CONST_K * T)/(n * m_A2 * omega_11.at(0, 0));

    binary_diff.zeros(n_vibr_levels_total+3);

    if (num_atoms==0) { // diatomic molecule mixture, pp. 158 Nagnibeda & Kustova, 2009.

      for (k=0; k<molecules[0].num_vibr_levels[0]; k++) {
        binary_diff[k] = D_A2 * (1./((molecules[0].mass * this_n_vl_mol[0][k])/rho) -1.); // DA2iA2i
      }
      binary_diff[molecules[0].num_vibr_levels[0]] = -D_A2; // DA2A2
    } else { // binary mixture, eq. 5.108-5.109
      for (k=0; k<molecules[0].num_vibr_levels[0]; k++) {

        // molecule + molecule at the same vibrational level, DA2iA2i
        // TODO bring constant terms out of loop
        binary_diff[k] = D_AA2 * (m_A/(rho/n)) * (m_A/(rho/n)) * ( (Y_A/D_A2) + (2./D_AA2) *
                         ( (1./((molecules[0].mass * this_n_vl_mol[0][k])/rho) - Y_A - 1. ))) / ( (Y_A2/(2.*D_A2)) + (Y_A/D_AA2) );
      }

      // molecule + molecule at different vibrational level, DA2A2 = DA2iA2k
      binary_diff[molecules[0].num_vibr_levels[0]] = D_AA2*(m_A/(rho/n))*(m_A/(rho/n))*
                                                     ((Y_A/D_A2)-(2./D_AA2)*(Y_A+1.))/((Y_A2/(2.*D_A2))+(Y_A/D_AA2));

      // molecule + atom, DAA2
      binary_diff[molecules[0].num_vibr_levels[0] + 1] = -D_AA2 * (m_A2 * m_A) / ((rho/n) * (rho/n));

      // atom + atom, DAA
      binary_diff[molecules[0].num_vibr_levels[0] + 2] = D_AA2 * ( (m_A2 * m_A) / ((rho/n) * (rho/n))) * ((1./Y_A) -1.);
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // general function to compute all transport coefficients
  void kappa::Mixture::compute_transport_coefficients(double T, const std::vector<arma::vec> &n_vl_molecule,
                                                      const arma::vec &n_atom, double n_electrons, kappa::models_omega model,
                                                      double perturbation)
  {

    #ifdef KAPPA_STRICT_CHECKS
    std::string error_string;
    if (T <= 0) {
      error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (perturbation < 0) {
      error_string = "Negative perturbation: perturbation=" + std::to_string(perturbation);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    check_n_vl_molecule(n_vl_molecule);
    check_n_atom(n_atom);
    if (n_electrons < 0) {
      error_string = "Number density of electrons is negative: n_electrons=" + std::to_string(n_electrons);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    cache_on = true;

    this_n_atom = n_atom;
    this_n_vl_mol = n_vl_molecule;
    this_n_molecules = compute_n_molecule(n_vl_molecule);
    this_total_n = compute_n(this_n_molecules, n_atom, n_electrons);
    this_total_dens = compute_density(this_n_molecules, n_atom, n_electrons);

    int i, j, addcounter=0;

    if (is_ionized) {
      addcounter = 1;
      this_n_electrons = (1 - perturbation) * n_electrons + this_total_n * perturbation /
                         (n_vibr_levels_total + num_atoms + addcounter);
    }

    for (i=0; i<num_molecules; i++) {
      for (j=0; j<molecules[i].num_vibr_levels[0]; j++) {
        this_n_vl_mol[i][j] = (1 - perturbation) * n_vl_molecule[i][j] + this_total_n * perturbation /
                              (n_vibr_levels_total + num_atoms + addcounter);
      }
    }

    for (i=0; i<num_atoms; i++) {
      this_n_atom[i] = (1 - perturbation) * n_atom[i] + this_total_n * perturbation /
                       (n_vibr_levels_total + num_atoms + addcounter);
    }

    // Check that total number density is conserved
    double n_dens_cons=0.;
    for (i=0; i<num_molecules; i++) {
      for (j=0; j<molecules[i].num_vibr_levels[0]; j++) {
        n_dens_cons += this_n_vl_mol[i][j];
      }
    }

    for (i=0; i<num_atoms; i++) {
      n_dens_cons += this_n_atom[i];
    }
    std::string error_string;
    if ((n_dens_cons/this_total_n) < 0.99 || (n_dens_cons/this_total_n) > 1.01) {
      error_string = "Total number density is NOT conserved!" + std::to_string(n_dens_cons/this_total_n);
      throw kappa::IncorrectValueException(error_string.c_str());
    }

    if (all_rigid_rotators) {
      compute_c_rot_rigid_rot(T);
      compute_full_crot_rigid_rot(T);
    } else {
      compute_c_rot(T);
      compute_full_crot(T);
    }
    this_ctr = c_tr(this_n_molecules, n_atom, n_electrons);

    // computation of bracket integrals
    compute_omega11(T, model);
    compute_omega12(T, model);
    compute_omega13(T, model);
    compute_omega22(T, model);

    if (is_ionized) {
      double dl = debye_length(T, n_vl_molecule, n_atom, n_electrons);
      compute_omega11(T, dl, model);
      compute_omega12(T, dl, model);
      compute_omega13(T, dl, model);
      compute_omega22(T, dl, model);
    }

    // computation of rotational relaxation times
    compute_rot_rel_times(T, this_total_n, model);

    // extraction of transport coefficients' values
    th_cond = thermal_conductivity(T, model);
    sh_visc = shear_viscosity(T, model);
    b_visc = bulk_viscosity(T, model);
    thermodiffusion(T, model);
    diffusion(T);
    // lite_diffusion(T);
    binary_diffusion(T, model);
    cache_on = false;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void kappa::Mixture::compute_transport_coefficients(double T, const std::vector<arma::vec> &n_vl_molecule,
                                                      const arma::vec &n_atom, kappa::models_omega model, double perturbation)
  {
    compute_transport_coefficients(T, n_vl_molecule, n_atom, 0.0, model, perturbation);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void kappa::Mixture::compute_transport_coefficients(double T, const std::vector<arma::vec> &n_vl_molecule, double n_electrons,
                                                      kappa::models_omega model, double perturbation)
  {
    compute_transport_coefficients(T, n_vl_molecule, empty_n_atom, n_electrons, model, perturbation);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void kappa::Mixture::compute_transport_coefficients(double T, const std::vector<arma::vec> &n_vl_molecule,
                                                      kappa::models_omega model, double perturbation)
  {
    compute_transport_coefficients(T, n_vl_molecule, empty_n_atom, 0.0, model, perturbation);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void kappa::Mixture::compute_transport_coefficients(double T, const arma::vec &n, kappa::models_omega model, double perturbation)
  {

    // TODO: Test function!
    #ifdef KAPPA_STRICT_CHECKS
    std::string error_string;
    if (T <= 0) {
      error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (perturbation < 0) {
      error_string = "Negative perturbation: perturbation=" + std::to_string(perturbation);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (!all_rigid_rotators) {
      error_string = "Molecules in mixture are not rigid rotators, must pass vector of vibrational level populations";
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    check_n(n);
    #endif

    cache_on = true;

    if (num_molecules > 0) {
      this_n_molecules = n.subvec(0, num_molecules-1);
    }
    if (num_atoms > 0) {
      this_n_atom = n.subvec(num_molecules, num_molecules+num_atoms-1);
    }
    if (is_ionized) {
      this_n_electrons = n[n.n_elem - 1];
    }

    this_total_n = compute_n(n);
    this_total_dens = compute_density(n);

    int i, addcounter=0;

    if (is_ionized) {
      addcounter = 1;
      this_n_electrons = (1 - perturbation) * this_n_electrons + this_total_n * perturbation /
                         (num_molecules + num_atoms + addcounter);
    }

    for (i=0; i<num_molecules; i++) {
      this_n_molecules[i] = (1 - perturbation) *  this_n_molecules[i]  + this_total_n * perturbation /
                            (num_molecules + num_atoms + addcounter);
    }

    for (i=0; i<num_atoms; i++) {
      this_n_atom[i] = (1 - perturbation) * this_n_atom[i] + this_total_n * perturbation /
                       (num_molecules + num_atoms + addcounter);
    }

    // TODO: Check that total number density is conserved!

    compute_c_rot_rigid_rot(T);
    compute_full_crot_rigid_rot(T);

    this_ctr = c_tr(n);

    compute_omega11(T, model);
    compute_omega12(T, model);
    compute_omega13(T, model);
    compute_omega22(T, model);

    if (is_ionized) {
      double dl = debye_length(T, n);
      compute_omega11(T, dl, model);
      compute_omega12(T, dl, model);
      compute_omega13(T, dl, model);
      compute_omega22(T, dl, model);
    }

    compute_rot_rel_times(T, this_total_n, model);

    th_cond = thermal_conductivity(T, model);
    sh_visc = shear_viscosity(T, model);
    b_visc = bulk_viscosity(T, model);
    thermodiffusion(T, model);
    diffusion(T);
    // lite_diffusion(T);
    binary_diffusion(T, model);
    cache_on = false;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Mixture::get_thermal_conductivity() {
    return th_cond;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Mixture::get_shear_viscosity() {
    return sh_visc;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Mixture::get_bulk_viscosity() {
    return b_visc;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  arma::vec kappa::Mixture::get_thermodiffusion() {
    return th_diff;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  arma::mat kappa::Mixture::get_diffusion() {
    return diff;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  arma::mat kappa::Mixture::get_lite_diffusion() {
    return lite_diff;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  arma::vec kappa::Mixture::get_binary_diffusion() {
    return binary_diff;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
