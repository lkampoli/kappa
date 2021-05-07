/*!
   \file ApproximationOneT.cpp
 
   \brief Chapter 4
   \brief Implementation of the One-Temperature Model for Chemically Non-equilibrium Gas Mixtures.
 
   \see Class ApproximationOneT
 */

#include "approximation_onet.hpp"

kappa::ApproximationOneT::ApproximationOneT():kappa::Approximation() {} 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationOneT::Z_int(double T, const arma::vec &electron_energy, const arma::Col<unsigned long> &statistical_weight, int n_electron_levels) {

    double res=0.0;
    for (int e=0; e<n_electron_levels; e++) {
      res += statistical_weight[e] * exp(-electron_energy[e] / (K_CONST_K * T));
    }
    return res;
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationOneT::Z_int(double T, 
					 const arma::vec &electron_energy, 
					 const arma::Col<unsigned long> &statistical_weight, 
					 int n_electron_levels,
					 const std::vector<arma::vec> &vibr_energy, 
					 const std::vector<int> &num_vibr_levels,
					 const std::vector<std::vector<arma::vec>> &rot_energy, 
					 const std::vector<std::vector<int>> &num_rot_levels, 
					 int rot_symmetry) {

    double res=0.0;
    for (int e=0; e<n_electron_levels; e++) {
      for (int i=0; i<num_vibr_levels[e]; i++) {
        res += statistical_weight[e] * arma::dot(2*arma::linspace<arma::vec>(0, num_rot_levels[e][i]-1, num_rot_levels[e][i]) + 1,
                                       arma::exp(-(electron_energy[e] + vibr_energy[e][i] + rot_energy[e][i]) / (K_CONST_K * T)));
      }
    }
    return res / rot_symmetry;
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Internal partition function for an atom
  double kappa::ApproximationOneT::Z_int(double T, const kappa::Atom &atom, double Delta_E) {

    if (Delta_E == -1) {
      return atom.statistical_weight[0];
    } else {
      return Z_int(T, atom.electron_energy, atom.statistical_weight, p_max_electron_level(atom.electron_energy, atom.ionization_potential, Delta_E));
    }
  }
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Internal partition function for a molecule
  double kappa::ApproximationOneT::Z_int(double T, const kappa::Molecule &molecule, int num_electron_levels) {
	
    if (num_electron_levels == -1) {
      return Z_int(T, molecule.electron_energy, 
   		      molecule.statistical_weight, 
		      molecule.num_electron_levels,
		      molecule.vibr_energy, 
                      molecule.num_vibr_levels, 
		      molecule.rot_energy, 
		      molecule.num_rot_levels, 
		      molecule.rot_symmetry);
     } else {
       return Z_int(T, molecule.electron_energy, 
   		       molecule.statistical_weight, 
		       num_electron_levels,
		       molecule.vibr_energy, 
		       molecule.num_vibr_levels, 
		       molecule.rot_energy, 
		       molecule.num_rot_levels, 
		       molecule.rot_symmetry);
     }
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationOneT::avg_energy(double T, const arma::vec &electron_energy, const arma::Col<unsigned long> &statistical_weight, int n_electron_levels) {

    double res=0.0;
    for (int e=0; e<n_electron_levels; e++) {
      res += statistical_weight[e] * electron_energy[e] * exp(-electron_energy[e] / (K_CONST_K * T));
    }
    return res / Z_int(T, electron_energy, statistical_weight, n_electron_levels);
  }   

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Averaged internal energy for an atom
  double kappa::ApproximationOneT::avg_energy(double T, const kappa::Atom &atom, double Delta_E) {

    if (Delta_E == -1) {
      return 0.0;
    } else {
      return avg_energy(T, atom.electron_energy, atom.statistical_weight, p_max_electron_level(atom.electron_energy, atom.ionization_potential, Delta_E));
    }
  }
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationOneT::avg_energy(	double T, 
						const arma::vec &electron_energy, 
						const arma::Col<unsigned long> &statistical_weight, 
						int n_electron_levels,
                       				const std::vector<arma::vec> &vibr_energy, 
						const std::vector<int> &num_vibr_levels,
                       				const std::vector<std::vector<arma::vec>> &rot_energy, 
						const std::vector<std::vector<int>> &num_rot_levels, 
						int rot_symmetry) {

    double res=0.0;
    for (int e=0; e<n_electron_levels; e++) {
      for (int i=0; i<num_vibr_levels[e]; i++) {
        res += statistical_weight[e] * arma::dot((2*arma::linspace<arma::vec>(0, num_rot_levels[e][i]-1, num_rot_levels[e][i])+1) 
                                     % (electron_energy[e] + vibr_energy[e][i] + rot_energy[e][i]),
            		    	       arma::exp(-(electron_energy[e] + vibr_energy[e][i] + rot_energy[e][i]) / (K_CONST_K * T)));
      }
    }
    return res / (rot_symmetry * Z_int(	T, 
 					electron_energy, 
					statistical_weight, 
					n_electron_levels, 
					vibr_energy, 
					num_vibr_levels,
                       		   	rot_energy, 
					num_rot_levels, 
					rot_symmetry));
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Averaged internal energy for a molecule
  double kappa::ApproximationOneT::avg_energy(double T, const kappa::Molecule &molecule, int num_electron_levels) {
	
    if (num_electron_levels == -1) {
      return avg_energy(T, molecule.electron_energy, 
  		           molecule.statistical_weight, 
			   molecule.num_electron_levels,
			   molecule.vibr_energy, 
			   molecule.num_vibr_levels, 
		           molecule.rot_energy, 
			   molecule.num_rot_levels, 
			   molecule.rot_symmetry);
    } else {
      return avg_energy(T, molecule.electron_energy, 
			   molecule.statistical_weight, 
			   num_electron_levels,
			   molecule.vibr_energy, 
			   molecule.num_vibr_levels, 
			   molecule.rot_energy, 
			   molecule.num_rot_levels, 
			   molecule.rot_symmetry);
    }
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationOneT::avg_energy_sq(double T, const arma::vec &electron_energy, const arma::Col<unsigned long> &statistical_weight, int n_electron_levels) {
    double res=0.0;
    for (int e=0; e<n_electron_levels; e++) {
      res += statistical_weight[e] * electron_energy[e] * electron_energy[e] * exp(-electron_energy[e] / (K_CONST_K * T));
    }
    return res / Z_int(T, electron_energy, statistical_weight, n_electron_levels);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationOneT::avg_energy_sq(	double T, 
							const arma::vec &electron_energy, 
							const arma::Col<unsigned long> &statistical_weight, 
							int n_electron_levels,
                          				const std::vector<arma::vec> &vibr_energy, 
							const std::vector<int> &num_vibr_levels,
                          				const std::vector<std::vector<arma::vec>> &rot_energy, 
							const std::vector<std::vector<int>> &num_rot_levels, 
							int rot_symmetry) {
    double res=0.0;
    for (int e=0; e<n_electron_levels; e++) {
      for (int i=0; i<num_vibr_levels[e]; i++) {
        res += statistical_weight[e] * arma::dot((2*arma::linspace<arma::vec>(0, num_rot_levels[e][i]-1, num_rot_levels[e][i])+1) 
                 % (electron_energy[e] + vibr_energy[e][i] + rot_energy[e][i])
		 % (electron_energy[e] + vibr_energy[e][i] + rot_energy[e][i]),
		 arma::exp(-(electron_energy[e] + vibr_energy[e][i] + rot_energy[e][i]) / (K_CONST_K * T)));
      }
    }
   return res / (rot_symmetry * Z_int(T, electron_energy, statistical_weight, n_electron_levels, vibr_energy, num_vibr_levels, rot_energy, num_rot_levels, rot_symmetry));
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Averaged squared internal energy of an atom
  double kappa::ApproximationOneT::avg_energy_sq(double T, const kappa::Atom &atom, double Delta_E) {
	
    if (Delta_E == -1) {
      return 0.0;
    } else {
      return avg_energy_sq(T, atom.electron_energy, atom.statistical_weight, p_max_electron_level(atom.electron_energy, atom.ionization_potential, Delta_E));
    }
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Averaged squared internal energy of an molecule
  double kappa::ApproximationOneT::avg_energy_sq(double T, const kappa::Molecule &molecule, int num_electron_levels) {

    if (num_electron_levels == -1) {
      return avg_energy_sq(T, molecule.electron_energy, 
	 		      molecule.statistical_weight, 
			      molecule.num_electron_levels,
	   	              molecule.vibr_energy, 
			      molecule.num_vibr_levels, 
			      molecule.rot_energy, 
			      molecule.num_rot_levels, 
			      molecule.rot_symmetry);
    } else {
      return avg_energy_sq(T, molecule.electron_energy, 
     			      molecule.statistical_weight, 
			      num_electron_levels,
			      molecule.vibr_energy, 
			      molecule.num_vibr_levels, 
			      molecule.rot_energy, 
			      molecule.num_rot_levels, 
			      molecule.rot_symmetry);
    }
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationOneT::c_int(double T, double mass, const arma::vec &electron_energy, const arma::Col<unsigned long> &statistical_weight, int n_electron_levels) {

    double ae = avg_energy(T, electron_energy, statistical_weight, n_electron_levels);
    return (avg_energy_sq(T, electron_energy, statistical_weight, n_electron_levels) - ae * ae) / (K_CONST_K * T * T * mass);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationOneT::c_int(	double T, 
						double mass, 
						const arma::vec &electron_energy, 
						const arma::Col<unsigned long> &statistical_weight, 
						int n_electron_levels,
                          			const std::vector<arma::vec> &vibr_energy, 
						const std::vector<int> &num_vibr_levels,
                          			const std::vector<std::vector<arma::vec>> &rot_energy, 
						const std::vector<std::vector<int>> &num_rot_levels, 
						int rot_symmetry) {

    double ae = avg_energy(T, electron_energy, statistical_weight, n_electron_levels, vibr_energy, num_vibr_levels, rot_energy, num_rot_levels, rot_symmetry);
    return (avg_energy_sq(T, electron_energy, statistical_weight, n_electron_levels, vibr_energy, num_vibr_levels, 
            rot_energy, num_rot_levels, rot_symmetry) - ae * ae) / (K_CONST_K * T * T * mass);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  // Specific heat capacity of internal degrees of freedom of an atom
  double kappa::ApproximationOneT::c_int(double T, const kappa::Atom &atom, double Delta_E) {
	
    if (Delta_E == -1) {
      return 0.0;
    } else {
      return c_int(T, atom.mass, atom.electron_energy, atom.statistical_weight, p_max_electron_level(atom.electron_energy, atom.ionization_potential, Delta_E));
    }
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Specific heat capacity of internal degrees of freedom of an molecule
  double kappa::ApproximationOneT::c_int(double T, const kappa::Molecule &molecule, int num_electron_levels) {
	
    if (num_electron_levels == -1) {
      return c_int(T, molecule.mass, 
                      molecule.electron_energy, 
                      molecule.statistical_weight, 
                      molecule.num_electron_levels,
   	              molecule.vibr_energy, 
                      molecule.num_vibr_levels, 
		      molecule.rot_energy, 
		      molecule.num_rot_levels, 
		      molecule.rot_symmetry);
    } else {
      return c_int(T, molecule.mass, 
    		      molecule.electron_energy, 
		      molecule.statistical_weight, 
		      num_electron_levels,
		      molecule.vibr_energy, 
		      molecule.num_vibr_levels, 
		      molecule.rot_energy, 
		      molecule.num_rot_levels, 
		      molecule.rot_symmetry);
    }
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
