/*
   \file approximation_multit.cpp
 */

#include "approximation_multit.hpp"

kappa::ApproximationMultiT::ApproximationMultiT() : kappa::Approximation() { } 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // non-private rotational partition function
  double kappa::ApproximationMultiT::Z_rot(double T, const kappa::Molecule &molecule) {

    // private rotational partition function
    return p_Z_rot(T, molecule.rot_energy[0][0], molecule.num_rot_levels[0][0], molecule.rot_symmetry);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::avg_vibr_energy(double T, double T1, const kappa::Molecule & molecule) {

    int i_star = max_i(T, T1, molecule);
    if (i_star >= molecule.num_vibr_levels[0]-1) {
      return p_avg_vibr_energy(T, T1, molecule.vibr_energy[0] - molecule.vibr_energy[0][0], molecule.num_vibr_levels[0]); 
    } else {
      return p_avg_vibr_energy(T, T1, molecule.vibr_energy[0].subvec(0, i_star) - molecule.vibr_energy[0][0], i_star + 1); }
    }
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::avg_vibr_energy_sq(double T, double T1, const kappa::Molecule & molecule) {

    int i_star = max_i(T, T1, molecule);
    if (i_star >= molecule.num_vibr_levels[0] - 1) {
      return p_avg_vibr_energy_sq(T, T1, molecule.vibr_energy[0] - molecule.vibr_energy[0][0], molecule.num_vibr_levels[0]);
    } else {
      return p_avg_vibr_energy_sq(T, T1, molecule.vibr_energy[0].subvec(0, i_star) - molecule.vibr_energy[0][0], i_star + 1);
    }
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::avg_vibr_i(double T, double T1, const kappa::Molecule & molecule) {

    int i_star = max_i(T, T1, molecule);
    if (i_star >= molecule.num_vibr_levels[0] - 1) {
      return p_avg_vibr_i(T, T1, molecule.vibr_energy[0] - molecule.vibr_energy[0][0], molecule.num_vibr_levels[0]);
    } else {
      return p_avg_vibr_i(T, T1, molecule.vibr_energy[0].subvec(0, i_star) - molecule.vibr_energy[0][0], i_star + 1);
    }
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::avg_vibr_i_sq(double T, double T1, const kappa::Molecule & molecule) {

    int i_star = max_i(T, T1, molecule);
    if (i_star >= molecule.num_vibr_levels[0] - 1) {
     return p_avg_vibr_i_sq(T, T1, molecule.vibr_energy[0] - molecule.vibr_energy[0][0], molecule.num_vibr_levels[0]);
    } else {
      return p_avg_vibr_i_sq(T, T1, molecule.vibr_energy[0].subvec(0, i_star) - molecule.vibr_energy[0][0], i_star + 1);
    }
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::avg_vibr_i_energy(double T, double T1, const kappa::Molecule & molecule) {

    int i_star = max_i(T, T1, molecule);
    if (i_star >= molecule.num_vibr_levels[0] - 1) {
      return p_avg_vibr_i_energy(T, T1, molecule.vibr_energy[0] - molecule.vibr_energy[0][0], molecule.num_vibr_levels[0]);
    } else {
      return p_avg_vibr_i_energy(T, T1, molecule.vibr_energy[0].subvec(0, i_star) - molecule.vibr_energy[0][0], i_star + 1);
    }
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::c_rot(double T, double T1, const kappa::Molecule & molecule) {

   return p_c_rot(T, molecule.mass, molecule.rot_energy[0][0], molecule.num_rot_levels[0][0], molecule.rot_symmetry);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::c_vibr_T(double T, double T1, const kappa::Molecule & molecule) {

    int i_star = max_i(T, T1, molecule);
    if (i_star >= molecule.num_vibr_levels[0] - 1) {
      return p_c_vibr_T(T, T1, molecule.mass, molecule.vibr_energy[0] - molecule.vibr_energy[0][0], molecule.num_vibr_levels[0]);
    } else {
      return p_c_vibr_T(T, T1, molecule.mass, molecule.vibr_energy[0].subvec(0, i_star) - molecule.vibr_energy[0][0], i_star + 1);
    }
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::c_vibr_T1(double T, double T1, const kappa::Molecule & molecule) {

    int i_star = max_i(T, T1, molecule);
    if (i_star >= molecule.num_vibr_levels[0] - 1) {
      return p_c_vibr_T1(T, T1, molecule.mass, molecule.vibr_energy[0] - molecule.vibr_energy[0][0], molecule.num_vibr_levels[0]);	
    } else {
      return p_c_vibr_T1(T, T1, molecule.mass, molecule.vibr_energy[0].subvec(0, i_star) - molecule.vibr_energy[0][0], i_star + 1);	
    }
  }  

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::c_W_T(double T, double T1, const kappa::Molecule & molecule) {

    int i_star = max_i(T, T1, molecule);
    if (i_star >= molecule.num_vibr_levels[0] - 1) {
      return p_c_W_T(T, T1, molecule.mass, molecule.vibr_energy[0] - molecule.vibr_energy[0][0], molecule.num_vibr_levels[0]);
    } else {
      return p_c_W_T(T, T1, molecule.mass, molecule.vibr_energy[0].subvec(0, i_star) - molecule.vibr_energy[0][0], i_star + 1);
    }
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::c_W_T1(double T, double T1, const kappa::Molecule & molecule) {

    int i_star = max_i(T, T1, molecule);
    if (i_star >= molecule.num_vibr_levels[0] - 1) {
      return p_c_W_T1(T, T1, molecule.mass, molecule.vibr_energy[0] - molecule.vibr_energy[0][0], molecule.num_vibr_levels[0]);
    } else {
      return p_c_W_T1(T, T1, molecule.mass, molecule.vibr_energy[0].subvec(0, i_star) - molecule.vibr_energy[0][0], i_star + 1);
    }
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::k_diss(double T, double T1, kappa::Molecule const & molecule, kappa::Interaction const & interaction, kappa::models_k_diss model) {

    double tmp = 0.0;
    int max_level = max_i(T, T1, molecule);
    for (int i = 0; i<max_level; i++) {
      tmp += k_diss(T, molecule, interaction, i, model) * exp(((-(molecule.vibr_energy[0][i] - molecule.vibr_energy[0][0]) 
     				                               + (molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0]) * i) / (K_CONST_K * T))
					                     + (-(molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0]) * i / (K_CONST_K * T1)));
    }
    return tmp / Z_vibr(T, T1, molecule);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::k_diss_i(double T, double T1, kappa::Molecule const & molecule, kappa::Interaction const & interaction, kappa::models_k_diss model) {

    double tmp = 0.0;
    int max_level = max_i(T, T1, molecule);
    for (int i = 0; i<max_level; i++) {
      tmp += i*k_diss(T, molecule, interaction, i, model) * exp(((-(molecule.vibr_energy[0][i] - molecule.vibr_energy[0][0])
   			                                         + (molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0]) * i) / (K_CONST_K * T))
							       + (-(molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0]) * i / (K_CONST_K * T1)));
    }
    return tmp / Z_vibr(T, T1, molecule);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::k_diss_ve(double T, double T1, kappa::Molecule const & molecule, kappa::Interaction const & interaction, kappa::models_k_diss model) {

    double tmp = 0.0;
    int max_level = max_i(T, T1, molecule);
    for (int i = 0; i<max_level; i++) {
      tmp += (molecule.vibr_energy[0][i] - molecule.vibr_energy[0][0])*k_diss(T, molecule, interaction, i, model) 
 				 * exp(((-(molecule.vibr_energy[0][i] - molecule.vibr_energy[0][0]) + (molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0])*i) / (K_CONST_K * T))
				      + (-(molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0]) * i / (K_CONST_K * T1)));
    }
    return tmp / Z_vibr(T, T1, molecule);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::k_VT_i(double T, double T1, kappa::Molecule const & molecule, kappa::Interaction const & interaction, int max_delta_i, kappa::models_k_vt model) {

    double tmp = 0.0;
    int max_level = max_i(T, T1, molecule);
    for (int i=0; i <= max_level; i++) {
      for (int delta_i = std::max(-i, -max_delta_i); delta_i <= std::min(max_level - i, max_delta_i); delta_i++) {
 	if (delta_i != 0) {
 	  tmp += i*( exp(((-(molecule.vibr_energy[0][delta_i + i] - molecule.vibr_energy[0][0])
			  + (molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0]) * delta_i) / (K_CONST_K * T))
			+ (-(molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0]) * delta_i / (K_CONST_K * T1)))
			* k_VT(T, molecule, interaction, i + delta_i, -delta_i, model) 
			- exp(((-(molecule.vibr_energy[0][i] - molecule.vibr_energy[0][0])
			+ (molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0]) * i) / (K_CONST_K * T))
		      + (-(molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0]) * i / (K_CONST_K * T1)))
		      * k_VT(T, molecule, interaction, i, delta_i, model) );
	}
      }
    }
    return tmp;
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::k_VT_ve(double T, double T1, kappa::Molecule const & molecule, kappa::Interaction const & interaction, int max_delta_i, kappa::models_k_vt model) {

    double tmp = 0.0;
    int max_level = max_i(T, T1, molecule);
    for (int i = 0; i <= max_level; i++) {
      for (int delta_i = std::max(-i, -max_delta_i); delta_i <= std::min( max_level - i, max_delta_i ); delta_i++) {
	if (delta_i != 0) {
  	 tmp += (molecule.vibr_energy[0][i] - molecule.vibr_energy[0][0])
                                   *( exp(((-(molecule.vibr_energy[0][delta_i + i] - molecule.vibr_energy[0][0])
				           + (molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0]) * delta_i) / (K_CONST_K * T))
				 	 + (-(molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0]) * delta_i / (K_CONST_K * T1)))
					 * k_VT(T, molecule, interaction, i + delta_i, -delta_i, model)
					- exp(((-(molecule.vibr_energy[0][i] - molecule.vibr_energy[0][0])
					+ (molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0])*i) / (K_CONST_K * T))
				      + (-(molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0]) * i / (K_CONST_K * T1)))
				      * k_VT(T, molecule, interaction, i, delta_i, model) );
	}
      }
    }
    return tmp;;
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::Z_vibr(double T, double T1, const kappa::Molecule &molecule) {

    if ((molecule.anharmonic_spectrum == false) || (T >= T1)) {
      return p_Z_vibr(T, T1, molecule.vibr_energy[0] - molecule.vibr_energy[0][0], molecule.num_vibr_levels[0]);
    } else {
      int i_star = max_i(T, T1, molecule);
      if (i_star < 1) {
 	std::string error_string = "No Treanor distribution possible for such values of T, T1";
	throw kappa::ModelParameterException(error_string.c_str());
      } else if (i_star >= molecule.num_vibr_levels[0] - 1) {
	return p_Z_vibr(T, T1, molecule.vibr_energy[0] - molecule.vibr_energy[0][0], molecule.num_vibr_levels[0]);
      } else {
	return p_Z_vibr(T, T1, molecule.vibr_energy[0].subvec(0, i_star) - molecule.vibr_energy[0][0], i_star + 1);
      }
    }
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  // Returns the maximum vibrational level possible for the molecule
  int kappa::ApproximationMultiT::max_i(double T, double T1, const kappa::Molecule &molecule) {
    
    if ((molecule.anharmonic_spectrum == false) || (T >= T1)) {
      return molecule.num_vibr_levels[0]; 
    } else {
      int i_star = p_max_i(T, T1, molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0], K_CONST_C * molecule.vibr_frequency[0], molecule.vibr_we_xe[0] / molecule.vibr_frequency[0]);
      if (i_star < 1) {
        std::string error_string = "No Treanor distribution possible for such values of T, T1";
        throw kappa::ModelParameterException(error_string.c_str());	
      } else if (i_star >= molecule.num_vibr_levels[0] - 1) {
        return molecule.num_vibr_levels[0] - 1;	
      } else {
        return p_max_i(T, T1, molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0], K_CONST_C * molecule.vibr_frequency[0], molecule.vibr_we_xe[0] / molecule.vibr_frequency[0]);
      }
    }
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::p_Z_vibr(	double T, double T1, const arma::vec &vibr_energy, int num_vibr_levels) {

    // eqn. 3.13
    return arma::sum(arma::exp(
         ( (-vibr_energy + (vibr_energy[1])*arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels))/ (K_CONST_K * T) ) +
 	   (               -vibr_energy[1]* arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels) / (K_CONST_K * T1)) ));
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::p_avg_vibr_energy(double T, double T1, const arma::vec & vibr_energy, int num_vibr_levels) {
		
    return (arma::dot(vibr_energy, arma::exp(((-vibr_energy
 	           + (vibr_energy[1])* arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)) / (K_CONST_K * T))
		   + (-vibr_energy[1] * arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)) / (K_CONST_K * T1))))
		   / p_Z_vibr(T, T1, vibr_energy, num_vibr_levels);
  }
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::p_avg_vibr_energy_sq(double T, double T1, const arma::vec & vibr_energy, int num_vibr_levels) {

    return (arma::dot(vibr_energy%vibr_energy, arma::exp(((-vibr_energy
 		    + vibr_energy[1] * arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)) / (K_CONST_K * T))
  		   + (-vibr_energy[1] * arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)) / (K_CONST_K * T1))))
		   / p_Z_vibr(T, T1, vibr_energy, num_vibr_levels);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::p_avg_vibr_i(double T, double T1, const arma::vec & vibr_energy, int num_vibr_levels) {

    return (arma::dot(arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels), arma::exp(((-vibr_energy
	  		  + (vibr_energy[1])* arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)) / (K_CONST_K * T))
			 + (-vibr_energy[1] * arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)) / (K_CONST_K * T1))))
			/ p_Z_vibr(T, T1, vibr_energy, num_vibr_levels);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::p_avg_vibr_i_sq(double T, double T1, const arma::vec & vibr_energy, int num_vibr_levels) {

    return (arma::dot(arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)%
                      arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels),
                      arma::exp(((-vibr_energy + (vibr_energy[1]) * 
                      arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)) / (K_CONST_K * T)) + 
                      (-vibr_energy[1] * arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)) / (K_CONST_K * T1))))
	              / p_Z_vibr(T, T1, vibr_energy, num_vibr_levels);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::p_avg_vibr_i_energy(double T, double T1, const arma::vec & vibr_energy, int num_vibr_levels) {

    return (arma::dot(arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)% vibr_energy, arma::exp(((-vibr_energy
			+ (vibr_energy[1])*arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)) / (K_CONST_K * T))
			+ (-vibr_energy[1] * arma::linspace<arma::vec>(0, num_vibr_levels - 1, num_vibr_levels)) / (K_CONST_K * T1))))
			/ p_Z_vibr(T, T1, vibr_energy, num_vibr_levels);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::p_c_vibr_T(double T, double T1, double mass, const arma::vec & vibr_energy, int num_vibr_levels) {

    double tmp = p_avg_vibr_energy(T, T1, vibr_energy, num_vibr_levels);
    return (p_avg_vibr_energy_sq(T, T1, vibr_energy, num_vibr_levels) - (tmp * tmp)
 			              + vibr_energy[1] * p_avg_vibr_i(T, T1, vibr_energy, num_vibr_levels) 
                                      * p_avg_vibr_energy(T, T1, vibr_energy, num_vibr_levels) 
			    	      - vibr_energy[1] * p_avg_vibr_i_energy(T, T1, vibr_energy, num_vibr_levels) )/(K_CONST_K * K_CONST_K * T * T);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::p_c_vibr_T1(double T, double T1, double mass, const arma::vec & vibr_energy, int num_vibr_levels) {

   return ( vibr_energy[1]*(p_avg_vibr_i_energy(T, T1, vibr_energy, num_vibr_levels)
	                  - p_avg_vibr_i(T, T1, vibr_energy, num_vibr_levels) * 
 			    p_avg_vibr_energy(T, T1, vibr_energy, num_vibr_levels)) )/(K_CONST_K * K_CONST_K * T1 * T1);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::p_c_W_T(double T, double T1, double mass, const arma::vec & vibr_energy, int num_vibr_levels) {

    double tmp = p_avg_vibr_i(T, T1, vibr_energy, num_vibr_levels);
    return (vibr_energy[1] * (p_avg_vibr_i_energy(T, T1, vibr_energy, num_vibr_levels) 
 	  - vibr_energy[1]* p_avg_vibr_i_sq(T, T1, vibr_energy, num_vibr_levels)
	  - p_avg_vibr_i(T, T1, vibr_energy, num_vibr_levels)*p_avg_vibr_energy(T, T1, vibr_energy, num_vibr_levels) 
	  + vibr_energy[1]*tmp*tmp)) / (K_CONST_K*K_CONST_K*T*T);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::ApproximationMultiT::p_c_W_T1(double T, double T1, double mass, const arma::vec & vibr_energy, int num_vibr_levels) {

    double tmp = p_avg_vibr_i(T, T1, vibr_energy, num_vibr_levels);
    return (((vibr_energy[1]) * (vibr_energy[1])) * (p_avg_vibr_i_sq(T, T1, vibr_energy, num_vibr_levels) - tmp*tmp))/(K_CONST_K * K_CONST_K * T1 * T1);
  } 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int kappa::ApproximationMultiT::p_max_i(double T, double T1, double vibr_energy1, double vibr_frequency, double alpha) {

    int p_i = -1;
    p_i = (vibr_energy1*T) / (2 * alpha*vibr_frequency*T1*K_CONST_H) + 0.5;
    return p_i;
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
