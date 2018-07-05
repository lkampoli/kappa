/*!
    \file approximation_multit.hpp
 */

#ifndef kappa_approximation_multit_hpp
#define kappa_approximation_multit_hpp

#include "approximation.hpp"
#include <algorithm>

namespace kappa {

class ApproximationMultiT : public Approximation {

using Approximation::k_diss;

 protected:

  double p_Z_vibr(double T, double T1, const arma::vec & vibr_energy, int num_vibr_levels);
  double p_avg_vibr_energy(double T, double T1, const arma::vec & vibr_energy, int num_vibr_levels);
  double p_avg_vibr_energy_sq(double T, double T1, const arma::vec & vibr_energy, int num_vibr_levels);
  double p_avg_vibr_i(double T, double T1, const arma::vec & vibr_energy, int num_vibr_levels);
  double p_avg_vibr_i_sq(double T, double T1, const arma::vec & vibr_energy, int num_vibr_levels);
  double p_avg_vibr_i_energy(double T, double T1, const arma::vec & vibr_energy, int num_vibr_levels);
  double p_c_vibr_T(double T, double T1, double mass,  const arma::vec & vibr_energy, int num_vibr_levels);
  double p_c_vibr_T1(double T, double T1, double mass, const arma::vec & vibr_energy, int num_vibr_levels);
  double p_c_W_T(double T, double T1, double mass, const arma::vec & vibr_energy, int num_vibr_levels);
  double p_c_W_T1(double T, double T1, double mass, const arma::vec & vibr_energy, int num_vibr_levels);
  int p_max_i(double T, double T1, double vibr_energy1, double vibr_frequency, double alpha);

 public:

  ApproximationMultiT(); 

  double Z_rot(double T, const kappa::Molecule &molecule);
  double Z_vibr(double T, double T1, const kappa::Molecule &molecule);
  int max_i(double T, double T1, const kappa::Molecule & molecule);
  double avg_vibr_energy(double T, double T1, const kappa::Molecule &molecule);
  double avg_vibr_energy_sq(double T, double T1, const kappa::Molecule &molecule);
  double avg_vibr_i(double T, double T1, const kappa::Molecule &molecule); 
  double avg_vibr_i_sq(double T, double T1, const kappa::Molecule &molecule);
  double avg_vibr_i_energy(double T, double T1, const kappa::Molecule &molecule);
  double c_rot(double T, double T1, const kappa::Molecule &molecule);
  double c_vibr_T(double T, double T1, const kappa::Molecule &molecule);
  double c_vibr_T1(double T, double T1, const kappa::Molecule &molecule);
  double c_W_T(double T, double T1, const kappa::Molecule &molecule);
  double c_W_T1(double T, double T1, const kappa::Molecule &molecule);
  
  double k_diss(double T, 
		double T1, 
		kappa::Molecule const &molecule, 
		kappa::Interaction const &interaction,
		kappa::models_k_diss model=kappa::models_k_diss::model_k_diss_vss_thresh_cmass_vibr);
                double k_diss_i(double T, 
		double T1, 
		kappa::Molecule const & molecule, 
		kappa::Interaction const & interaction, 
		kappa::models_k_diss model=kappa::models_k_diss::model_k_diss_vss_thresh_cmass_vibr);

  double k_diss_ve(double T, 
  		   double T1, 
  		   kappa::Molecule const & molecule, 
  		   kappa::Interaction const & interaction, 
  		   kappa::models_k_diss model=kappa::models_k_diss::model_k_diss_vss_thresh_cmass_vibr);
  
  double k_VT_i(double T, 
  		double T1, 
  		kappa::Molecule const &molecule, 
  		kappa::Interaction const &interaction, 
  		int max_delta_i,
  		kappa::models_k_vt model = kappa::models_k_vt::model_k_vt_vss_fho);
  
  double k_VT_ve(double T, 
  		 double T1, 
  		 kappa::Molecule const &molecule, 
  		 kappa::Interaction const &interaction, 
  		 int max_delta_i,
  		 kappa::models_k_vt model = kappa::models_k_vt::model_k_vt_vss_fho);

};
};
#endif /* kappa_approximation_sts_hpp */
