/*!
    \file models.h
    \brief String values of models
 */

#ifndef kappa_models_h
#define kappa_models_h

namespace kappa {

enum class models_prob_vv {	model_prob_vv_fho};
enum class models_prob_vt {	model_prob_vt_fho};

enum class models_prob_diss {	model_prob_diss_thresh_cmass_vibr, 
				model_prob_diss_thresh_vibr, 
				model_prob_diss_thresh_cmass, 
				model_prob_diss_thresh};

enum class models_cs_elastic {	model_cs_el_rs,	    // The model of hard spheres (RS) in combination with the FHO model
				model_cs_el_vss};                   // The model of spheres of variable diameter (VSS) in combination with model FHO

enum class models_cs_vv {	model_cs_vv_rs_fho, 
				model_cs_vv_vss_fho};

enum class models_cs_vt {	model_cs_vt_rs_fho,		  
				model_cs_vt_vss_fho}; 		 

enum class models_cs_diss {	model_cs_diss_rs_thresh_cmass_vibr, // RS model in combination with a threshold model that takes into account vibrational and translational energy along the center of mass line
				model_cs_diss_rs_thresh_vibr,       // RS model in combination with a threshold model that takes into account vibrational and total translational energy
				model_cs_diss_rs_thresh_cmass,      // RS model in combination with a threshold model that takes into account translational energy along the center of mass line
				model_cs_diss_rs_thresh, 	    // RS model in combination with a threshold model that takes into account total translational energy
				model_cs_diss_vss_thresh_cmass_vibr,// VSS model in combination with a threshold model that takes into account vibrational and translational energy along the center of mass line 
				model_cs_diss_vss_thresh_vibr, 	    // VSS model in combination with a threshold model that takes into account vibrational and total translational energy
				model_cs_diss_vss_thresh_cmass,     // VSS model in combination with a threshold model that takes into account translational energy along the center of mass line
				model_cs_diss_vss_thresh,           // VSS model in combination with a threshold model that takes into account translational energy
				model_cs_diss_ilt};	  	    // a model based on the Laplace transform from the results of quasi-classical trajectory calculations in the Phys4Entry database

enum class models_k_vv {	model_k_vv_rs_fho, 
				model_k_vv_vss_fho, 
				model_k_vv_ssh, 
				model_k_vv_billing};

enum class models_k_vt {	model_k_vt_rs_fho, 
				model_k_vt_vss_fho, 
				model_k_vt_ssh, 
				model_k_vt_phys4entry,
				model_k_vt_billing};

enum class models_k_exch {	model_k_exch_arrh_scanlon, 
				model_k_exch_arrh_park, 
				model_k_exch_warnatz, 
				model_k_exch_rf, // Rusanov-Friedman model (doesn't take into account the electr. excitation of the reagent molecule and the vibr. excitation of the product molecule)
				model_k_exch_polak,
			  	model_k_exch_maliat_D6k_arrh_scanlon,   // Modified Aliat model with data from Scanlon et al. with \f$ U=\frac{D}{6k} \f$
				model_k_exch_maliat_3T_arrh_scanlon, 
				model_k_exch_maliat_infty_arrh_scanlon,
			  	model_k_exch_maliat_D6k_arrh_park, 
				model_k_exch_maliat_3T_arrh_park, 
				model_k_exch_maliat_infty_arrh_park};

enum class models_k_diss {	model_k_diss_rs_thresh_cmass_vibr, 
				model_k_diss_rs_thresh_vibr, 
				model_k_diss_rs_thresh_cmass, 
				model_k_diss_rs_thresh,
				model_k_diss_vss_thresh_cmass_vibr, 
				model_k_diss_vss_thresh_vibr, 
				model_k_diss_vss_thresh_cmass, 
				model_k_diss_vss_thresh,
				model_k_diss_arrh_scanlon, 
				model_k_diss_arrh_park,
				model_k_diss_tm_D6k_arrh_scanlon,  	// The Treanor-Marrone model with the data from Scanlon et al. with \f$ U=\frac{D}{6k} \f$
				model_k_diss_tm_3T_arrh_scanlon,	// The Treanor-Marrone model with the data from Scanlon et al. with \f$ U=3T \f$
				model_k_diss_tm_infty_arrh_scanlon,	// The Treanor-Marrone model with the data from Scanlon et al. with \f$ U=\infty \f$
				model_k_diss_tm_D6k_arrh_park,  
				model_k_diss_tm_3T_arrh_park, 
				model_k_diss_tm_infty_arrh_park,
				model_k_diss_phys4entry, 
				model_k_diss_ilt};

enum class models_omega {	model_omega_rs, 
				model_omega_vss, 
				model_omega_bornmayer, 
				model_omega_lennardjones, 
				model_omega_esa};
}
#endif /* models_h */
