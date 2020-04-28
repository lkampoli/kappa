/*!
    \file approximation.cpp 
 */

#include "approximation.hpp"

kappa::Approximation::Approximation() {}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Gets the maximum electronic level before ionization
  int kappa::Approximation::p_max_electron_level(const arma::vec &electron_energy, double ionization_potential, double Delta_E) {

    int i=0;
    while (electron_energy[i+1] <= ionization_potential - Delta_E) {
      i++;
    } 
    return i;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Computes the rotational partition function
  double kappa::Approximation::p_Z_rot(double T, const arma::vec &rot_energy, int num_rot_levels, int rot_symmetry) {
    return arma::dot(2 * arma::linspace<arma::vec>(0, num_rot_levels - 1, num_rot_levels) + 1, arma::exp(-rot_energy / (K_CONST_K * T))) / rot_symmetry;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Computes the equilibrium vibrational partition function
  double kappa::Approximation::p_Z_vibr_eq(double T, const arma::vec &vibr_energy) {
    return arma::sum(arma::exp(-vibr_energy / (K_CONST_K * T)));
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Computes the electronic partition function
  double kappa::Approximation::p_Z_electron(double T, const arma::vec &electron_energy, const arma::Col<unsigned long> &statistical_weight, int n_electron_levels) {

    double res=0.0;
    for (int e=0; e<n_electron_levels; e++) {
      res += statistical_weight[e] * exp(-electron_energy[e] / (K_CONST_K * T));
    }
    return res;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::Z_diss(	double T, 
					double U, 
					const arma::vec &electron_energy, 
					const arma::Col<unsigned long> &statistical_weight, 
					int n_electron_levels,
                                    	const std::vector<arma::vec> &vibr_energy, 
					int i, 
					int e) {
    double res = 0.;
    for (int j=0; j<n_electron_levels; j++) {
      res += statistical_weight[j] * exp(electron_energy[j] / U) * p_Z_vibr_eq(-U, vibr_energy[j]) / p_Z_vibr_eq(T, vibr_energy[j]);
    }
    return p_Z_electron(T, electron_energy, statistical_weight, n_electron_levels) * exp((vibr_energy[e][i] + electron_energy[e]) * (1. / T + 1. / U) / K_CONST_K) / res;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::Z_diss(	double T, 
					const arma::vec &electron_energy, 
					const arma::Col<unsigned long> &statistical_weight, 
					int n_electron_levels,
                                    	const std::vector<arma::vec> &vibr_energy, 
					const std::vector<int> &num_vibr_levels, 
					int i, 
					int e) {
    double res = 0.;
    for (int j=0; j<n_electron_levels; j++) {
      res += statistical_weight[j] * num_vibr_levels[j] / p_Z_vibr_eq(T, vibr_energy[j]);
    }
    return p_Z_electron(T, electron_energy, statistical_weight, n_electron_levels) * exp((vibr_energy[e][i] + electron_energy[e]) * (1. / T) / K_CONST_K) / res;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::Z_diss(double T, double U, const arma::vec &vibr_energy, int i) {
    return p_Z_vibr_eq(T, vibr_energy) * exp(vibr_energy[i] * (1. / T + 1. / U) / K_CONST_K) / p_Z_vibr_eq(-U, vibr_energy);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::Z_diss(double T, const arma::vec &vibr_energy, int num_vibr_levels, int i) {
    return p_Z_vibr_eq(T, vibr_energy) * exp(vibr_energy[i] / (K_CONST_K * T)) / num_vibr_levels;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::C_aliat(	double T, 
					const arma::vec &electron_energy, 
					const arma::Col<unsigned long> &statistical_weight, 
                                        int n_electron_levels, 
					const std::vector<arma::vec> &vibr_energy, 
					const std::vector<int> &num_vibr_levels,
                                       	double vibr_energy_product, 
					double activation_energy, 
					int i, 
					int e, 
					double U) { // U != infty
    double tmp=0, res=0;
    double kT = K_CONST_K * T, kU = K_CONST_K * U;
    double i_e;

    for (int ee=0; ee<n_electron_levels; ee++) {
      tmp = 0;
      for (int k = 0; k < num_vibr_levels[ee]; k++) {
        i_e = vibr_energy[ee][k] + electron_energy[ee];
        if (i_e > activation_energy + vibr_energy_product) {
          tmp += statistical_weight[ee] * exp((activation_energy + vibr_energy_product - i_e) / kT);
        } else {
          tmp += statistical_weight[ee] * exp(-(activation_energy + vibr_energy_product - i_e) / kU);
        }
      }
      res += tmp / p_Z_vibr_eq(T, vibr_energy[ee]);
    }
    res = 1./res;
    i_e = vibr_energy[e][i] + electron_energy[e];
    if (i_e > activation_energy + vibr_energy_product) {
      res *= exp((activation_energy + vibr_energy_product) / kT);
    } else {
      res *= exp(-(activation_energy + vibr_energy_product) / kU) * exp(i_e * (1./kT + 1./kU));
    }

    return res * p_Z_electron(T, electron_energy, statistical_weight, n_electron_levels);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::C_aliat(	double T, 
					const arma::vec &electron_energy, 
					const arma::Col<unsigned long> &statistical_weight, 
					int n_electron_levels,
                                     	const std::vector<arma::vec> &vibr_energy, 	
					const std::vector<int> &num_vibr_levels,
                                     	double vibr_energy_product, 
					double activation_energy, 
					int i, 
					int e) { // U = infty
    double tmp=0, res=0;
    double kT = K_CONST_K * T;
    double i_e;

    for (int ee=0; ee<n_electron_levels; ee++) {
      tmp = 0;
      for (int k = 0; k < num_vibr_levels[ee]; k++) {
        i_e = vibr_energy[ee][k] + electron_energy[ee];
        if (i_e > activation_energy + vibr_energy_product) {
          tmp += statistical_weight[ee] * exp((activation_energy + vibr_energy_product - i_e) / kT);
        } else {
          tmp += statistical_weight[ee];
        }
      }
      res += tmp / p_Z_vibr_eq(T, vibr_energy[ee]);
    }
    res = 1./res;
    i_e = vibr_energy[e][i] + electron_energy[e];
    if (i_e > activation_energy + vibr_energy_product) {
      res *= exp((activation_energy + vibr_energy_product) / kT);
    } else {
      res *= exp(i_e / kT);
    }

    return res * p_Z_electron(T, electron_energy, statistical_weight, n_electron_levels);
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Returns the maximum vibrational level possible for the molecule
  int kappa::Approximation::max_i(double T, const kappa::Molecule &molecule) {

    if (molecule.anharmonic_spectrum == false) {
      return molecule.num_vibr_levels[0]; 
    } else {
      int i_star = p_max_i(T, molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0], K_CONST_C * molecule.vibr_frequency[0], molecule.vibr_we_xe[0] / molecule.vibr_frequency[0]);
      if (i_star < 1) {
        std::string error_string = "No Treanor distribution possible for such values of T, T1";
        throw kappa::ModelParameterException(error_string.c_str());      
      } else if (i_star >= molecule.num_vibr_levels[0] - 1) {
        return molecule.num_vibr_levels[0] - 1;  
      } else {
        return p_max_i(T, molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0], K_CONST_C * molecule.vibr_frequency[0], molecule.vibr_we_xe[0] / molecule.vibr_frequency[0]);
      }
    }
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  // averaged rotational energy
  double kappa::Approximation::p_avg_rot_energy(double T, const arma::vec &rot_energy, int num_rot_levels, int rot_symmetry) {

    return arma::dot((2 * arma::linspace<arma::vec>(0, num_rot_levels - 1, num_rot_levels) + 1) % rot_energy, arma::exp(-rot_energy / (K_CONST_K * T))) 
          / (rot_symmetry * p_Z_rot(T, rot_energy, num_rot_levels, rot_symmetry));
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // squared averaged rotational energy 
  double kappa::Approximation::p_avg_rot_energy_sq(double T, const arma::vec &rot_energy, int num_rot_levels, int rot_symmetry) {

    return arma::dot((2 * arma::linspace<arma::vec>(0, num_rot_levels - 1, num_rot_levels) + 1) % rot_energy % rot_energy, arma::exp(-rot_energy / (K_CONST_K * T))) 
           / (rot_symmetry * p_Z_rot(T, rot_energy, num_rot_levels, rot_symmetry));
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::p_c_tr(double T, double mass) {
    return 1.5 * K_CONST_K / mass;
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::p_c_rot(double T, double mass, const arma::vec &rot_energy, int num_rot_levels, int rot_symmetry) {

    double are = p_avg_rot_energy(T, rot_energy, num_rot_levels, rot_symmetry);

    // eq. 15, On the applicability of simplified state-to-state models of transport coefficients
    return (p_avg_rot_energy_sq(T, rot_energy, num_rot_levels, rot_symmetry) - are * are) / (K_CONST_K * T * T * mass);

    // eq. 1.24 -- Hyp. rigid rotator model and kT >> h*h/(8π*π*Ic) 
    //return K_CONST_K * mass;

  }  

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::vel_avg_vv(double rel_vel, double coll_mass, double Delta_E_vibr) {

    // vibr_energy_before1 - vibr_energy_after1 + vibr_energy_before2 - vibr_energy_after2
    double rel_vel_after_sq = Delta_E_vibr * (2.0 / coll_mass) + rel_vel * rel_vel;
    if (rel_vel_after_sq < 0) {
      return -1;
    } else {
      return 0.5 * (rel_vel + sqrt(rel_vel_after_sq));
    }
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::vel_avg_vt(double rel_vel, double coll_mass, double Delta_E_vibr) {

    double rel_vel_after_sq = Delta_E_vibr * (2.0 / coll_mass) + rel_vel * rel_vel;
    if (rel_vel_after_sq < 0) {
     return -1;
    } else {
      return 0.5 * (rel_vel + sqrt(rel_vel_after_sq));
    }
  }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Vibrational energy transfer rates using a Forced Harmonic Oscillator (FHO) model.
  // Journal of Thermodynamics and Heat Transfer, Vol. 12, No. 1, January-March 1998.
  // I. V. Adamovich, S. Macheret, J. W. Rich, C. E. Treanor.
  double kappa::Approximation::probability_VV_FHO(	double rel_vel, 
 							double coll_mass, 
	                                        	double Delta_E_vibr,
			                        	int i, 
							int k, 
							int delta_i,
			                        	double omega1, 
							double omega2, 
							double omega, 
							double alpha_FHO) {

    double avg_vel = vel_avg_vv(rel_vel, coll_mass, Delta_E_vibr);
    if (avg_vel <= 0.0) {
      return 0.0;
    }

    double ns12;
    int abs_delta;

    // eq. 12, 14
    if (delta_i > 0) { // i -> i + delta_i>0, k -> k - delta_i
      abs_delta = delta_i;
      ns12 = (kappa::factorial_table[i + delta_i] / kappa::factorial_table[i]) * (kappa::factorial_table[k] / kappa::factorial_table[k - delta_i]);
    } else { // i -> i + delta_i<0, k -> k - delta_i
      abs_delta = -delta_i;
      ns12 = (kappa::factorial_table[i] / kappa::factorial_table[i + delta_i]) * (kappa::factorial_table[k - delta_i] / kappa::factorial_table[k]);
    }

    // eq. 9
    double rhoxi = 2.0 * omega / (alpha_FHO * rel_vel);
    if (rhoxi < 0.00001) {
      rhoxi = 0.25 * K_CONST_SVV_FHO * alpha_FHO * avg_vel * alpha_FHO * avg_vel / (omega1 * omega2);
    } else {
      rhoxi = K_CONST_SVV_FHO * omega * omega / (omega1 * omega2 * pow(sinh(rhoxi), 2.0));
    }

    // eq. 17
    return pow(rhoxi, abs_delta) * ns12 * exp(-2.0 * rhoxi * pow(ns12, 1.0 / abs_delta) / (abs_delta + 1.0)) / 
                                          (kappa::factorial_table[abs_delta] * kappa::factorial_table[abs_delta]);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Vibrational energy transfer rates using a Forced Harmonic Oscillator (FHO) model.
  // Journal of Thermodynamics and Heat Transfer, Vol. 12, No. 1, January-March 1998.
  // I. V. Adamovich, S. Macheret, J. W. Rich, C. E. Treanor.
  double kappa::Approximation::probability_VT_FHO(	double rel_vel, 
                                              		double coll_mass, 
							double osc_mass, 
							double Delta_E_vibr,
		                                	int i,	 
							int delta_i,
		                                	double omega,
							double ram,
							double alpha_FHO,
							double E_FHO, 
							double svt_FHO) {

    double avg_vel = vel_avg_vt(rel_vel, coll_mass, Delta_E_vibr);
   
    if (avg_vel <= 0.0) {
      return 0.0;
    }

    double ns;
    int s;

    if (delta_i > 0) {
      s = delta_i;
      ns = kappa::factorial_table[i + delta_i] / kappa::factorial_table[i];
    } else {
      s = -delta_i;
      ns = kappa::factorial_table[i] / kappa::factorial_table[i + delta_i];
    }

    double phi = (2. / K_CONST_PI) * atan(sqrt((2. * E_FHO) / (coll_mass * avg_vel * avg_vel)));
    double eps = (cosh((1. + phi) * K_CONST_PI * omega / (alpha_FHO * avg_vel)));
    eps *= K_CONST_PI * ram * coll_mass / (sinh(2. * K_CONST_PI * omega / (alpha_FHO * avg_vel)) * alpha_FHO);
    eps = eps * eps;
    eps *= 8. * svt_FHO * omega / (osc_mass * K_CONST_H);

    return ns * pow(eps, s) * exp(- 2. * pow(ns, 1. / s) * eps / (s + 1.)) / (kappa::factorial_table[s] * kappa::factorial_table[s]);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::min_vel_diss(double coll_mass, double diss_energy, double vibr_energy) {
    return sqrt(2 * (diss_energy - vibr_energy) / coll_mass);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::min_vel_diss_ILT_N2N(double coll_mass, double vibr_energy, double i) {

    double c1;
    if (i <= 8) {
      c1 = 1.786e-18;
    } else if (i<=34) {
      c1 = 1.71e-18;
    } else if (i<=52) {
      c1 = 1.68e-18;
    } else {
      c1 = 1.66e-18;
    }
    return sqrt(2 * (c1 - vibr_energy) / coll_mass);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::min_vel_diss_ILT_O2O(double coll_mass, double i) {

    double c1 = 0.3867 * i * i * i - 2.7425 * i * i - 1901.9 * i + 61696;
    return sqrt(2 * c1 * K_CONST_K / coll_mass);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::convert_vibr_ladder_N2(double vibr_energy) {

    double te = 1000 * vibr_energy / K_CONST_EV;
    return 0.58142 - 0.0039784 * te + 1.3127e-5 * te * te - 9.877e-9 * te * te * te + 3.8664e-12 * te * te * te * te - 8.4268e-16 * te * te * te * te * te
           + 1.0326e-19 * te * te * te * te * te * te - 6.649e-24 * te * te * te * te * te * te * te + 1.7506e-28 * te * te * te * te * te * te * te * te; 
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::convert_vibr_ladder_O2(double vibr_energy) {

    double te = vibr_energy / K_CONST_EV;
    return 0.4711499685975582 + 2.23383944 * te + 2.17784653 * te * te + 0.61279172 * te * te * te - 1.38809821 * te * te * te * te + 0.66364689 * 
           te * te * te * te * te - 0.13147542 * te * te * te * te * te * te + 0.00950554 * te * te * te * te * te * te * te; 
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  double kappa::Approximation::p_probability_diss(double rel_vel, double coll_mass, double diss_energy, double vibr_energy, bool center_of_mass) {

    double energy = vibr_energy + rel_vel * rel_vel * coll_mass / 2;

    if (energy < diss_energy) {
      return 0.0;
    } else if (center_of_mass) {
      return 1.0 - diss_energy / energy;
    } else {
      return 1.0;
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::crosssection_elastic_RS(double diameter) {
    return K_CONST_PI * diameter * diameter;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::crosssection_elastic_VSS(double rel_vel, double coll_mass, double vss_c, double vss_omega) {
    return vss_c * pow(coll_mass * rel_vel * rel_vel / (2 * K_CONST_K), -vss_omega);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::crosssection_elastic_VSS(double rel_vel, double vss_c_cs, double vss_omega) {
    return vss_c_cs * pow(rel_vel, 1. - 2. * vss_omega);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::crosssection_VT_FHO_RS(	double rel_vel, 
							double coll_mass, 
							double diameter, 
							double osc_mass,
                                                    	double Delta_E_vibr, 
							int i, 
							int delta_i,
 							double omega,
                                                    	double ram, 
							double alpha_FHO,
							double E_FHO, 
							double svt_FHO) {
    if (delta_i < 0) {
      return probability_VT_FHO(rel_vel, coll_mass, osc_mass, Delta_E_vibr, i, delta_i, omega, ram, alpha_FHO, E_FHO, svt_FHO) * crosssection_elastic_RS(diameter);
    } else {
      return 0.; // TODO rewrite
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::crosssection_VT_FHO_VSS(	double rel_vel, 
							double coll_mass, 	
							double vss_c_cs, 
							double vss_omega, 
							double osc_mass,
                                                      	double Delta_E_vibr, 
							int i, 
							int delta_i,
                                                      	double omega, 
							double ram, 
							double alpha_FHO, 
							double E_FHO, 
							double svt_FHO) {
    if (delta_i < 0) {
      return probability_VT_FHO(rel_vel, coll_mass, osc_mass, Delta_E_vibr, i, delta_i, omega, ram, alpha_FHO, E_FHO, svt_FHO) 
             * crosssection_elastic_VSS(rel_vel, vss_c_cs, vss_omega);
    } else {
      return 0.; // TODO rewrite
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::crosssection_diss_RS(double rel_vel, double coll_mass, double diameter, double diss_energy, double vibr_energy, bool center_of_mass) {
   return crosssection_elastic_RS(diameter) * p_probability_diss(rel_vel, coll_mass, diss_energy, vibr_energy, center_of_mass);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::crosssection_diss_VSS(double rel_vel, double coll_mass, double vss_c_cs, double vss_omega, double diss_energy, double vibr_energy, bool center_of_mass) {
    return crosssection_elastic_VSS(rel_vel, vss_c_cs, vss_omega) * p_probability_diss(rel_vel, coll_mass, diss_energy, vibr_energy, center_of_mass);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::crosssection_diss_ILT_N2N(double rel_vel, double coll_mass, double vibr_energy, double i) {

    double t = coll_mass * rel_vel * rel_vel / (2. * K_CONST_K);
    double c1;
    if (i <= 8) {
      c1 = 1.786e-18;
    } else if (i<=34) {
      c1 = 1.71e-18;
    } else if (i<=52) {
      c1 = 1.68e-18;
    } else {
      c1 = 1.66e-18;
    }

    if (t < (c1 - vibr_energy) / K_CONST_K) {
      return 0;
    }

    return sqrt(K_CONST_PI * coll_mass / (8 * K_CONST_K)) * 7.16e-2 * (5.24e-19 * i * i * i - 7.41e-17 * i * i + 6.42e-15 * i + 7.3e-14) 
           * pow(t - (c1 - vibr_energy) / K_CONST_K, 0.25) / (0.91 * t);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::crosssection_diss_ILT_O2O(double rel_vel, double coll_mass, double vibr_energy, double i) {

    double t = coll_mass * rel_vel * rel_vel / (2 * K_CONST_K);
    double c1 = 0.3867 * i * i * i - 2.7425 * i * i - 1901.9 * i + 61696;
    double c2;

    if (t < c1) {
      return 0;
    }

    if (i <= 31) {
      c2 = 1.63e-9 * i * i * i - 1.25e-7 * i * i + 3.24e-6 * i + 7.09e-5;
    } else if (i<=37) {
      c2 = - 6.67e-6 * i * i + 4.65e-4 * i - 7.91e-3;
    } else {
      c2 = 7.83e-7 * i * i * i * i - 1.31e-4 * i * i * i + 8.24e-3 * i * i - 0.23 * i + 2.4049;
    }
  
    return sqrt(K_CONST_PI * coll_mass / (8 * K_CONST_K)) * 1.53e-10 * c2 * pow(t - c1, 0.4) / (0.89 * t);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // The influence of vibrational state-resolvent transport coefficients on the wave propagation in diatomic gases.
  double kappa::Approximation::p_Z_coll(double T, double n, double coll_mass, double diameter) {

   // collision frequency, eq. 82a
   return 4 * K_CONST_PI * n * diameter * diameter * sqrt(K_CONST_K * T / (2 * K_CONST_PI * coll_mass));
  } 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Transport properties of a reacting gas mixture with strong vibrational and chemical nonequilibrium, Kustova and Nagnibeda, 1998.
  double kappa::Approximation::integral_VT_FHO_RS(	double T, 
							int degree, 
							double coll_mass, 
							double diameter, 
							double osc_mass,
							double Delta_E_vibr, 
							int i, 
							int delta_i, 
 							double omega,
							double ram, 
							double alpha_FHO,
							double E_FHO, 
							double svt_FHO) {

    double conversion = sqrt(2 * K_CONST_K * T / coll_mass);

    auto integrand = [T, degree, coll_mass, diameter, osc_mass, Delta_E_vibr, i, delta_i, omega, ram, alpha_FHO, E_FHO, svt_FHO, conversion]
    (double g) {return pow(g, 2 * degree + 3) *
    crosssection_VT_FHO_RS(conversion * g, coll_mass, diameter, osc_mass, Delta_E_vibr, i, delta_i, omega, ram, alpha_FHO, E_FHO, svt_FHO) * exp(-g * g); };

    return sqrt(K_CONST_K * T / (2 * K_CONST_PI * coll_mass)) * integrate_semi_inf(integrand, 0);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::integral_VT_FHO_VSS(	double T, 
							int degree, 
							double coll_mass, 
							double vss_c_cs, 
							double vss_omega, 
							double osc_mass,
							double Delta_E_vibr, 
							int i, 
							int delta_i,
							double omega, 
							double ram, 						
							double alpha_FHO,
							double E_FHO,
							double svt_FHO) {

    double conversion = sqrt(2 * K_CONST_K * T / coll_mass);
 
    auto integrand = [T, degree, coll_mass, vss_c_cs, vss_omega, osc_mass, Delta_E_vibr, i, delta_i, omega, ram, alpha_FHO, E_FHO, svt_FHO, conversion]
    (double g) {return pow(g, 2 * degree + 3) *
    crosssection_VT_FHO_VSS(conversion * g, coll_mass, vss_c_cs, vss_omega, osc_mass, Delta_E_vibr, i, delta_i, omega, ram, alpha_FHO, E_FHO, svt_FHO) * exp(-g * g);};
  
    return sqrt(K_CONST_K * T / (2 * K_CONST_PI * coll_mass)) * integrate_semi_inf(integrand, 0);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::integral_diss_RS(double T, 
						int degree, 
						double coll_mass, 
						double diameter, 
						double diss_energy, 
						double vibr_energy, 
						bool center_of_mass) {

    double conversion = sqrt(2 * K_CONST_K * T / coll_mass);

    auto integrand = [T, degree, coll_mass, diameter, diss_energy, vibr_energy, center_of_mass, conversion](double g) 
                      {return pow(g, 2 * degree + 3) * 
                      crosssection_diss_RS(conversion * g, coll_mass, diameter, diss_energy, vibr_energy, center_of_mass) * exp(-g * g); };

    return sqrt(K_CONST_K * T / (2 * K_CONST_PI * coll_mass)) * integrate_semi_inf(integrand, min_vel_diss(coll_mass, diss_energy, vibr_energy) / conversion);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::integral_diss_VSS(	double T, 
							int degree, 
							double coll_mass, 	
							double vss_c_cs, 
							double vss_omega, 
							double diss_energy, 
							double vibr_energy, 
							bool center_of_mass) {

    double conversion = sqrt(2 * K_CONST_K * T / coll_mass);
    auto integrand = [T, degree, coll_mass, vss_c_cs, vss_omega, diss_energy, vibr_energy, center_of_mass, conversion](double g) {return pow(g, 2 * degree + 3) *
    crosssection_diss_VSS(conversion * g, coll_mass, vss_c_cs, vss_omega, diss_energy, vibr_energy, center_of_mass) * exp(-g * g); };

    return sqrt(K_CONST_K * T / (2 * K_CONST_PI * coll_mass)) * integrate_semi_inf(integrand, min_vel_diss(coll_mass, diss_energy, vibr_energy) / conversion);
  }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::integral_diss_ILT_N2N(double T, int degree, double coll_mass, double vibr_energy, double i) {

    double conversion = sqrt(2 * K_CONST_K * T / coll_mass);
    auto integrand = [T, degree, coll_mass, vibr_energy, i, conversion](double g) {return pow(g, 2 * degree + 3) *
    crosssection_diss_ILT_N2N(conversion * g, coll_mass, vibr_energy, i) * exp(-g * g); };

    return sqrt(K_CONST_K * T / (2 * K_CONST_PI * coll_mass)) * integrate_semi_inf(integrand, min_vel_diss_ILT_N2N(coll_mass, vibr_energy, i) / conversion);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::integral_diss_ILT_O2O(double T, int degree, double coll_mass, double vibr_energy, double i) {

    double conversion = sqrt(2 * K_CONST_K * T / coll_mass);
    auto integrand = [T, degree, coll_mass, vibr_energy, i, conversion](double g) {return pow(g, 2 * degree + 3) *
    crosssection_diss_ILT_O2O(conversion * g, coll_mass, vibr_energy, i) * exp(-g * g); };

    return sqrt(K_CONST_K * T / (2 * K_CONST_PI * coll_mass)) * integrate_semi_inf(integrand, min_vel_diss_ILT_O2O(coll_mass, i) / conversion);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // The influence of vibrational state-resolvent transport coefficients on the wave propagation in diatomic gases.
  // average probability of the VT transition M(1) + M -> M(0) + M
  double kappa::Approximation::P_SSH_VT_10(double T, double coll_mass, double omega_e, double epsilon, double diameter, double r_e) {

    double kT = K_CONST_K * T;
    // eq. 86a
    double alpha = 17.5 / diameter; // diameter is the distance at which the LJ potential is zero
    double omega = 2 * K_CONST_PI * K_CONST_C * omega_e; // angular frequency of the oscillator
    // eq. 86b
    double chi = pow(K_CONST_PI * K_CONST_PI * coll_mass * omega * omega / (2 * alpha * alpha * kT), 0.3333333); // TODO 4 or 2?
    // eq. 86c
    double r = diameter * pow(0.5 * sqrt(1 + chi * kT / epsilon) + 0.5, -0.1666666); // coordinate of the turning point
    // eq. 85
    double Z0 = alpha * r_e * alpha * r_e * exp( -(3 * alpha * r_e * r_e) / (8 * r) ); // steric factor
    // eq. 83
    return 1.294 * pow(r / diameter, 2) / (Z0 * alpha * alpha * K_CONST_HBAR * (1 + 1.1 * epsilon / kT)) * 
           4 * K_CONST_PI * K_CONST_PI * coll_mass * omega * sqrt(4 * K_CONST_PI * chi / 3) * 
           exp(-3*chi + (K_CONST_HBAR * omega / (2 * kT)) + epsilon/kT); // TODO second 2 or 4?
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // The influence of vibrational state-resolvent transport coefficients on the wave propagation in diatomic gases.
  // average probability of the VV transition M(1) + M(0) -> M(0) + M(1)
  double kappa::Approximation::P_SSH_VV_01(double T, double omega_e, double epsilon, double osc_mass, double diameter, double r_e) {

    // eq. 86a
    double alpha = 17.5 / diameter; // diameter is the distance at which the LJ potential is zero
    double lambda = 0.5; // for diatomic homonuclear molecules
    double omega = 2 * K_CONST_PI * K_CONST_C * omega_e; //angular frequency of the oscillator
    // eq. 84
    return pow(lambda, 4) * 4 * K_CONST_K * T / osc_mass * alpha * alpha / omega / omega; 
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // The influence of vibrational state-resolvent transport coefficients on the wave propagation in diatomic gases.
  // Compute rate coefficients of the VT transition using Schwartz-Slavsky-Herzfeld (SSH) theory for the harmonic oscillator model
  double kappa::Approximation::k_VT_SSH(double T, int i, double coll_mass, double diameter, double omega_e, double epsilon, double r_e) {

    // eq. 80b
    return p_Z_coll(T, 1., coll_mass, diameter) * i * P_SSH_VT_10(T, coll_mass, omega_e, epsilon, diameter, r_e);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // The influence of vibrational state-resolvent transport coefficients on the wave propagation in diatomic gases.
  // Compute rate coefficients of the VV transition using Schwartz-Slavsky-Herzfeld (SSH) theory
  double kappa::Approximation::k_VV_SSH(double T, int i, int k, double coll_mass, double diameter, double omega_e, double epsilon, double r_e) {

    // eq. 81
    return p_Z_coll(T, 1., coll_mass, diameter) * (i+1)*(k+1) * P_SSH_VV_01(T, coll_mass, omega_e, epsilon, diameter, r_e);
  } 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // The influence of vibrational state-resolvent transport coefficients on the wave propagation in diatomic gases.
  // Compute rate coefficients of the VT transition using Schwartz-Slavsky-Herzfeld (SSH) theory for the anharmonic oscillator model
  double kappa::Approximation::k_VT_SSH(double T, 
 					int i, 			// vibrational level
					double coll_mass, 	// collisional mass
					double diameter, 	// diameter is the distance at which the (LJ) potential is zero
					double omega_e,         
					double epsilon, 	// depth of the potential well
					double r_e, 		// internuclear distance of molecules
					double Delta_E_vibr, 	
					double vibr_energy_1) {

    double alpha = 17.5 / diameter;
    double gamma0 = K_CONST_PI * sqrt(coll_mass / (2 * K_CONST_K * T)) / (alpha * K_CONST_HBAR);
    double gammai = gamma0 * (vibr_energy_1 - 2 * i * Delta_E_vibr);
    gamma0 *= vibr_energy_1;
    double delta_VT = Delta_E_vibr / vibr_energy_1;
  
    if (gammai >= 20) {
      delta_VT *= 4 * pow(gamma0, 2./3.);
    } else {
      delta_VT *= (4./3) * gamma0;
    }
  
    return k_VT_SSH(T, i, coll_mass, diameter, omega_e, epsilon, r_e) * exp(i * delta_VT) * exp(-i * Delta_E_vibr / (K_CONST_K * T));
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::k_VT_FHO_RS(	double T, 
						double coll_mass, 
						double diameter, 
						double osc_mass,
                                         	double Delta_E_vibr, 
						int i, 
						int delta_i, 
						double omega,
						double ram, 
						double alpha_FHO,
						double E_FHO,
						double svt_FHO) { 

    // Transport properties of a reacting gas mixture with strong vibrational and chemical nonequilibrium. Kustova and Nagnibeda, 1998. Eq. 65.
    return 8 * integral_VT_FHO_RS(T, 0, coll_mass, diameter, osc_mass, Delta_E_vibr, i, delta_i, omega, ram, alpha_FHO, E_FHO, svt_FHO);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::k_VT_FHO_VSS(	double T, 
						double coll_mass, 
						double vss_c_cs, 
						double vss_omega, 
						double osc_mass,
                                         	double Delta_E_vibr, 
						int i, 
						int delta_i,
                                          	double omega,  
						double ram,
						double alpha_FHO,
						double E_FHO, 
						double svt_FHO) {

    return 8 * integral_VT_FHO_VSS(T, 0, coll_mass, vss_c_cs, vss_omega, osc_mass, Delta_E_vibr, i, delta_i, omega, ram, alpha_FHO, E_FHO, svt_FHO);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // VV rate coefficients in N2
  double kappa::Approximation::k_VV_Billing_N2N2(double T, int i, int k) {

    double tmp = exp(-6.8 * abs(i - 1 - k) / sqrt(T));
    return 2.5e-20 * i * (k + 1) * pow(T / 300., 1.5) * tmp * (1.5 - 0.5 * tmp);
  
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // VT rate coefficients in N2-N2
  double kappa::Approximation::k_VT_Billing_N2N2(double T, int i) { 

    // eq. 6.9 FIXME why 1e-6?
    return 1e-6 * i * exp(-3.24093 - 140.69597 * pow(T, -0.2)) * exp((0.26679 - 6.99237e-5 * T + 4.70073e-9 * T * T) * (i-1));
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // VT rate coefficients in N2-N
  double kappa::Approximation::k_VT_Billing_N2N(double T, int i, int delta_i) { 

    // eq. 6.10 FIXME why 1e-6?
    return 1e-6 * exp(-25.708 - 5633.1543 / T + (0.1554 - 111.3426 / T) * delta_i + (0.0054 - 2.189 / T) * delta_i * delta_i
               + i * (0.0536 + 122.4835 / T - (0.0013 - 4.2365 / T) * delta_i + (-1.197e-4 + 0.0807 / T) * delta_i * delta_i));
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Phys4Entary database for VT exchange
  double kappa::Approximation::k_VT_N2N_p4e(double T, int i, int delta_i) {

    double a[5]={0,0,0,0,0};
    int j=0;

    delta_i = -delta_i;
    if (delta_i==1) {
      for (j=0;j<=4;j++) {
        a[j]=p4e_n2n_vt_bij1[j][0]+p4e_n2n_vt_bij1[j][1]*i+p4e_n2n_vt_bij1[j][2]*(i*i)+p4e_n2n_vt_bij1[j][3]*(i*i*i)+p4e_n2n_vt_bij1[j][4]*log(i);
      }
      return exp(a[0]+a[1]/T+a[2]/(pow(T,2))+a[3]/(pow(T,3))+a[4]*log(T));
    } else if (delta_i>=2 && delta_i<=5){
      for (j=0;j<=4;j++) {
        a[j] = (p4e_n2n_vt_bij2[j][0]+p4e_n2n_vt_cij2[j][0]*delta_i)
             + (p4e_n2n_vt_bij2[j][1]+p4e_n2n_vt_cij2[j][1]*delta_i)*i
             + (p4e_n2n_vt_bij2[j][2]+p4e_n2n_vt_cij2[j][2]*delta_i)*(i*i)
             + (p4e_n2n_vt_bij2[j][3]+p4e_n2n_vt_cij2[j][3]*delta_i)*(i*i*i)
             + (p4e_n2n_vt_bij2[j][4]+p4e_n2n_vt_cij2[j][4]*delta_i)*log(i);
      }
      return 0.000001 * exp(a[0]+a[1]/T+a[2]/(pow(T,2))+a[3]/(pow(T,3))+a[4]*log(T));
    } else if (delta_i>=6 && delta_i<=10) {
      for (j=0;j<=4;j++) {
        a[j] = (p4e_n2n_vt_bij3[j][0]+p4e_n2n_vt_cij3[j][0]*delta_i)
             + (p4e_n2n_vt_bij3[j][1]+p4e_n2n_vt_cij3[j][1]*delta_i)*i
             + (p4e_n2n_vt_bij3[j][2]+p4e_n2n_vt_cij3[j][2]*delta_i)*(i*i)
             + (p4e_n2n_vt_bij3[j][3]+p4e_n2n_vt_cij3[j][3]*delta_i)*(i*i*i)
             + (p4e_n2n_vt_bij3[j][4]+p4e_n2n_vt_cij3[j][4]*delta_i)*log(i);
      }
      return 0.000001 * exp(a[0]+a[1]/T+a[2]/(pow(T,2))+a[3]/(pow(T,3))+a[4]*log(T));
    } else if (delta_i>=11 && delta_i<=20){
      for (j=0;j<=4;j++) {
        a[j] = (p4e_n2n_vt_bij4[j][0]+p4e_n2n_vt_cij4[j][0]*delta_i)
             + (p4e_n2n_vt_bij4[j][1]+p4e_n2n_vt_cij4[j][1]*delta_i)*i
             + (p4e_n2n_vt_bij4[j][2]+p4e_n2n_vt_cij4[j][2]*delta_i)*(i*i)
             + (p4e_n2n_vt_bij4[j][3]+p4e_n2n_vt_cij4[j][3]*delta_i)*(i*i*i)
             + (p4e_n2n_vt_bij4[j][4]+p4e_n2n_vt_cij4[j][4]*delta_i)*log(i);
      }
      return 0.000001 * exp(a[0]+a[1]/T+a[2]/(pow(T,2))+a[3]/(pow(T,3))+a[4]*log(T));
    } else if (delta_i>=21 && delta_i<=30){
      for (j=0;j<=4;j++) {
        a[j] = (p4e_n2n_vt_bij5[j][0]+p4e_n2n_vt_cij5[j][0]*delta_i)
             + (p4e_n2n_vt_bij5[j][1]+p4e_n2n_vt_cij5[j][1]*delta_i)*i
             + (p4e_n2n_vt_bij5[j][2]+p4e_n2n_vt_cij5[j][2]*delta_i)*(i*i)
             + (p4e_n2n_vt_bij5[j][3]+p4e_n2n_vt_cij5[j][3]*delta_i)*(i*i*i)
             + (p4e_n2n_vt_bij5[j][4]+p4e_n2n_vt_cij5[j][4]*delta_i)*log(i);
      }
      return 0.000001 * exp(a[0]+a[1]/T+a[2]/(pow(T,2))+a[3]/(pow(T,3))+a[4]*log(T));
    }  else {
       for (j=0;j<=3;j++) {
         a[j] = (p4e_n2n_vt_bij6[j][0]+p4e_n2n_vt_cij6[j][0]*delta_i)
              + (p4e_n2n_vt_bij6[j][1]+p4e_n2n_vt_cij6[j][0]*delta_i)*i
              + (p4e_n2n_vt_bij6[j][2]+p4e_n2n_vt_cij6[j][0]*delta_i)*(i*i);
       };
       return 0.000001 * exp(a[0]+a[1]/T+a[2]/(pow(T,4))+a[3]/log(T));
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::k_VT_O2O_p4e(double T, int k, int delta_i) {

    double Ai[3];
    double Bij[3][5];
    int i,j;

    delta_i = -delta_i;
    
    if (delta_i==1) {
      for (i=0;i<=2;i++)
        for (j=0;j<=4;j++)
          Bij[i][j] = p4e_o2o_vt_Cijk1[i][j][0]
                    + p4e_o2o_vt_Cijk1[i][j][1]*log(delta_i)
                    + p4e_o2o_vt_Cijk1[i][j][2]*delta_i*exp(-delta_i)
                    + p4e_o2o_vt_Cijk1[i][j][3]*delta_i*delta_i;
     } else if (delta_i>1 && delta_i<=10) {
       for (i=0;i<=2;i++)
         for (j=0;j<=4;j++)
           Bij[i][j] = p4e_o2o_vt_Cijk2[i][j][0]
                     + p4e_o2o_vt_Cijk2[i][j][1] * log(delta_i)
                     + p4e_o2o_vt_Cijk2[i][j][2] * delta_i * exp(-delta_i)
                     + p4e_o2o_vt_Cijk2[i][j][3] * delta_i * delta_i;
     } else if (delta_i>10 && delta_i<=20) {
       for (i=0;i<=2;i++)
         for (j=0;j<=4;j++)
           Bij[i][j] = p4e_o2o_vt_Cijk3[i][j][0]
                     + p4e_o2o_vt_Cijk3[i][j][1] * log(delta_i)
                     + p4e_o2o_vt_Cijk3[i][j][2] * delta_i * exp(-delta_i)
                     + p4e_o2o_vt_Cijk3[i][j][3] * delta_i * delta_i;
     } else if (delta_i>20 && delta_i<=30) {
       for (i=0;i<=2;i++)
         for (j=0;j<=4;j++)
           Bij[i][j] = p4e_o2o_vt_Cijk4[i][j][0]
                     + p4e_o2o_vt_Cijk4[i][j][1] * log(delta_i)
                     + p4e_o2o_vt_Cijk4[i][j][2] * delta_i * exp(-delta_i)
                     + p4e_o2o_vt_Cijk4[i][j][3] * delta_i * delta_i;
     } else {
       return -1.;
     }   

     for (i=0;i<=2;i++) {
       Ai[i] = Bij[i][0] + Bij[i][1] * log(k) + (Bij[i][2] + Bij[i][3] * k + Bij[i][4] * k * k) / (1E+21 + exp(k));
     }
     return 0.000001 * (1./27 + exp(Ai[0] + Ai[1] / log(T) + Ai[2] * log(T)));
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Empirical Arrhenius law
  // energy -> dissociation energy
  // A, n parameters from [?] <- different models
  double kappa::Approximation::k_Arrhenius(double T, double arrhenius_A, double arrhenius_n, double energy) {
    return arrhenius_A * pow(T, arrhenius_n) * exp(-energy / (K_CONST_K * T));
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::k_exch_WRFP(double T, double vibr_energy, double activation_energy, double alpha_exch, double beta_exch, double A_exch, double n_exch) {

    if (activation_energy > alpha_exch * vibr_energy) {
      return A_exch * pow(T, n_exch) * exp(-(activation_energy - alpha_exch * vibr_energy) / (K_CONST_K * T * beta_exch)); 
    } else {
      return A_exch * pow(T, n_exch);
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::k_diss_RS(double T, double coll_mass, double diameter, double diss_energy, double vibr_energy, bool center_of_mass) {
    return 8. * integral_diss_RS(T, 0, coll_mass, diameter, diss_energy, vibr_energy, center_of_mass);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::k_diss_VSS(double T, double coll_mass, double vss_c_cs, double vss_omega, double diss_energy, double vibr_energy, bool center_of_mass) {
    return 8. * integral_diss_VSS(T, 0, coll_mass, vss_c_cs, vss_omega, diss_energy, vibr_energy, center_of_mass);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::k_diss_ILT_N2N(double T, double coll_mass, double vibr_energy, double i) {
    return 8. * integral_diss_ILT_N2N(T, 0, coll_mass, vibr_energy, i);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::k_diss_ILT_O2O(double T, double coll_mass, double vibr_energy, double i) {
    return 8. * integral_diss_ILT_O2O(T, 0, coll_mass, vibr_energy, i);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::diss_rate_N2N_p4e(double T, double i) {

    double a[5]={0,0,0,0,0};
    for (int j=0; j<=4; j++){
      a[j] = p4e_n2n_diss_bjk[0][j]
           + p4e_n2n_diss_bjk[1][j] * i
           + p4e_n2n_diss_bjk[2][j] * (i * i)
           + p4e_n2n_diss_bjk[3][j] * (i * i * i);
    }
    return 0.000001 * exp(a[0] + a[1] / T + a[2] / (T * T) + a[3] / (T * T * T) + a[4] * log(T));
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::diss_rate_O2O_p4e(double T, double i) {

    double a[3]={0,0,0};
    for (int j=0;j<=2;j++) {
      a[j] = p4e_o2o_diss_bjk[j][0]
          + p4e_o2o_diss_bjk[j][1]*i
          + p4e_o2o_diss_bjk[j][2]*(i*i)
          + p4e_o2o_diss_bjk[j][3]*(i*i*i)
          + p4e_o2o_diss_bjk[j][4]*(i*i*i*i)
          + p4e_o2o_diss_bjk[j][5]*(i*i*i*i*i)
          + p4e_o2o_diss_bjk[j][6]*(i*i*i*i*i*i)
          + p4e_o2o_diss_bjk[j][7]*(i*i*i*i*i*i*i);
    }
    return 0.000001 * exp(a[0] + a[1] / T + a[2] * log(T));
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::p_rot_collision_number_parker(double T, double xi_inf, double epsilon) {

    double kTe = T * K_CONST_K / epsilon; // epsilon/k=97.5 for N2, 107.4 for O2, k=1.3807e-23
    return xi_inf /  ( 1 + pow(K_CONST_PI, 1.5) * pow(kTe, -0.5) / 2 + (K_CONST_PI * K_CONST_PI / 4 + 2) / kTe + pow(K_CONST_PI, 1.5) * pow(kTe, -1.5));
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::p_vibr_relaxation_time_MW(double T, double n, double char_vibr_temp, double coll_mass) {

    // 9.8692326671601e-6 is the conversion factor from pascals to atmospheres
    double md = coll_mass * K_CONST_NA * 1000; // reduced molar mass, g/mol , N_Av avogadro number
    return pow(10.0, ((5e-4 * pow(md, 0.5)
           * pow(char_vibr_temp, 1.333333333) 
           * (pow(T, -0.3333333) - 0.015 * pow(md, 0.25))) - 8)) / (n * K_CONST_K * T * 9.8692326671601e-6);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::p_vibr_relaxation_time_Park_corr(double T, double n, double coll_mass, double crosssection) {
    return 1. / (sqrt(4 * K_CONST_K * T / (K_CONST_PI * coll_mass)) * n * crosssection);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
  // Analytical expressions for the integrals Ωcd (l,r) can be obtained only for the rigid sphere model
  // Returns the \Omega^{(l,r)}-integral for a rigid sphere potential for any l > 0 and r > 0
  double kappa::Approximation::omega_integral_RS(double T, int l, int r, double diameter, double coll_mass) {

    // pp. 145
    //std::cout << T << " " << l << " " << r << " " << diameter << " "  << coll_mass << std::endl;
    //std::cout << sqrt(T * K_CONST_K / (2 * K_CONST_PI * coll_mass)) * 0.5 * kappa::factorial_table[r + 1] 
    //             * (1.0 - 0.5 * (1.0 + pow(-1.0, l)) / (l + 1)) * K_CONST_PI * (diameter * diameter) << std::endl;

    return sqrt(T * K_CONST_K / (2 * K_CONST_PI * coll_mass)) * 0.5 * kappa::factorial_table[r + 1] 
                  * (1.0 - 0.5 * (1.0 + pow(-1.0, l)) / (l + 1)) * K_CONST_PI * (diameter * diameter);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Returns the \Omega^{(l,r)}-integral for the VSS potential
  double kappa::Approximation::omega_integral_VSS(double T, int l, int r, double coll_mass, double vss_c_cs, double vss_omega, double vss_alpha) {
    
    double mult = 1;
    switch (l) {
      case 1:
        mult = 1 / (1 + vss_alpha);
        break;
      case 2:
        mult = 2 * vss_alpha / ((2 + vss_alpha) * (1 + vss_alpha));
        break;
      case 3:
        mult = (3 * vss_alpha * vss_alpha + vss_alpha + 2) / ((3 + vss_alpha) * (2 + vss_alpha) * (1 + vss_alpha));
        break;
      case 4:
        mult = (4 * vss_alpha * vss_alpha * vss_alpha + 12 * vss_alpha * vss_alpha + 32 * vss_alpha) 
                            / ((4 + vss_alpha) * (3 + vss_alpha) * (2 + vss_alpha) * (1 + vss_alpha));
        break;
      case 5:
        mult = (5 * vss_alpha * vss_alpha * vss_alpha * vss_alpha + 30 * vss_alpha * vss_alpha * vss_alpha + 115 * vss_alpha * vss_alpha + 90 * vss_alpha + 120) / 
              ((5 + vss_alpha) * (4 + vss_alpha) * (3 + vss_alpha) * (2 + vss_alpha) * (1 + vss_alpha));
        break;
      default:
        return -1;
        break; 
    }
  
    return mult * sqrt(T * K_CONST_K / (2 * K_CONST_PI * coll_mass)) * tgamma(2.5 + r - vss_omega) * vss_c_cs * pow(2 * T * K_CONST_K / coll_mass, 0.5 - vss_omega);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Returns the \Omega^{(l,r)}-integral for the Lennard-Jones potential
  double kappa::Approximation::omega_integral_LennardJones(double T, int l, int r, double epsilon) {
	
    if (r==1 && l==1) {
      double log_KT_eps = log(K_CONST_K * T / epsilon) + 1.4;
      // Kustova-Nagnibeda eqn. 5.99 and tab. 5.5
      return 1. / ( - 0.16845 						- 
   	   	      0.02258 / (log_KT_eps * log_KT_eps) 		+ 
                      0.19779 / log_KT_eps 				+ 
     		      0.64373 * log_KT_eps 				- 
		      0.09267 * log_KT_eps * log_KT_eps 		+
		      0.00711 * log_KT_eps * log_KT_eps * log_KT_eps )	;
    } else if (r==2 && l==2) {
      double log_KT_eps = log(K_CONST_K * T / epsilon) + 1.5;
      return 1. / ( - 0.40811 						- 
		       0.05086 / (log_KT_eps * log_KT_eps) 		+ 
                       0.34010 / log_KT_eps 				+ 
                       0.70375 * log_KT_eps 				- 
                       0.10699 * log_KT_eps * log_KT_eps 		+ 
                       0.00763 * log_KT_eps * log_KT_eps * log_KT_eps )	;	
    } else {
      return -1.;
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Coefficients for Born-Mayer potential
  double kappa::Approximation::omega_integral_Born_Mayer(double T, int l, int r, double beta, double phi_zero, double diameter, double epsilon) {

    if ((r==1 && l==1) || (r==2 && l==2)) {
      double Born_Mayer_A[6] = { 0 }; 
      double log_v_star = log(phi_zero / (10 * epsilon));
      double KT_eps = K_CONST_K * T / epsilon;
        
      for (int i=0; i<6; i++) {
	  
        // coefficients for Born-Mayer potential, Kustova-Nagnibeda (5.100)
        Born_Mayer_A[i] = Born_Mayer_coeff_array[i][0] + pow(log_v_star  / (beta * diameter), -2) * 
                         (Born_Mayer_coeff_array[i][1] + 
                          Born_Mayer_coeff_array[i][2] / log_v_star +	
                          Born_Mayer_coeff_array[i][3] / (log_v_star * log_v_star) );

      }

      if (r==1 && l==1) { 
	return pow(log(phi_zero / (K_CONST_K * T)) / (beta * diameter), 2) 	* 
               ( 0.89 + Born_Mayer_A[0] / (KT_eps * KT_eps) 			+ 
                        Born_Mayer_A[1] / pow(KT_eps, 4) 			+ 
                        Born_Mayer_A[2] / pow(KT_eps, 6) )			;
      } else { // r == 2 and l == 2 
        KT_eps = log(KT_eps);
	return pow(log(phi_zero / (K_CONST_K * T)) / (beta * diameter), 2) 	* 
                       ( 1.04 + Born_Mayer_A[3] / (KT_eps * KT_eps) 		+ 
    				Born_Mayer_A[4] / (KT_eps * KT_eps * KT_eps) 	+ 
 				Born_Mayer_A[5] / pow(KT_eps, 4));
      }
    } else {
      return -1.;
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::omega_integral_ESA_hydrogen(double T, double diameter, double a1, double a2, double a3, double a4, double a5, double a6, double a7) {

    double x = log(T);
    double exp_1 = exp((x - a3) / a4);
    double exp_2 =  exp((x - a6) / a7);
    return ((a1 + a2 * x) * exp_1 / (exp_1 + 1 / exp_1) + a5 * exp_2 / (exp_2 + 1 / exp_2)) / (diameter * diameter * 1e20);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::omega_integral_ESA_H2H2(double T, int l, int r, double diameter) {
 
    if (l == 1 && r == 1) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 24.00841090, -1.61027553, 3.88885724, -8.89043396, 0.44260972, 8.88408687, -1.05402226);
    } else if (l == 1 && r == 2) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 23.02146328, -1.70509850, 3.88885724, -10.46929121, 0.36330166, 8.26405726, -1.02331842);
    } else if (l == 1 && r == 3) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 21.17218602, -1.57714612, 3.88885724, -9.72209606, 0.59112956, 8.15580488, -1.46063978);
    } else if (l == 1 && r == 4) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 20.05416161, -1.51326919, 3.88885724, -9.38278743, 0.70004430, 8.00952510, -1.57063623);
    } else if (l == 1 && r == 5) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 19.06639058, -1.45577823, 3.88885724, -9.14716131, 0.81250655, 7.85268967, -1.66995743);
    } else if (l == 2 && r == 2) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 27.54387526, -1.98253166, 3.88885724, -12.91940775, 0.34707960, 8.72131306, -0.88296275);
    } else if (l == 2 && r == 3) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 26.22527642, -1.94538819, 3.88885724, -13.40557645, 0.40398208, 8.42662474, -0.96878644);
    } else if (l == 2 && r == 4) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 24.59185702, -1.83729737, 3.88885724, -12.78050876, 0.62739891, 8.27557505, -1.33071440);
    } else if (l == 3 && r == 3) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 24.57128293, -1.80855250, 3.88885724, -11.86035430, 0.36590658, 8.38682707, -1.00746362);
    } else {
      return -1;
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::omega_integral_ESA_H2H(double T, int l, int r, double diameter) {
 
    if (l == 1 && r == 1) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 12.49063970, -1.14704753, 8.76515503, -3.52921496, 0.32874932, 12.77040465, -3.04802967);
    } else if (l == 1 && r == 2) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 12.02124035, -1.19514025, 8.76515503, -3.45192920, 0.45922882, 12.77040465, -2.29080329);
    } else if (l == 1 && r == 3) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 11.69204285, -1.24240232, 8.76515503, -3.49608019, 0.63354264, 12.77040465, -2.29080329);
    } else if (l == 1 && r == 4) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 11.45792771, -1.29677120, 8.76515503, -3.64478512, 0.85298582, 12.77040465, -2.29080329);
    } else if (l == 1 && r == 5) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 11.00483923, -1.27212994, 8.76515503, -3.51537463, 0.85298582, 12.77040465, -2.29080329);
    } else if (l == 2 && r == 2) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 7.45048892, -1.43326160, 9.59201391, -1.35066206, 7.15859874, 9.88881724, -1.39484886);
    } else if (l == 2 && r == 3) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 10.84507417, -1.42859529, 9.20889644, -1.29890434, 3.37747184, 9.83307970, -1.30321649);
    } else if (l == 2 && r == 4) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 11.55088396, -1.41480945, 8.98739895, -1.39880703, 2.32276221, 9.89142509, -1.26804718);
    } else if (l == 3 && r == 3) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, -15.25288758, -1.39293852, 9.59147724, -1.62599901, 28.71128123, 9.68396961, -1.63186985);
    } else {
      return -1;
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::omega_integral_ESA_HH(double T, int l, int r, double diameter) {
 
    if (l == 1 && r == 1) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 15.09506044, -1.25710008, 9.57839369, -3.80371463, 0.98646613, 9.25705877, -0.93611707);
    } else if (l == 1 && r == 2) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 14.14566908, -1.17057105, 9.02830724, -3.00779776, 0.74653903, 9.10299040, -0.68184353);
    } else if (l == 1 && r == 3) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 13.39722075, -1.09886403, 8.50097335, -2.86025395, 0.85345727, 8.90666490, -0.67571329);
    } else if (l == 1 && r == 4) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 12.97073246, -1.06479185, 8.18885522, -2.78105132, 0.89401865, 8.73403138, -0.65658782);
    } else if (l == 1 && r == 5) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 12.69248000, -1.04857945, 7.97861283, -2.73621289, 0.90816787, 8.57840253, -0.63732002);
    } else if (l == 2 && r == 2) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 22.08948804, -1.85066626, 8.50932055, -7.66943974, 0.77454531, 9.69545318, -0.62104466);
    } else if (l == 2 && r == 3) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 17.94703897, -1.42488999, 7.66669340, -4.76239721, 1.26783524, 9.53716768, -0.73914215);
    } else if (l == 2 && r == 4) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 18.78590499, -1.59291967, 7.97734302, -5.66814860, 1.01816360, 9.32328437, -0.60882006);
    } else if (l == 3 && r == 3) {
      return kappa::Approximation::omega_integral_ESA_hydrogen(T, diameter, 13.82986524, -1.01454290, 7.48970759, -3.27628187, 2.08225623, 9.21388055, -1.32086596);
    } else {
      return -1;
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::omega_integral_ESA_nn(double T, int l, int r, double beta, double epsilon_zero, double r_e) {

    double x = log(K_CONST_K * T / epsilon_zero);
    if (l == 1 && r == 1) {
      double a[7] = { 	(7.884756e-1) - (2.438494e-2) * beta, 
			-2.952759e-1 - (1.744149e-3) * beta, 
			(5.020892e-1) + (4.316985e-2) * beta, 
			(-9.042460e-1) - (4.017103e-2) * beta,
			(-3.373058) + (2.458538e-1) * beta - 
			(4.850047e-3) * beta * beta, 
			(4.161981) + (2.202737e-1) * beta - (1.718010e-2) * beta * beta, 
			(2.462523) + (3.231308e-1) * beta - (2.281072e-2) * beta * beta };

      // Optimization: vectorized product and store coefficients seperately?
      double exp_1 = exp((x - a[2]) / a[3]);
      double exp_2 =  exp((x - a[5]) / a[6]);
      return exp((a[0] + a[1] * x) * exp_1 / (exp_1 + 1 / exp_1) + a[4] * exp_2 / (exp_2 + 1 / exp_2));
    } else if (l == 1 && r == 2) {
      double a[7] = { (7.123565e-1) + (-2.688875e-2) * beta, 
                      (-2.9105e-1) + (-2.065175e-3) * beta, 
                      (4.187065e-2) + (4.060236e-2) * beta, 
                      (-9.287685e-1) + (-2.342270e-2) * beta,
		      (-3.598542) + (2.54512e-1) * beta + 
                      (-4.685966e-3) * beta * beta, 
                      (3.934824) + (2.699944e-1) * beta + (-2.009886e-2) * beta * beta, 
                       2.578084 + (3.449024e-1) * beta + (-2.2927e-2) * beta * beta };

      double exp_1 = exp((x - a[2]) / a[3]);
      double exp_2 =  exp((x - a[5]) / a[6]);
      return exp((a[0] + a[1] * x) * exp_1 / (exp_1 + 1 / exp_1) + a[4] * exp_2 / (exp_2 + 1 / exp_2));
    } else if (l == 1 && r == 3) {
      double a[7] = { (6.606022e-1) + (-2.831448e-2) * beta, 
 		      (-2.8709e-1) + (-2.232827e-3) * beta, 
                      (-2.51969e-1) + (3.778211e-2) * beta, 
                      (-9.173046e-1) + (-1.864476e-2) * beta,
		      (-3.776812) + (2.552528e-1) * beta + (-4.23722e-3) * beta * beta,
                      (3.768103) + (3.155025e-1) * beta + (-2.218849e-2) * beta * beta, 
                       2.695440 + (3.597998e-1) * beta + (-2.267102e-2) * beta * beta };

      double exp_1 = exp((x - a[2]) / a[3]);
      double exp_2 =  exp((x - a[5]) / a[6]);
      return exp((a[0] + a[1] * x) * exp_1 / (exp_1 + 1 / exp_1) + a[4] * exp_2 / (exp_2 + 1 / exp_2));
    } else if (l == 1 && r == 4) {
      double a[7] = { (6.268016e-1) + (-2.945078e-2) * beta, 
                      (-2.830834e-1) + (-2.361273e-3) * beta, 
                      (-4.559927e-1) + (3.705640e-2) * beta, 
                      (-9.334638e-1) + (-1.797329e-2) * beta,
		      (-3.947019) + (2.446843e-1) * beta + (-3.176374e-3) * beta * beta, 
                      (3.629926) + (3.761272e-1) * beta + (-2.451016e-2) * beta * beta, 
                       2.824905 + (3.781709e-1) * beta + (-2.251978e-2) * beta * beta };

      double exp_1 = exp((x - a[2]) / a[3]);
      double exp_2 =  exp((x - a[5]) / a[6]);
      return exp((a[0] + a[1] * x) * exp_1 / (exp_1 + 1 / exp_1) + a[4] * exp_2 / (exp_2 + 1 / exp_2));
    } else if (l == 1 && r == 5) {
      double a[7] = { (5.956859e-1) + (-2.915893e-2) * beta, 
                      (-2.804989e-1) + (-2.298968e-3) * beta, 
                      (-5.965551e-1) + (3.724395e-2) * beta, 
                      (-8.946001e-1) + (-2.550731e-2) * beta,
		      (-4.076798) + (1.983892e-1) * beta + (-5.014065e-3) * beta * beta, 
                      (3.458362) + (4.770695e-1) * beta + (-2.678054e-2) * beta * beta, 
                       2.982260 + (4.014572e-1) * beta + (-2.142580e-2) * beta * beta };
      double exp_1 = exp((x - a[2]) / a[3]);
      double exp_2 =  exp((x - a[5]) / a[6]);
      return exp((a[0] + a[1] * x) * exp_1 / (exp_1 + 1 / exp_1) + a[4] * exp_2 / (exp_2 + 1 / exp_2));
    } else if (l == 2 && r == 2) {
      double a[7] = { (7.898524e-1) - (2.114115e-2) * beta, 
                      (-2.998325e-1) - (1.243977e-3) * beta, 
                      (7.077103e-1) + (3.583907e-2) * beta, 
                      (-8.946857e-1) - (2.473947e-2) * beta,
		      (-2.958969) + (2.303358e-1) * beta - (5.226562e-3) * beta * beta, 
                      (4.348412) + (1.920321e-1) * beta - (1.496557e-2) * beta * beta, 
                      (2.205440) + (2.567027e-1) * beta - (1.861359e-2) * beta * beta };
      double exp_1 = exp((x - a[2]) / a[3]);
      double exp_2 =  exp((x - a[5]) / a[6]);
      return exp((a[0] + a[1] * x) * exp_1 / (exp_1 + 1 / exp_1) + a[4] * exp_2 / (exp_2 + 1 / exp_2));
    } else if (l == 2 && r == 3) {
      double a[7] = { (7.269006e-1) - (2.233866e-2) * beta, 
                      (-2.972304e-1) - (1.392888e-3) * beta, 
                      (3.904230e-1) + (3.231655e-2) * beta, (-9.442201e-1) - (1.494805e-2) * beta,
		      (-3.137828) + (2.347767e-1) * beta - (4.963979e-3) * beta * beta, 
                      (4.190370) + (2.346004e-1) * beta - (1.718963e-2) * beta * beta, 
                      (2.319751) + (2.700236e-1) * beta - (1.854217e-2) * beta * beta };
      double exp_1 = exp((x - a[2]) / a[3]);
      double exp_2 =  exp((x - a[5]) / a[6]);
      return exp((a[0] + a[1] * x) * exp_1 / (exp_1 + 1 / exp_1) + a[4] * exp_2 / (exp_2 + 1 / exp_2));
    } else if (l == 2 && r == 4) {
      double a[7] = { (6.829159e-1) - (2.233866e-2) * beta, 
                      (-2.943232e-1) - (1.514322e-3) * beta, 
                      (1.414623e-1) + (3.075351e-2) * beta, 
                      (-9.720228e-1) - (1.038869e-2) * beta,
	  	      (-3.284219) + (2.243767e-1) * beta - (3.913041e-3) * beta * beta, 
                      (4.011692) + (3.005083e-1) * beta - (2.012373e-2) * beta * beta, 
                      (2.401249) + (2.943600e-1) * beta - (1.884503e-2) * beta * beta };
      double exp_1 = exp((x - a[2]) / a[3]);
      double exp_2 =  exp((x - a[5]) / a[6]);
      return exp((a[0] + a[1] * x) * exp_1 / (exp_1 + 1 / exp_1) + a[4] * exp_2 / (exp_2 + 1 / exp_2));
    } else if (l == 3 && r == 3) {
      double a[7] = { (7.468781e-1) - (2.518134e-2) * beta, 
                      (-2.947438e-1) - (1.811571e-3) * beta,
                      (2.234096e-1) + (3.681114e-2) * beta, 
                      (-9.974591e-1) - (2.670805e-2) * beta,
		      (-3.381787) + (2.372932e-1) * beta - (4.239629e-3) * beta * beta, 
                      (4.094540) + (2.756466e-1) * beta - (2.009227e-2) * beta * beta, 
                      (2.476087) + (3.300898e-1) * beta - (2.223317e-2) * beta * beta };
      double exp_1 = exp((x - a[2]) / a[3]);
      double exp_2 =  exp((x - a[5]) / a[6]);
      return exp((a[0] + a[1] * x) * exp_1 / (exp_1 + 1 / exp_1) + a[4] * exp_2 / (exp_2 + 1 / exp_2));
    } else if (l == 4 && r == 4) {
      double a[7] = { (7.365470e-1) - (2.242357e-2) * beta, 
                      (-2.968650e-1) - (1.396696e-3) * beta, 
                      (3.747555e-1) + (2.847063e-2) * beta, 
                      (-9.944036e-1) - (1.378926e-2) * beta,
		      (-3.136655) + (2.176409e-1) * beta - (3.899247e-3) * beta * beta,
                      (4.145871) + (2.855836e-1) * beta - (1.939452e-2) * beta * beta, 
                      (2.315532) + (2.842981e-1) * beta - (1.874462e-2) * beta * beta };
      double exp_1 = exp((x - a[2]) / a[3]);
      double exp_2 =  exp((x - a[5]) / a[6]);
      return exp((a[0] + a[1] * x) * exp_1 / (exp_1 + 1 / exp_1) + a[4] * exp_2 / (exp_2 + 1 / exp_2));
    } else {
      return -1;
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::omega_integral_ESA_cn(double T, int l, int r, double beta, double epsilon_zero, double r_e) {

    // double m = 4; // parameter for neutral-ion interactions ESA STR 256 Table 2 [dimensionless]
    double x = log(K_CONST_K * T / epsilon_zero);
    if (l == 1 && r == 1) {
      double a[7] = { 	(9.851755e-1) - (2.870704e-2) * beta, 
			-(4.737800e-1) - (1.370344e-3) * beta, 
			(7.080799e-1) + (4.575312e-3) * beta, 
			(-1.239439) - (3.683605e-2) * beta,
			(-4.638467) + (4.418904e-1) * beta - (1.220292e-2) * beta * beta, 
			(3.841835) + (3.277658e-1) * beta - (2.660275e-2) * beta * beta, 
			(2.317342) + (3.912768e-1) * beta - (3.136223e-2) * beta * beta} ;

      double exp_1 = exp((x - a[2]) / a[3]);
      double exp_2 =  exp((x - a[5]) / a[6]);
      return exp((a[0] + a[1] * x) * exp_1 / (exp_1 + 1 / exp_1) + a[4] * exp_2 / (exp_2 + 1 / exp_2));
    } else if (l == 1 && r == 2) {
      double a[7] = { (8.361751e-1) - (3.201292e-2) * beta, 
		      -(4.707355e-1) - (1.783284e-3) * beta, 
		      (1.771157e-1) + (1.172773e-2) * beta, 
		      (-1.094937) - (3.115598e-2) * beta,
	              (-4.976384) + (4.708074e-1) * beta - (1.283818e-2) * beta * beta, 
 		      (3.645873) + (3.699452e-1) * beta - (2.988684e-2) * beta * beta, 
 		      (2.428864) + (4.267351e-1) * beta - (3.278874e-2) * beta * beta };

      double exp_1 = exp((x - a[2]) / a[3]);
      double exp_2 =  exp((x - a[5]) / a[6]);
      // return exp((a[0] + a[1] * x) * exp((x - a[2]) / a[3]) / (exp((x - a[2]) / a[3]) + exp((a[2] - x) / a[3])) + a[4] * exp((x - a[5]) / a[6]) / (exp((x - a[5]) / a[6]) + exp((a[5] - x) / a[6]))) * sigma_sqr;
      return exp((a[0] + a[1] * x) * exp_1 / (exp_1 + 1 / exp_1) + a[4] * exp_2 / (exp_2 + 1 / exp_2));
    } else if (l == 1 && r == 3) {
      double a[7] = { (7.440562e-1) - (3.453851e-2) * beta,
                     -(4.656306e-1) - (2.097901e-3) * beta, 
		      (-1.465752e-1) + (1.446209e-2) * beta, 
		      (-1.080410) - (2.712029e-2) * beta,
		      (-5.233907) + (4.846691e-1) * beta - (1.280346e-2) * beta * beta, 
		      (3.489814) + (4.140270e-1) * beta - (3.250138e-2) * beta * beta, 
		      (2.529678) + (4.515088e-1) * beta - (3.339293e-2) * beta * beta };

      double exp_1 = exp((x - a[2]) / a[3]);
      double exp_2 =  exp((x - a[5]) / a[6]);
      // return exp((a[0] + a[1] * x) * exp((x - a[2]) / a[3]) / (exp((x - a[2]) / a[3]) + exp((a[2] - x) / a[3])) + a[4] * exp((x - a[5]) / a[6]) / (exp((x - a[5]) / a[6]) + exp((a[5] - x) / a[6]))) * sigma_sqr;
      return exp((a[0] + a[1] * x) * exp_1 / (exp_1 + 1 / exp_1) + a[4] * exp_2 / (exp_2 + 1 / exp_2));
    } else if (l == 1 && r == 4) {
      double a[7] = { (6.684360e-1) - (3.515695e-2) * beta,
                     -(4.622014e-1) - (2.135808e-3) * beta, 
		      (-3.464990e-1) + (1.336362e-2) * beta, 
		      (-1.054374) - (3.149321e-2) * beta,
		      (-5.465789) + (4.888443e-1) * beta - (1.228090e-2) * beta * beta, 
		      (3.374614) + (4.602468e-1) * beta - (3.463073e-2) * beta * beta, 
		      (2.648622) + (4.677409e-1) * beta - (3.339297e-2) * beta * beta };

      double exp_1 = exp((x - a[2]) / a[3]);
      double exp_2 =  exp((x - a[5]) / a[6]);
      // return exp((a[0] + a[1] * x) * exp((x - a[2]) / a[3]) / (exp((x - a[2]) / a[3]) + exp((a[2] - x) / a[3])) + a[4] * exp((x - a[5]) / a[6]) / (exp((x - a[5]) / a[6]) + exp((a[5] - x) / a[6]))) * sigma_sqr;
      return exp((a[0] + a[1] * x) * exp_1 / (exp_1 + 1 / exp_1) + a[4] * exp_2 / (exp_2 + 1 / exp_2));
    } else if (l == 1 && r == 5) {
      double a[7] = { (6.299083e-1) - (3.720000e-2) * beta,
		     -(4.560921e-1) - (2.395779e-3) * beta,
		      (-5.228598e-1) + (1.594610e-2) * beta,
		      (-1.124725) - (2.862354e-2) * beta,
		      (-5.687354) + (4.714789e-1) * beta - (1.056602e-2) * beta * beta,
		      (3.267709) + (5.281419e-1) * beta - (3.678869e-2) * beta * beta,
		      (2.784725) + (4.840700e-1) * beta - (3.265127e-2) * beta * beta };

      double exp_1 = exp((x - a[2]) / a[3]);
      double exp_2 =  exp((x - a[5]) / a[6]);
      // return exp((a[0] + a[1] * x) * exp((x - a[2]) / a[3]) / (exp((x - a[2]) / a[3]) + exp((a[2] - x) / a[3])) + a[4] * exp((x - a[5]) / a[6]) / (exp((x - a[5]) / a[6]) + exp((a[5] - x) / a[6]))) * sigma_sqr;
      return exp((a[0] + a[1] * x) * exp_1 / (exp_1 + 1 / exp_1) + a[4] * exp_2 / (exp_2 + 1 / exp_2));
    } else if (l == 2 && r == 2) {
      double a[7] = { (9.124518e-1) + (-2.398461e-2) * beta,
    	              (-4.697184e-1) + (-7.809681e-4) * beta,
		      (1.031053) + (4.069668e-3) * beta,
		      (-1.090782) + (-2.413508e-2) * beta,
		      (-4.127243) + (4.302667e-1) * beta + (-1.352874e-2) * beta * beta,
		      (4.059078) + (2.597379e-1) * beta + (-2.169951e-2) * beta * beta,
		      (2.086906) + (2.920310e-1) * beta + (-2.560437e-2) * beta * beta };

      double exp_1 = exp((x - a[2]) / a[3]);
      double exp_2 =  exp((x - a[5]) / a[6]);
      // return exp((a[0] + a[1] * x) * exp((x - a[2]) / a[3]) / (exp((x - a[2]) / a[3]) + exp((a[2] - x) / a[3])) + a[4] * exp((x - a[5]) / a[6]) / (exp((x - a[5]) / a[6]) + exp((a[5] - x) / a[6]))) * sigma_sqr;
      return exp((a[0] + a[1] * x) * exp_1 / (exp_1 + 1 / exp_1) + a[4] * exp_2 / (exp_2 + 1 / exp_2));
    } else if (l == 2 && r == 3) {
      double a[7] = { (8.073459e-1) + (-2.581232e-2) * beta,
  		      (-4.663682e-1) + (-1.030271e-3) * beta,
		      (6.256342e-1) + (4.086881e-3) * beta,
		      (-1.063437) + (-1.235489e-2) * beta,
		      (-4.365989) + (4.391454e-1) * beta + (-1.314615e-2) * beta * beta,
		      (3.854346) + (3.219224e-1) * beta + (-2.587493e-2) * beta * beta,
		      (2.146207) + (3.325620e-1) * beta + (-2.686959e-2) * beta * beta };

      double exp_1 = exp((x - a[2]) / a[3]);
      double exp_2 =  exp((x - a[5]) / a[6]);
      // return exp((a[0] + a[1] * x) * exp((x - a[2]) / a[3]) / (exp((x - a[2]) / a[3]) + exp((a[2] - x) / a[3])) + a[4] * exp((x - a[5]) / a[6]) / (exp((x - a[5]) / a[6]) + exp((a[5] - x) / a[6]))) * sigma_sqr;
      return exp((a[0] + a[1] * x) * exp_1 / (exp_1 + 1 / exp_1) + a[4] * exp_2 / (exp_2 + 1 / exp_2));
    } else if (l == 2 && r == 4) {
      double a[7] = { (7.324117e-1) + (-2.727580e-2) * beta,
		      (-4.625614e-1) + (-1.224292e-3) * beta,
		      (3.315871e-1) + (7.216776e-3) * beta,
		      (-1.055706) + (-8.585500e-3) * beta,
		      (-4.571022) + (4.373660e-1) * beta + (-1.221457e-2) * beta * beta,
		      (3.686006) + (3.854493e-1) * beta + (-2.937568e-2) * beta * beta,
		      (2.217893) + (3.641196e-1) * beta + (-2.763824e-2) * beta * beta };

      double exp_1 = exp((x - a[2]) / a[3]);
      double exp_2 =  exp((x - a[5]) / a[6]);
      // return exp((a[0] + a[1] * x) * exp((x - a[2]) / a[3]) / (exp((x - a[2]) / a[3]) + exp((a[2] - x) / a[3])) + a[4] * exp((x - a[5]) / a[6]) / (exp((x - a[5]) / a[6]) + exp((a[5] - x) / a[6]))) * sigma_sqr;
      return exp((a[0] + a[1] * x) * exp_1 / (exp_1 + 1 / exp_1) + a[4] * exp_2 / (exp_2 + 1 / exp_2));
    } else if (l == 3 && r == 3) {
      double a[7] = { (8.402943e-1) + (-2.851694e-2) * beta,
		      (-4.727437e-1) + (-1.328784e-3) * beta,
		      (4.724228e-1) + (7.706027e-3) * beta,
		      (-1.213660) + (-3.456656e-2) * beta,
		      (-4.655574) + (4.467685e-1) * beta + (-1.237864e-2) * beta * beta,
		      (3.817178) + (3.503180e-1) * beta + (-2.806506e-2) * beta * beta,
		      (2.313186) + (3.889828e-1) * beta + (-3.120619e-2) * beta * beta };

      double exp_1 = exp((x - a[2]) / a[3]);
      double exp_2 =  exp((x - a[5]) / a[6]);
      // return exp((a[0] + a[1] * x) * exp((x - a[2]) / a[3]) / (exp((x - a[2]) / a[3]) + exp((a[2] - x) / a[3])) + a[4] * exp((x - a[5]) / a[6]) / (exp((x - a[5]) / a[6]) + exp((a[5] - x) / a[6]))) * sigma_sqr;
      return exp((a[0] + a[1] * x) * exp_1 / (exp_1 + 1 / exp_1) + a[4] * exp_2 / (exp_2 + 1 / exp_2));
    } else if (l == 4 && r == 4) {
      double a[7] = { (8.088842e-1) + (-2.592379e-2) * beta,
		      (-4.659483e-1) + (-1.041599e-3) * beta,
		      (6.092981e-1) + (1.428402e-3) * beta,
		      (-1.113323) + (-1.031574e-2) * beta,
		      (-4.349145) + (4.236246e-1) * beta + (-1.210668e-2) * beta * beta,
		      (3.828467) + (3.573461e-1) * beta + (-2.759622e-2) * beta * beta,
		      (2.138075) + (3.388072e-1) * beta + (-2.669344e-2) * beta * beta };

      double exp_1 = exp((x - a[2]) / a[3]);
      double exp_2 =  exp((x - a[5]) / a[6]);
      // return exp((a[0] + a[1] * x) * exp((x - a[2]) / a[3]) / (exp((x - a[2]) / a[3]) + exp((a[2] - x) / a[3])) + a[4] * exp((x - a[5]) / a[6]) / (exp((x - a[5]) / a[6]) + exp((a[5] - x) / a[6]))) * sigma_sqr;
      return exp((a[0] + a[1] * x) * exp_1 / (exp_1 + 1 / exp_1) + a[4] * exp_2 / (exp_2 + 1 / exp_2));
    } else {
      return -1;
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::omega_integral_ESA_cn_corr(double T, double d1, double d2, double d3, double beta, double r_e) {

    double x = log(T);
    double x0 = 0.7564 * pow(beta, 0.064605); //  ESA STR 256 (5.6); Table 4; x0 = xi_1 * beta^xi_2
    double sigma_sqr = x0 * x0 * r_e * r_e / 1e20; // value to obtain dimensional \Omega^{(l,r)}-integral [m]; 

    return (d1 + d2*x + d3*x*x) / (sigma_sqr);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::omega_integral_ESA_cc(double T, int l, int r, int charge1, int charge2, double debye_length) {

    double nl = 1. / (1 - (1 + pow(-.1, l)) / (2 * l + 2));
    double Tstar = fabs(debye_length * 4 * K_CONST_PI * K_CONST_E0 * K_CONST_K * T / (K_CONST_ELEMENTARY_CHARGE * K_CONST_ELEMENTARY_CHARGE * charge1 * charge2));

    return l * nl * log((4 * Tstar / (K_CONST_EULER * K_CONST_EULER)) * exp(cc_ar[r-1] - cc_cl[l-1]) + 1) / (Tstar * Tstar * r * (r + 1));
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::omega_integral_ESA_ne(double T, double g1, double g2, double g3, double g4, double g5, double g6, double g7, double g8,
                                                     double g9, double g10, double diameter) {

    double x = log(T);
    double exp_1 = exp(x - g1 / g2);

    return (g3 * pow(x, g6) * exp_1 / (exp_1 + 1. / exp_1) + g7 * exp(-(x - g8) * (x - g8) / (g9 * g9)) + g4 + g10 * pow(x, g5)) / (diameter * diameter);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // compute the backward and forward reaction rate coefficients for VT exchange (private)
  double kappa::Approximation::p_k_bf_VT(double T, 
  				         double vibr_energy_before, 
					 double vibr_energy_after, 
					 const arma::vec &rot_energy_before, 
					 int num_rot_levels_before, 
                                     	 const arma::vec &rot_energy_after, 
					 int num_rot_levels_after, 
					 int rot_symmetry, 
					 bool is_rigid_rotator) {

    if (is_rigid_rotator) {
      return exp((vibr_energy_after - vibr_energy_before) / (K_CONST_K * T));
    } else {
      // simplified for VT eq. 2.74 in book Kustova & Nagnibeda, 2003/2009.
      return exp((vibr_energy_after - vibr_energy_before) / (K_CONST_K * T)) * 
             p_Z_rot(T, rot_energy_before, num_rot_levels_before, rot_symmetry) /  
             p_Z_rot(T, rot_energy_after, num_rot_levels_after, rot_symmetry);
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // compute the backward and forward reaction rate coefficients for VV exchange (private)
  double kappa::Approximation::p_k_bf_VV(double T, 
  					 double vibr_energy1_before, 
					 double vibr_energy1_after, 
					 double vibr_energy2_before, 
					 double vibr_energy2_after,
                                     	 const arma::vec &rot_energy1_before, 
					 int num_rot_levels1_before, 
					 const arma::vec &rot_energy1_after, 
					 int num_rot_levels1_after, 
					 int rot_symmetry1, 
                                     	 const arma::vec &rot_energy2_before, 
					 int num_rot_levels2_before, 
					 const arma::vec &rot_energy2_after, 
					 int num_rot_levels2_after, 
					 int rot_symmetry2) { 

    // eq. 2.74 in book Kustova & Nagnibeda, 2003/2009.
    return exp((vibr_energy1_after - vibr_energy1_before + vibr_energy2_after - vibr_energy2_before) / (K_CONST_K * T)) *
           (p_Z_rot(T, rot_energy1_before, num_rot_levels1_before, rot_symmetry1) * 
 	   p_Z_rot(T, rot_energy2_before, num_rot_levels2_before, rot_symmetry2)) /
           (p_Z_rot(T, rot_energy1_after,  num_rot_levels1_after,  rot_symmetry1) * 
	   p_Z_rot(T, rot_energy2_after,  num_rot_levels2_after,  rot_symmetry2));
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the backward and forward reaction rate coefficients for chemical exchange (private)
  double kappa::Approximation::p_k_bf_exch(	double T, 
						double molecule_before_mass, 
						double atom_before_mass, 
						double molecule_after_mass, 
						double atom_after_mass,
                                       		double diss_energy_before, 
						double diss_energy_after, 
						double vibr_energy_before, 
						double vibr_energy_after,
                                       		const arma::vec &rot_energy_before, 
						int num_rot_levels_before, 
						int rot_symmetry_before,
                                      		const arma::vec &rot_energy_after, 
						int num_rot_levels_after, 
						int rot_symmetry_after) { 

    // look at eq. 36 in system.pdf by Olga Kunova or eq. 2.78 in book Kustova & Nagnibeda, 2003/2009.
    return pow(molecule_before_mass * atom_before_mass / (molecule_after_mass * atom_after_mass), 1.5) * 
           (p_Z_rot(T, rot_energy_before, num_rot_levels_before, rot_symmetry_before) / 
     	   p_Z_rot(T, rot_energy_after, num_rot_levels_after, rot_symmetry_after)) * 
           exp((vibr_energy_after - vibr_energy_before) / (K_CONST_K * T)) * 
           exp((diss_energy_before - diss_energy_after) / (K_CONST_K * T));
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // compute the backward and forward reaction rate coefficients for dissociation-recombination (private)
  double kappa::Approximation::p_k_bf_diss(	double T, 
						double molecule_mass, 
						double atom1_mass, 
						double atom2_mass, 
						double diss_energy, 
						double vibr_energy,
                                       		const arma::vec &rot_energy, 
						int num_rot_levels, 
						int rot_symmetry) { 

    // look at eq. 37 in system.pdf by Olga Kunova or eq. 2.79 in book Kustova & Nagnibeda, 2003/2009.
    return pow(molecule_mass / (atom1_mass * atom2_mass), 1.5) * 
 	   K_CONST_H * K_CONST_H * K_CONST_H * pow(2 * K_CONST_PI * K_CONST_K * T, -1.5) * 
	   p_Z_rot(T, rot_energy, num_rot_levels, rot_symmetry) *
	   exp((diss_energy - vibr_energy) / (K_CONST_K * T));
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::debye_length(double T, const arma::vec &concentrations, const arma::vec &charges) {

    #ifdef KAPPA_STRICT_CHECKS
    if (T <= 0) {
      std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif
  
    return sqrt(K_CONST_E0 * K_CONST_K * T / (arma::dot(concentrations, charges % charges)));
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int kappa::Approximation::max_electron_level(kappa::Atom const &atom, double Delta_E) {

    #ifdef KAPPA_STRICT_CHECKS
    if (Delta_E <= 0) {
      std::string error_string = "Non-positive value of Delta_E specified: Delta_E=" + std::to_string(Delta_E);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    return p_max_electron_level(atom.electron_energy, atom.ionization_potential, Delta_E);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::Z_rot(double T, const kappa::Molecule &molecule, int i, int e) {

    #ifdef KAPPA_STRICT_CHECKS
    if (T <= 0) {
      std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (i < 0) {
      std::string error_string = "Negative vibrational level: i=" + std::to_string(i);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e < 0) {
      std::string error_string = "Negative electronic level: e=" + std::to_string(e);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e > molecule.num_electron_levels-1) {
      std::string error_string = "Electron level larger than max. level: e=" + std::to_string(e) + ", max. level=" + std::to_string(molecule.num_electron_levels-1);
      throw kappa::IncorrectValueException(error_string.c_str());
    }	
    if (i > molecule.num_vibr_levels[e]-1) {
      std::string error_string = "Vibrational level larger than max. level: i=" + std::to_string(i) + ", max. level=" + std::to_string(molecule.num_vibr_levels[e]-1);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif
  
    return p_Z_rot(T, molecule.rot_energy[e][i], molecule.num_rot_levels[e][i], molecule.rot_symmetry);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::Z_vibr_eq(double T, const kappa::Molecule &molecule, int e) {

    #ifdef KAPPA_STRICT_CHECKS
    if (T <= 0) {
      std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e < 0) {
      std::string error_string = "Negative electronic level: e=" + std::to_string(e);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e > molecule.num_electron_levels-1) {
      std::string error_string = "Electron level larger than max. level: e=" + std::to_string(e) + ", max. level=" + std::to_string(molecule.num_electron_levels-1);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    return p_Z_vibr_eq(T, molecule.vibr_energy[e]);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::Z_electron(double T, const kappa::Molecule &molecule, int num_electron_levels) {

    #ifdef KAPPA_STRICT_CHECKS
    if (T <= 0) {
      std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if ((num_electron_levels < -1) || (num_electron_levels == 0)) {
      std::string error_string = "Negative or 0 value of num. of electronic levels (can only be -1 or >0): num_electron_levels=" + std::to_string(num_electron_levels);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (num_electron_levels > molecule.num_electron_levels) {
      std::string error_string = "Number of electron levels passed to function larger than number of levels: num_electron_levels:"
                             + std::to_string(num_electron_levels) + ", number of levels" + std::to_string(molecule.num_electron_levels);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif
    if (num_electron_levels == -1) {
      return p_Z_electron(T, molecule.electron_energy, molecule.statistical_weight, molecule.num_electron_levels);
    } else {
      return p_Z_electron(T, molecule.electron_energy, molecule.statistical_weight, num_electron_levels);
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::avg_rot_energy(double T, const kappa::Molecule &molecule, int i, int e) {

    #ifdef KAPPA_STRICT_CHECKS
    if (T <= 0) {
      std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (i < 0) {
      std::string error_string = "Negative vibrational level: i=" + std::to_string(i);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e < 0) {
      std::string error_string = "Negative electronic level: e=" + std::to_string(e);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e > molecule.num_electron_levels-1) {
      std::string error_string = "Electron level larger than max. level: e=" + std::to_string(e) + ", max. level=" + std::to_string(molecule.num_electron_levels-1);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (i > molecule.num_vibr_levels[e]-1) {
      std::string error_string = "Vibrational level larger than max. level: i=" + std::to_string(i) + ", max. level=" + std::to_string(molecule.num_vibr_levels[e]-1);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    return p_avg_rot_energy(T, molecule.rot_energy[e][i], molecule.num_rot_levels[e][i], molecule.rot_symmetry);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::avg_rot_energy_sq(double T, const kappa::Molecule &molecule, int i, int e) {

    #ifdef KAPPA_STRICT_CHECKS
    if (T <= 0) {
      std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (i < 0) {
      std::string error_string = "Negative vibrational level: i=" + std::to_string(i);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e < 0) {
      std::string error_string = "Negative electronic level: e=" + std::to_string(e);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e > molecule.num_electron_levels-1) {
      std::string error_string = "Electron level larger than max. level: e=" + std::to_string(e) + ", max. level=" + std::to_string(molecule.num_electron_levels-1);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (i > molecule.num_vibr_levels[e]-1) {
      std::string error_string = "Vibrational level larger than max. level: i=" + std::to_string(i) + ", max. level=" + std::to_string(molecule.num_vibr_levels[e]-1);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    return p_avg_rot_energy_sq(T, molecule.rot_energy[e][i], molecule.num_rot_levels[e][i], molecule.rot_symmetry);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::c_tr(double T, const kappa::Particle &particle) {

    #ifdef KAPPA_STRICT_CHECKS
    if (T <= n) {
      std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    return p_c_tr(T, particle.mass);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::c_rot(double T, const kappa::Molecule &molecule, int i, int e) {

    #ifdef KAPPA_STRICT_CHECKS
    if (T <= 0) {
      std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (i < 0) {
      std::string error_string = "Negative vibrational level: i=" + std::to_string(i);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e < 0) {
      std::string error_string = "Negative electronic level: e=" + std::to_string(e);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e > molecule.num_electron_levels-1) {
      std::string error_string = "Electron level > than maximum level: e=" + std::to_string(e) + ", maximum level=" + std::to_string(molecule.num_electron_levels-1);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (i > molecule.num_vibr_levels[e]-1) {
      std::string error_string = "Vibrational level > than maximum level: i=" + std::to_string(i) + ", maximum level=" + std::to_string(molecule.num_vibr_levels[e]-1);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    return p_c_rot(T, molecule.mass, molecule.rot_energy[e][i], molecule.num_rot_levels[e][i], molecule.rot_symmetry);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // LC: computation of vibrational specific heat (simplified formula)
  double kappa::Approximation::c_vibr_approx(double T, const kappa::Molecule & molecule) {

    double char_vibr_temp = molecule.characteristic_vibr_temperatures[0];
    return K_CONST_K/molecule.mass * (char_vibr_temp/T)*(char_vibr_temp/T)*exp(char_vibr_temp/T)/( (exp(char_vibr_temp/T) -1.) * (exp(char_vibr_temp/T) -1.) );
  } 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::probability_VV(	double rel_vel, 	
 						kappa::Molecule const &molecule1, 	
						kappa::Molecule const &molecule2, 
						kappa::Interaction const &interaction,
						int i, int k, int delta_i, int e1, int e2, 
						kappa::models_prob_vv model) {

    #ifdef KAPPA_STRICT_CHECKS
    if (rel_vel < 0.0) {
      std::string error_string = "Negative relative velocity: rel_vel=" + std::to_string(rel_vel);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (delta_i == 0) {
      std::string error_string = "No change in vibrational levels of molecules after VV exchange: delta_i=0";
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e1 < 0) {
      std::string error_string = "Negative electronic level of first molecule: e1=" + std::to_string(e1);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e2 < 0) {
      std::string error_string = "Negative electronic level of second molecule: e2=" + std::to_string(e2);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if ((i < 0) || (k < 0) || (i + delta_i < 0) || (k - delta_i < 0)) {
      std::string error_string = "One of the vibrational levels before or after the VV exchange is negative: i, k, i', k'=";
      error_string += std::to_string(i) + ", " + std::to_string(k) + ", " + std::to_string(i + delta_i) + ", " + std::to_string(k - delta_i);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (((interaction.particle1_name != molecule1.name) || (interaction.particle2_name != molecule2.name)) && 
        ((interaction.particle1_name != molecule2.name) || (interaction.particle2_name != molecule1.name))) {

      std::string error_string = "Interaction is for a different set of particles, the interaction is for: " + interaction.particle1_name + " + " 
                                                                                                             + interaction.particle2_name;
      error_string += "; the particles passed to the function are " + molecule1.name + ", " + molecule2.name;
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    if (model == models_prob_vv::model_prob_vv_fho) {

      int abs_delta = abs(delta_i);
      double omega1 = abs(molecule1.vibr_energy[e1][i] - molecule1.vibr_energy[e1][i + delta_i]) / abs_delta;
      double omega2 = abs(molecule2.vibr_energy[e2][k] - molecule2.vibr_energy[e2][k - delta_i]) / abs_delta;
 
      return probability_VV_FHO(rel_vel, interaction.collision_mass, molecule1.vibr_energy[e1][i] 		- 
                                                                     molecule1.vibr_energy[e1][i + delta_i] 	+ 
                                                                     molecule2.vibr_energy[e2][k]               - 
                                                                     molecule2.vibr_energy[e2][k - delta_i],
                                                                     i, k, delta_i, omega1, omega2, abs(omega1 - omega2), interaction["alpha_FHO"]);
    } else {
      std::string error_string = "Unknown choice of VV exchange model";
      throw kappa::ModelParameterException(error_string.c_str());
    }
  }
 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::probability_VT(	double rel_vel, 
		 				kappa::Molecule const &molecule, 
						kappa::Interaction const &interaction,
                                             	int i, int delta_i, int e, 
						kappa::models_prob_vt model) {

    #ifdef KAPPA_STRICT_CHECKS
    if (rel_vel < 0.0) {
      std::string error_string = "Negative relative velocity: rel_vel=" + std::to_string(rel_vel);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e < 0) {
      std::string error_string = "Negative electronic level: e=" + std::to_string(e);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (delta_i == 0) {
      std::string error_string = "No change in vibrational levels of molecule after VT transitions";
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if ((i < 0) || (i + delta_i < 0)) {
      std::string error_string = "One of the vibrational levels before or after the VT transition is negative: i, i'=";
      error_string += std::to_string(i) + ", " + std::to_string(i + delta_i);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if ((interaction.particle1_name != molecule.name) && (interaction.particle2_name != molecule.name)) {
      std::string error_string = "Interaction is not for the molecule passed for the function, the interaction is for: " + interaction.particle1_name + " + " 
                                                                                                                         + interaction.particle2_name;
      error_string += "; the molecule passed to the function is " + molecule.name;
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    // FHO model for VT energy exchange
    if (model == models_prob_vt::model_prob_vt_fho) {

      // FIXME BUG with abs on linux -- workaround for the problematic abs
      double omega = 0.;
      if ( (molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i]) > 0. ) {
        omega = (molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i]) / (K_CONST_HBAR * abs(delta_i));
      } else {
        omega = - (molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i]) / (K_CONST_HBAR * abs(delta_i));
      }

      double res = probability_VT_FHO(	rel_vel, 
					interaction.collision_mass, 
					molecule.reduced_osc_mass,
		                        molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i],
 	                          	i, delta_i, 	
					omega, 
					molecule.mA_mAB, 
					interaction["alpha_FHO"], 
					interaction["E_FHO"], 
					interaction["SVT_FHO"]);

      if (molecule.rot_symmetry == 1) {

        res = 0.5 * (res + probability_VT_FHO(	rel_vel, 
  						interaction.collision_mass, 
						molecule.reduced_osc_mass,
                        			molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i],
                       				i, delta_i, 
						omega, 	
						molecule.mB_mAB, 
						interaction["alpha_FHO"], 
						interaction["E_FHO"], 
						interaction["SVT_FHO"]));
      }
      return res;
    } else {
      std::string error_string = "Unknown choice of VT transition model";
      throw kappa::ModelParameterException(error_string.c_str());
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::probability_diss(double rel_vel, 
						kappa::Molecule const &molecule, 
						kappa::Interaction const &interaction, 
						int i, 
						int e, 
						kappa::models_prob_diss model) {
    #ifdef KAPPA_STRICT_CHECKS
    if (rel_vel < 0.0) {
      std::string error_string = "Negative relative velocity: rel_vel=" + std::to_string(rel_vel);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e < 0) {
      std::string error_string = "Negative electronic level: e=" + std::to_string(e);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (i < 0) {
      std::string error_string = "Negative vibrational level: i=";
      error_string += std::to_string(i);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if ((interaction.particle1_name != molecule.name) && (interaction.particle2_name != molecule.name)) {
      std::string error_string = "Interaction is not for the molecule passed for the function, the interaction is for: " + interaction.particle1_name + " + " 
                                											 + interaction.particle2_name;
      error_string += "; the molecule passed to the function is " + molecule.name;
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    switch (model) {
      case models_prob_diss::model_prob_diss_thresh_cmass_vibr:
        return p_probability_diss(rel_vel, interaction.collision_mass, molecule.diss_energy[e], molecule.vibr_energy[e][i], true);
	break;
      case models_prob_diss::model_prob_diss_thresh_vibr:
	return p_probability_diss(rel_vel, interaction.collision_mass, molecule.diss_energy[e], molecule.vibr_energy[e][i], false);
	break;
      case models_prob_diss::model_prob_diss_thresh_cmass:
	return p_probability_diss(rel_vel, interaction.collision_mass, molecule.diss_energy[e], true);
	break;
      case models_prob_diss::model_prob_diss_thresh:
	return p_probability_diss(rel_vel, interaction.collision_mass, molecule.diss_energy[e], false);
	break;
      default:
        std::string error_string = "Unknown choice of dissociation probability model";
        throw kappa::ModelParameterException(error_string.c_str());
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::probability_diss(double rel_vel, 
						kappa::Molecule const &molecule, 
						kappa::Interaction const &interaction, 
						int i, 
						kappa::models_prob_diss model) {

    return probability_diss(rel_vel, molecule, interaction, i, 0, model);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::crosssection_elastic(double rel_vel, kappa::Interaction const &interaction, kappa::models_cs_elastic model) {

    #ifdef KAPPA_STRICT_CHECKS
    if (rel_vel < 0.0) {
      std::string error_string = "Negative relative velocity: rel_vel=" + std::to_string(rel_vel);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif
  
    switch (model) {
      case models_cs_elastic::model_cs_el_vss:
        if (interaction.vss_data) {
          return crosssection_elastic_VSS(rel_vel, interaction.vss_c_cs, interaction.vss_omega);
        } else {
          std::string error_string = "No VSS data found for " + interaction.particle1_name + " + " + interaction.particle2_name + " interaction";
          throw kappa::DataNotFoundException(error_string.c_str());
        }
	break;
      case models_cs_elastic::model_cs_el_rs:
	return crosssection_elastic_RS(interaction.collision_diameter);
	break;
      default:
        std::string error_string = "Unknown choice of dissociation probability model";
        throw kappa::ModelParameterException(error_string.c_str());
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::crosssection_VT(double rel_vel, kappa::Molecule const &molecule, kappa::Interaction const &interaction,
                                               int i, int delta_i, int e, kappa::models_cs_vt model) {

    #ifdef KAPPA_STRICT_CHECKS
    if (rel_vel < 0.0) {
      std::string error_string = "Negative relative velocity: rel_vel=" + std::to_string(rel_vel);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e < 0) {
      std::string error_string = "Negative electronic level: e=" + std::to_string(e);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (delta_i == 0) {
      std::string error_string = "No change in vibrational levels of molecule after VT transitions";
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if ((i < 0) || (i + delta_i < 0)) {
      std::string error_string = "One of the vibrational levels before or after the VT transition is negative: i, i'=";
      error_string += std::to_string(i) + ", " + std::to_string(i + delta_i);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if ((interaction.particle1_name != molecule.name) && (interaction.particle2_name != molecule.name)) {
      std::string error_string = "Interaction is not for the molecule passed for the function, the interaction is for: " + interaction.particle1_name + " + " + interaction.particle2_name;
      error_string += "; the molecule passed to the function is " + molecule.name;
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    switch(model) {
      case models_cs_vt::model_cs_vt_rs_fho: {
        double omega = abs(molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i]) / (K_CONST_HBAR * abs(delta_i));

        double res = crosssection_VT_FHO_RS(rel_vel, interaction.collision_mass, interaction.collision_diameter, molecule.reduced_osc_mass,
                                          molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i],
                                          i, delta_i, omega, molecule.mA_mAB, interaction["alpha_FHO"], interaction["E_FHO"], interaction["SVT_FHO"]);


        if (molecule.rot_symmetry == 1) {
          res = 0.5 * (res + crosssection_VT_FHO_RS(rel_vel, interaction.collision_mass, interaction.collision_diameter, molecule.reduced_osc_mass,
                       molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i],
                       i, delta_i, omega, molecule.mB_mAB, interaction["alpha_FHO"], interaction["E_FHO"], interaction["SVT_FHO"]));
        }
        return res;
        break;
      }

      case models_cs_vt::model_cs_vt_vss_fho: {
        double omega = abs(molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i]) / (K_CONST_HBAR * abs(delta_i));
        if (interaction.vss_data) {

          double res = crosssection_VT_FHO_VSS(rel_vel, interaction.collision_mass, interaction.vss_c_cs, interaction.vss_omega, molecule.reduced_osc_mass,
                                               molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i],
                                               i, delta_i, omega, molecule.mA_mAB, interaction["alpha_FHO"], interaction["E_FHO"], interaction["SVT_FHO"]);

          if (molecule.rot_symmetry == 1) {
            res = 0.5 * (res + crosssection_VT_FHO_VSS(rel_vel, interaction.collision_mass, interaction.vss_c_cs, interaction.vss_omega, molecule.reduced_osc_mass,
                  molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i],
                  i, delta_i, omega, molecule.mB_mAB, interaction["alpha_FHO"], interaction["E_FHO"], interaction["SVT_FHO"]));
          }
          return res;
        } else {
          std::string error_string = "No VSS data found for " + interaction.particle1_name + " + " + interaction.particle2_name + " interaction";
          throw kappa::DataNotFoundException(error_string.c_str());
        }
      
        break;
      }
      default:
        std::string error_string = "Unknown choice of VT transition model";
        throw kappa::ModelParameterException(error_string.c_str());
    }
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::crosssection_VT(double rel_vel, kappa::Molecule const &molecule, kappa::Interaction const &interaction,
                                              int i, int delta_i, kappa::models_cs_vt model) {

    return crosssection_VT(rel_vel, molecule, interaction, i, delta_i, 0, model);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::crosssection_diss(	double rel_vel, 
							kappa::Molecule const &molecule, 
							kappa::Interaction const &interaction, 
							int i, 
							int e, 
							kappa::models_cs_diss model) {

    #ifdef KAPPA_STRICT_CHECKS
    if (rel_vel < 0.0) {
      std::string error_string = "Negative relative velocity: rel_vel=" + std::to_string(rel_vel);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e < 0) {
      std::string error_string = "Negative electronic level: e=" + std::to_string(e);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (i < 0) {
      std::string error_string = "Negative vibrational level: i=";
      error_string += std::to_string(i);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if ((interaction.particle1_name != molecule.name) && (interaction.particle2_name != molecule.name)) {
      std::string error_string = "Interaction is not for the molecule passed for the function, the interaction is for: " + interaction.particle1_name + " + " + interaction.particle2_name;
      error_string += "; the molecule passed to the function is " + molecule.name;
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    switch (model) {
      case models_cs_diss::model_cs_diss_rs_thresh_cmass_vibr:
	return crosssection_diss_RS(rel_vel, interaction.collision_mass, interaction.collision_diameter, molecule.diss_energy[e], molecule.vibr_energy[e][i], true);
 	break;
      case models_cs_diss::model_cs_diss_rs_thresh_vibr:
	return crosssection_diss_RS(rel_vel, interaction.collision_mass, interaction.collision_diameter, molecule.diss_energy[e], molecule.vibr_energy[e][i], false);
	break;
      case models_cs_diss::model_cs_diss_rs_thresh_cmass:
	return crosssection_diss_RS(rel_vel, interaction.collision_mass, interaction.collision_diameter, molecule.diss_energy[e], true);
	break;
      case models_cs_diss::model_cs_diss_rs_thresh:
	return crosssection_diss_RS(rel_vel, interaction.collision_mass, interaction.collision_diameter, molecule.diss_energy[e], false);
	break;
      case models_cs_diss::model_cs_diss_vss_thresh_cmass_vibr:
        if (interaction.vss_data) {
          return crosssection_diss_VSS(rel_vel, interaction.collision_mass, interaction.vss_c_cs, interaction.vss_omega, molecule.diss_energy[e], molecule.vibr_energy[e][i], true);
        } else {
          std::string error_string = "No VSS data found for " + interaction.particle1_name + " + " + interaction.particle2_name + " interaction";
          throw kappa::DataNotFoundException(error_string.c_str());
        }
	break;
      case models_cs_diss::model_cs_diss_vss_thresh_vibr:
        if (interaction.vss_data) {
          return crosssection_diss_VSS(rel_vel, interaction.collision_mass, interaction.vss_c_cs, interaction.vss_omega, molecule.diss_energy[e], molecule.vibr_energy[e][i], false);
        } else {
          std::string error_string = "No VSS data found for " + interaction.particle1_name + " + " + interaction.particle2_name + " interaction";
          throw kappa::DataNotFoundException(error_string.c_str());
        }
	break;
      case models_cs_diss::model_cs_diss_vss_thresh_cmass:
        if (interaction.vss_data) {
          return crosssection_diss_VSS(rel_vel, interaction.collision_mass, interaction.vss_c_cs, interaction.vss_omega, molecule.diss_energy[e], true);
        } else {
          std::string error_string = "No VSS data found for " + interaction.particle1_name + " + " + interaction.particle2_name + " interaction";
          throw kappa::DataNotFoundException(error_string.c_str());
        }
	break;
      case models_cs_diss::model_cs_diss_vss_thresh:
        if (interaction.vss_data) {
          return crosssection_diss_VSS(rel_vel, interaction.collision_mass, interaction.vss_c_cs, interaction.vss_omega, molecule.diss_energy[e], false);
        } else {
          std::string error_string = "No VSS data found for " + interaction.particle1_name + " + " + interaction.particle2_name + " interaction";
          throw kappa::DataNotFoundException(error_string.c_str());
        }
	break;
      case models_cs_diss::model_cs_diss_ilt:
        if (e==0) {
          if (((interaction.particle1_name == "N2") && (interaction.particle2_name == "N")) || ((interaction.particle1_name == "N") && (interaction.particle2_name == "N2"))) {
            return crosssection_diss_ILT_N2N(rel_vel, interaction.collision_mass, molecule.vibr_energy[e][i], convert_vibr_ladder_N2(molecule.vibr_energy[e][i]));
          } else if (((interaction.particle1_name == "O2") && (interaction.particle2_name == "O")) || ((interaction.particle1_name == "O") && (interaction.particle2_name == "O2"))) {
            return crosssection_diss_ILT_O2O(rel_vel, interaction.collision_mass, molecule.vibr_energy[e][i], convert_vibr_ladder_O2(molecule.vibr_energy[e][i]));
          } else {
            std::string error_string = "No ILT model available for dissociation for the interaction: " + interaction.particle1_name + "+" + interaction.particle2_name;
            throw kappa::ModelParameterException(error_string.c_str());
          }
        } else {
          std::string error_string = "No ILT model available for dissociation from excited electronic states";
          throw kappa::ModelParameterException(error_string.c_str());
        }
        break;
      default:
        std::string error_string = "Unknown choice of dissociation probability model";
        throw kappa::ModelParameterException(error_string.c_str());
     }  
  }  

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::crosssection_diss(	double rel_vel, 	
							kappa::Molecule const &molecule, 
							kappa::Interaction const &interaction, 
							int i, 
							kappa::models_cs_diss model) {

    return crosssection_diss(rel_vel, molecule, interaction, i, 0, model);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::Z_coll(double T, double n, kappa::Interaction const &interaction) {

    #ifdef KAPPA_STRICT_CHECKS
    if (T <= 0) {
      std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (n <= 0) {
      std::string error_string = "Non-positive number density: n=" + std::to_string(n);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    return p_Z_coll(T, n, interaction.collision_mass, interaction.collision_diameter);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::rot_collision_number_parker(double T, kappa::Molecule const &molecule, kappa::Interaction const &interaction) {

    #ifdef KAPPA_STRICT_CHECKS
    if (T <= 0) {
      std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if ((interaction.particle1_name != molecule.name) && (interaction.particle2_name != molecule.name)) {
      std::string error_string = "Interaction is not for the molecule passed for the function, the interaction is for: " + interaction.particle1_name + 
 														   " + " + interaction.particle2_name;
      error_string += "; the molecule passed to the function is " + molecule.name;
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    return p_rot_collision_number_parker(T, molecule.parker_const, interaction.epsilon);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Rotational relaxation time by Parker
   double kappa::Approximation::rot_relaxation_time_parker(	double T, 
								double n, 
								kappa::Molecule const &molecule, 
								kappa::Interaction const &interaction, 
								kappa::models_omega model) {

     #ifdef KAPPA_STRICT_CHECKS
     if (T <= 0) {
       std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
       throw kappa::IncorrectValueException(error_string.c_str());
     }
     if (n <= 0) {
       std::string error_string = "Non-positive number density: n=" + std::to_string(n);
       throw kappa::IncorrectValueException(error_string.c_str());
     }
     if ((interaction.particle1_name != molecule.name) && (interaction.particle2_name != molecule.name)) {
       std::string error_string = "Interaction is not for the molecule passed for the function, the interaction is for: " + interaction.particle1_name 
  														  + " + " + interaction.particle2_name;
       error_string += "; the molecule passed to the function is " + molecule.name;
       throw kappa::IncorrectValueException(error_string.c_str());
     }
     #endif

     return p_rot_collision_number_parker(T, molecule.parker_const, interaction.epsilon) * K_CONST_PI * 0.15625 / 
	                                 (n * omega_integral(T, interaction, 2, 2, model, true));
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Vibrational relaxation time by Millikan-White
  double kappa::Approximation::vibr_relaxation_time_MW(double T, double n, kappa::Molecule const &molecule, kappa::Interaction const &interaction) {

    #ifdef KAPPA_STRICT_CHECKS
    if (T <= 0) {
      std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (n <= 0) {
      std::string error_string = "Non-positive number density: n=" + std::to_string(n);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if ((interaction.particle1_name != molecule.name) && (interaction.particle2_name != molecule.name)) {
      std::string error_string = "Interaction is not for the molecule passed for the function, the interaction is for: " + interaction.particle1_name + 
														   " + " + interaction.particle2_name;
      error_string += "; the molecule passed to the function is " + molecule.name;
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    return p_vibr_relaxation_time_MW(T, n, molecule.characteristic_vibr_temperatures[0], interaction.collision_mass);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::vibr_relaxation_time_Park_corr(double T, double n, kappa::Interaction const &interaction, double crosssection) {

    #ifdef KAPPA_STRICT_CHECKS
    if (T <= 0) {
      std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (n <= 0) {
      std::string error_string = "Non-positive number density: n=" + std::to_string(n);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (crosssection <= 0) {
      std::string error_string = "Non-positive crosssection: crosssection=" + std::to_string(n);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    return p_vibr_relaxation_time_Park_corr(T, n, interaction.collision_mass, crosssection);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::omega_integral(double T, kappa::Interaction const &interaction, int l, int r, kappa::models_omega model, bool dimensional) {
   
    #ifdef KAPPA_STRICT_CHECKS
    if (T <= 0) {
      std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if ((l <= 0) || (r <= 0)) {
      std::string error_string = "Non-positive number density l or degree: (l,r)=(" + std::to_string(l) + "," + std::to_string(r) + ")";
      throw kappa::IncorrectValueException(error_string.c_str());
    } 
    #endif

    switch (model) {

      // Rigid Sphere (RS) model
      case models_omega::model_omega_rs:
        if (dimensional == true) {
          return omega_integral_RS(T, l, r, interaction.collision_diameter, interaction.collision_mass);
        } else {
          return 1.0;
        }
        break;

      // VSS model
      case models_omega::model_omega_vss:
        if (interaction.vss_data) {
          if ((l>5)) {
            std::string error_string = "Incorrect value of l specified: can only be l<=5 for the VSS potential";
            throw kappa::ModelParameterException(error_string.c_str());
          } else {
            if (dimensional) {
              return omega_integral_VSS(T, l, r, interaction.collision_mass, interaction.vss_c_cs, interaction.vss_omega, interaction.vss_alpha);
            } else {
              return omega_integral_VSS(T, l, r, interaction.collision_mass, interaction.vss_c_cs, interaction.vss_omega, interaction.vss_alpha) / 
              omega_integral_RS(T, l, r, interaction.collision_diameter, interaction.collision_mass);
            }
          }
        } else {
          std::string error_string = "No VSS data found for " + interaction.particle1_name + " + " + interaction.particle2_name + " interaction";
          throw kappa::DataNotFoundException(error_string.c_str());
        }
        break;

      // Bracket integrals -- Lennard-Jones potential
      case models_omega::model_omega_lennardjones:
        if ((l == 1 && r == 1) || (l == 2 && r == 2)) {
          if (dimensional == true) {
	    return omega_integral_LennardJones(T, l, r, interaction.epsilon) * omega_integral_RS(T, l, r, interaction.collision_diameter, interaction.collision_mass);
          } else {
            return omega_integral_LennardJones(T, l, r, interaction.epsilon);
          }
        } else if ((l == 1 && r == 2) || (l == 2 && r == 3)) {

          double res =  (omega_integral_LennardJones(T+K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.epsilon) *
	                 omega_integral_RS(T+K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.collision_diameter, interaction.collision_mass) -
	                 omega_integral_LennardJones(T-K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.epsilon) *
	                 omega_integral_RS(T-K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.collision_diameter, interaction.collision_mass)) * T / (2 * K_CONST_OMEGA_D_STEP_SIZE) +
	                 (l + 1.5) * omega_integral_LennardJones(T, l, l, interaction.epsilon) * 
	                 omega_integral_RS(T, l, l, interaction.collision_diameter, interaction.collision_mass);

          if (dimensional == true) {
  	    return res;
	  } else {
	    return res / omega_integral_RS(T, l, r, interaction.collision_diameter, interaction.collision_mass);
	  }
        } else if ((l == 1 && r == 3) || (l == 2 && r == 4)) {

          double r0 = omega_integral_LennardJones(T, l, l, interaction.epsilon) * omega_integral_RS(T, l, l, interaction.collision_diameter, interaction.collision_mass);

          double rp1 = omega_integral_LennardJones(T+K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.epsilon)
                     * omega_integral_RS(T+K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.collision_diameter, interaction.collision_mass);

          double rm1 = omega_integral_LennardJones(T-K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.epsilon)
                     * omega_integral_RS(T-K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.collision_diameter, interaction.collision_mass);

          r0 = T * T * (rp1 - 2 * r0 + rm1) / (K_CONST_OMEGA_D_STEP_SIZE * K_CONST_OMEGA_D_STEP_SIZE) + T * (2 * l + 5) * (rp1 - rm1) / (2 * K_CONST_OMEGA_D_STEP_SIZE) + (l + 2.5) * (l + 1.5) * r0;

	  if (dimensional == true) {
            return r0;
          } else {
            return r0 / omega_integral_RS(T, l, r, interaction.collision_diameter, interaction.collision_mass);
          }
        } else if (l == 1 && r == 4) {

          double r0 = omega_integral_LennardJones(T, l, l, interaction.epsilon)
                    * omega_integral_RS(T, l, l, interaction.collision_diameter, interaction.collision_mass);

          double rp1 = omega_integral_LennardJones(T+K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.collision_mass)
                     * omega_integral_RS(T+K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.collision_diameter, interaction.collision_mass);

          double rm1 = omega_integral_LennardJones(T-K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.epsilon)
                     * omega_integral_RS(T-K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.collision_diameter, interaction.collision_mass);

          double rp2 = omega_integral_LennardJones(T+2*K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.epsilon)
                     * omega_integral_RS(T+2*K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.collision_diameter, interaction.collision_mass);

          double rm2 = omega_integral_LennardJones(T-2*K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.epsilon)
                     * omega_integral_RS(T-2*K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.collision_diameter, interaction.collision_mass);

          r0 = T * T * T * (rp2 - 2 * rp1 + 2 * rm1 - rm2) / (2 * K_CONST_OMEGA_D_STEP_SIZE * K_CONST_OMEGA_D_STEP_SIZE * K_CONST_OMEGA_D_STEP_SIZE)
             + T * T * (3 * l + 14.5) * (rp1 - 2 * r0 + rm1) / (K_CONST_OMEGA_D_STEP_SIZE * K_CONST_OMEGA_D_STEP_SIZE)
             + T * ((2 * l + 9) * (l + 4.5) + (l + 2.5) * (l + 1.5)) * (rp1 - rm1) / (2 * K_CONST_OMEGA_D_STEP_SIZE)
             + (l + 3.5) * (l + 2.5) * (l + 1.5) * r0;               

          if (dimensional == true) {
            return r0;
          } else {
              return r0 / omega_integral_RS(T, l, r, interaction.collision_diameter, interaction.collision_mass);
          }
        } else if (l == 1 && r == 5) {

          double r0 = omega_integral(T, interaction, 1, 4, model, true);

          double res = (omega_integral(T+K_CONST_OMEGA_D_STEP_SIZE, interaction, 1, 4, model, true) -
                        omega_integral(T-K_CONST_OMEGA_D_STEP_SIZE, interaction, 1, 4, model, true)) * T / (2 * K_CONST_OMEGA_D_STEP_SIZE) + 5.5 * r0;             

          if (dimensional == true) {
            return r0;
          } else {
            return r0 / omega_integral_RS(T, l, r, interaction.collision_diameter, interaction.collision_mass);
          }
        } else {
	  std::string error_string = "Incorrect values of (l, r) specified: (" + std::to_string(l) + ", " + std::to_string(r) + "); can only be (1, 1), (1, 2), (1, 3), (1, 4), (1, 5)"
                                   + "(2, 2), (2, 3), (2, 4) for the Lennard-Jones potential";
	  throw kappa::ModelParameterException(error_string.c_str());
        }
	break;

      // Bracket integrals -- Born-Meyer potential
      case models_omega::model_omega_bornmayer:
        if ((l == 1 && r == 1) || (l == 2 && r == 2)) {
	  if (dimensional == true) {
	    return omega_integral_Born_Mayer(T, l, r, interaction["beta Born-Mayer"], interaction["phi_zero Born-Mayer"], interaction.collision_diameter, interaction["Lennard-Jones epsilon"])
	         * omega_integral_RS(T, l, r, interaction.collision_diameter, interaction.collision_mass);
	  } else {
	    return omega_integral_Born_Mayer(T, l, r, interaction["beta Born-Mayer"], interaction["phi_zero Born-Mayer"], interaction.collision_diameter, interaction["Lennard-Jones epsilon"]);
	  }
	 } else if ((l == 1 && r == 2) || (l == 2 && r == 3)) {

	   double res =  (omega_integral_Born_Mayer(T+K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction["beta Born-Mayer"],
	                  							       interaction["phi_zero Born-Mayer"], 
										       interaction.collision_diameter, 
										       interaction["Lennard-Jones epsilon"])
	                * omega_integral_RS(T+K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.collision_diameter, interaction.collision_mass) -
	                  omega_integral_Born_Mayer(T-K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction["beta Born-Mayer"],
	                  							       interaction["phi_zero Born-Mayer"], 
										       interaction.collision_diameter, 
									               interaction["Lennard-Jones epsilon"])
	                * omega_integral_RS(T-K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.collision_diameter, interaction.collision_mass)) * T / 
			(2 * K_CONST_OMEGA_D_STEP_SIZE) + (l + 1.5) * omega_integral_Born_Mayer(T, l, l, interaction["beta Born-Mayer"],
	             	  										 interaction["phi_zero Born-Mayer"], 
													 interaction.collision_diameter, 
													 interaction["Lennard-Jones epsilon"])
	                * omega_integral_RS(T, l, l, interaction.collision_diameter, interaction.collision_mass);

	   if (dimensional == true) {
             return res;
           } else {
             return res / omega_integral_RS(T, l, r, interaction.collision_diameter, interaction.collision_mass);
           }
         } else if ((l == 1 && r == 3) || (l == 2 && r == 4)) {

           double r0 = omega_integral_Born_Mayer(T, l, l, interaction["beta Born-Mayer"], 
 						          interaction["phi_zero Born-Mayer"], 
							  interaction.collision_diameter, 
							  interaction["Lennard-Jones epsilon"])
	              * omega_integral_RS(T, l, l, interaction.collision_diameter, interaction.collision_mass);

           double rp1 = omega_integral_Born_Mayer(T+K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction["beta Born-Mayer"],
	                                                   			     interaction["phi_zero Born-Mayer"], 
										     interaction.collision_diameter, 
										     interaction["Lennard-Jones epsilon"])
	              * omega_integral_RS(T+K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.collision_diameter, interaction.collision_mass);

           double rm1 = omega_integral_Born_Mayer(T-K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction["beta Born-Mayer"],
	                                                   			     interaction["phi_zero Born-Mayer"], 
										     interaction.collision_diameter, 
										     interaction["Lennard-Jones epsilon"])
	              * omega_integral_RS(T-K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.collision_diameter, interaction.collision_mass);

	   r0 = T * T * (rp1 - 2 * r0 + rm1) / (K_CONST_OMEGA_D_STEP_SIZE * K_CONST_OMEGA_D_STEP_SIZE) + T * (2 * l + 5) * (rp1 - rm1) / (2 * K_CONST_OMEGA_D_STEP_SIZE) + (l + 2.5) * (l + 1.5) * r0;             	

           if (dimensional == true) {
             return r0;
           }  else {
              return r0 / omega_integral_RS(T, l, r, interaction.collision_diameter, interaction.collision_mass);
           }
         } else if (l == 1 && r == 4) {

           double r0 = omega_integral_Born_Mayer(T, l, l, interaction["beta Born-Mayer"],
 	                                                  interaction["phi_zero Born-Mayer"], 
	 					          interaction.collision_diameter, 
					                  interaction["Lennard-Jones epsilon"])
	                       			          * omega_integral_RS(T, l, l, interaction.collision_diameter, interaction.collision_mass);

           double rp1 = omega_integral_Born_Mayer(T+K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction["beta Born-Mayer"],
                                                      			             interaction["phi_zero Born-Mayer"], 
										     interaction.collision_diameter, 
										     interaction["Lennard-Jones epsilon"])
	                         						     * omega_integral_RS(T + K_CONST_OMEGA_D_STEP_SIZE, 
									             l, l, 
										     interaction.collision_diameter, 
										     interaction.collision_mass);

           double rm1 = omega_integral_Born_Mayer(T-K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction["beta Born-Mayer"],
	                                                   			      interaction["phi_zero Born-Mayer"], 
										      interaction.collision_diameter, 
										      interaction["Lennard-Jones epsilon"])
						                                      * omega_integral_RS(T-K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.collision_diameter, interaction.collision_mass);

           double rp2 = omega_integral_Born_Mayer(T+2*K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction["beta Born-Mayer"],
                                                  interaction["phi_zero Born-Mayer"], interaction.collision_diameter, interaction["Lennard-Jones epsilon"])
	                                          * omega_integral_RS(T+2*K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.collision_diameter, interaction.collision_mass);

           double rm2 = omega_integral_Born_Mayer(T-2*K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction["beta Born-Mayer"],
	                                          interaction["phi_zero Born-Mayer"], interaction.collision_diameter, interaction["Lennard-Jones epsilon"])
	                                          * omega_integral_RS(T-2*K_CONST_OMEGA_D_STEP_SIZE, l, l, interaction.collision_diameter, interaction.collision_mass);

           r0 = T * T * T * (rp2 - 2 * rp1 + 2 * rm1 - rm2) / (2 * K_CONST_OMEGA_D_STEP_SIZE * K_CONST_OMEGA_D_STEP_SIZE * K_CONST_OMEGA_D_STEP_SIZE)
	      + T * T * (3 * l + 14.5) * (rp1 - 2 * r0 + rm1) / (K_CONST_OMEGA_D_STEP_SIZE * K_CONST_OMEGA_D_STEP_SIZE)
	      + T * ((2 * l + 9) * (l + 4.5) + (l + 2.5) * (l + 1.5)) * (rp1 - rm1) / (2 * K_CONST_OMEGA_D_STEP_SIZE)
	      + (l + 3.5) * (l + 2.5) * (l + 1.5) * r0;               

           if (dimensional == true) {
             return r0;
           } else {
	     return r0 / omega_integral_RS(T, l, r, interaction.collision_diameter, interaction.collision_mass);
	   }
	 } else if (l == 1 && r == 5) {
	   double r0 = omega_integral(T, interaction, 1, 4, model, true);
	   double res = (omega_integral(T+K_CONST_OMEGA_D_STEP_SIZE, interaction, 1, 4, model, true) -
	                 omega_integral(T-K_CONST_OMEGA_D_STEP_SIZE, interaction, 1, 4, model, true)) * T / (2 * K_CONST_OMEGA_D_STEP_SIZE)
	                       + 5.5 * r0;           

           if (dimensional == true) {
             return r0;
           } else {
	     return r0 / omega_integral_RS(T, l, r, interaction.collision_diameter, interaction.collision_mass);
	   }
	 } else {
	   std::string error_string = "Incorrect values of (l, r) specified: (" + std::to_string(l) + ", " + std::to_string(r) + "); can only be (1, 1), (1, 2), (1, 3), (1, 4), (1, 5)"
                                    + "(2, 2), (2, 3), (2, 4) for the Born-Mayer potential";
	   throw kappa::ModelParameterException(error_string.c_str());
	 }
	 break;
 	
      // Bracket integrals -- ESA model
      case models_omega::model_omega_esa:
        if (interaction.interaction_type == interaction_types::interaction_neutral_neutral) {
          if ((interaction.particle1_name == "H2") && (interaction.particle2_name == "H2")) {
            if ((l == 1 && r == 1) || (l == 1 && r == 2) || (l == 1 && r == 3) || (l == 1 && r == 4) || (l == 1 && r == 5) || (l == 2 && r == 2) || (l == 2 && r == 3) || (l == 2 && r == 4) || (l == 3 && r == 3)) {
              if (dimensional == true) {
                return omega_integral_ESA_H2H2(T, l, r, interaction.collision_diameter) * omega_integral_RS(T, l, r, interaction.collision_diameter, interaction.collision_mass);
              } else {
                return omega_integral_ESA_H2H2(T, l, r, interaction.collision_diameter);
              }
            } else {
              std::string error_string = "Incorrect values of (l, r) specified: (" + std::to_string(l) + ", " + std::to_string(r) + ");" +
                " can only be (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (2, 2), (2, 3), (2, 4), (3, 3), for the ESA potential for H2+H2 interactions";
              throw kappa::ModelParameterException(error_string.c_str());
            }
          } else if (((interaction.particle1_name == "H2") && (interaction.particle2_name == "H")) || ((interaction.particle1_name == "H") && (interaction.particle2_name == "H2"))) {
            if ((l == 1 && r == 1) || (l == 1 && r == 2) || (l == 1 && r == 3) || (l == 1 && r == 4) || (l == 1 && r == 5) || (l == 2 && r == 2) || (l == 2 && r == 3) || (l == 2 && r == 4) || (l == 3 && r == 3)) {
              if (dimensional == true) {
                return omega_integral_ESA_H2H(T, l, r, interaction.collision_diameter) * omega_integral_RS(T, l, r, interaction.collision_diameter, interaction.collision_mass);
              } else {
                return omega_integral_ESA_H2H(T, l, r, interaction.collision_diameter);
              }
            } else {
              std::string error_string = "Incorrect values of (l, r) specified: (" + std::to_string(l) + ", " + std::to_string(r) + ");" +
                " can only be (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (2, 2), (2, 3), (2, 4), (3, 3), for the ESA potential for H2+H interactions";
              throw kappa::ModelParameterException(error_string.c_str());
            }
          } else if ((interaction.particle1_name == "H") && (interaction.particle2_name == "H")) {
            if ((l == 1 && r == 1) || (l == 1 && r == 2) || (l == 1 && r == 3) || (l == 1 && r == 4) || (l == 1 && r == 5) || (l == 2 && r == 2) || (l == 2 && r == 3) || (l == 2 && r == 4) || (l == 3 && r == 3)) {
              if (dimensional == true) {
                return omega_integral_ESA_H2H2(T, l, r, interaction.collision_diameter) * omega_integral_RS(T, l, r, interaction.collision_diameter, interaction.collision_mass);
              } else {
                return omega_integral_ESA_H2H2(T, l, r, interaction.collision_diameter);
              }
            } else {
              std::string error_string = "Incorrect values of (l, r) specified: (" + std::to_string(l) + ", " + std::to_string(r) + ");" +
                " can only be (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (2, 2), (2, 3), (2, 4), (3, 3), for the ESA potential for hydrogen H+H interactions";
              throw kappa::ModelParameterException(error_string.c_str());
            }
          } else if ((l == 1 && r == 1) || (l == 1 && r == 2) || (l == 1 && r == 3) || (l == 1 && r == 4) || (l == 1 && r == 5) || (l == 2 && r == 2) || (l == 2 && r == 3) || (l == 2 && r == 4) || (l == 3 && r == 3) || (l == 4 && r == 4)) {
  	    if (dimensional == true) {
  	      return omega_integral_ESA_nn(T, l, r, interaction["beta Bruno"], interaction["epsilon0 Bruno"], interaction["re Bruno"]) * omega_integral_RS(T, l, r, interaction.collision_diameter, interaction.collision_mass);
  	    } else {
  	      return omega_integral_ESA_nn(T, l, r, interaction["beta Bruno"], interaction["epsilon0 Bruno"], interaction["re Bruno"]);
  	    }
  	  } else {
  	    std::string error_string = "Incorrect values of (l, r) specified: (" + std::to_string(l) + ", " + std::to_string(r) + ");" +
  						" can only be (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (2, 2), (2, 3), (2, 4), (3, 3), (4, 4) for the ESA potential";
  	    throw kappa::ModelParameterException(error_string.c_str());
  	  }
        } else if (interaction.interaction_type == interaction_types::interaction_neutral_ion) {
  	  if ((l == 1 && r == 1) || (l == 1 && r == 2) || (l == 1 && r == 3) || (l == 1 && r == 4) || (l == 1 && r == 5) || (l == 2 && r == 2) || (l == 2 && r == 3) || (l == 2 && r == 4) || (l == 3 && r == 3) || (l == 4 && r == 4)) {

  	    double tmp = omega_integral_ESA_cn(T, l, r, interaction["beta Bruno"], interaction["epsilon0 Bruno"], interaction["re Bruno"]);
  	    double ch_exch = 0;

    	    if (((interaction.particle1_name == "N") && (interaction.particle2_name == "N+")) || ((interaction.particle1_name == "N+") && (interaction.particle2_name == "N"))) {
    	      if (l == 1 && r == 1) {
    	        ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 65.8423, -4.5492, 7.8608e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    	      } else
    		if (l == 1 && r == 2) {
    		  ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 64.3259, -4.4969, 7.8619e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    	        } else
    		  if (l == 1 && r == 3) {
    		    ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 63.2015, -4.4575, 7.8613e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    		} else
    		  if (l == 1 && r == 4) {
    		    ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 62.3098, -4.4260, 7.8606e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    		} else
    		  if (l == 1 && r == 5) {
    		    ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 61.5722, -4.3997, 7.8601e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    		}
    		tmp = sqrt(tmp*tmp + ch_exch*ch_exch);

    	    } else if (((interaction.particle1_name == "Ar") && (interaction.particle2_name == "Ar+")) || ((interaction.particle1_name == "Ar+") && (interaction.particle2_name == "Ar"))) {
    	      if (l == 1 && r == 1) {
    	        ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 68.6279, -4.3366, 6.8542e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    	      } else
    		if (l == 1 && r == 2) {
    	          ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 67.1824, -4.2908, 6.8533e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    		} else
    		  if (l == 1 && r == 3) {
    		    ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 66.1106, -4.2568, 6.8542e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    		  } else
    		    if (l == 1 && r == 4) {
    		      ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 65.2603, -4.2297, 6.8562e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    		    } else
    		      if (l == 1 && r == 5) {
    		        ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 64.5529, -4.2062, 6.8523e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    		      }
    		      tmp = sqrt(tmp*tmp + ch_exch*ch_exch);           

            } else if (((interaction.particle1_name == "O") && (interaction.particle2_name == "O+")) || ((interaction.particle1_name == "O+") && (interaction.particle2_name == "O"))) {
    	      if (l == 1 && r == 1) {
    	        ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 64.7044, -4.1936, 6.7974e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    	      } else
    	        if (l == 1 && r == 2) {
    		  ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 63.3075, -4.1485, 6.7991e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    		} else
    		  if (l == 1 && r == 3) {
    		    ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 62.2707, -4.1146, 6.8000e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    		  } else
    		    if (l == 1 && r == 4) {
    		      ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 61.4470, -4.0872, 6.7986e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    		    } else
    		      if (l == 1 && r == 5) {
    		        ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 60.7663, -4.0646, 6.7987e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    		      }
    		      tmp = sqrt(tmp*tmp + ch_exch*ch_exch);           

            } else if (((interaction.particle1_name == "C") && (interaction.particle2_name == "C+")) || ((interaction.particle1_name == "C+") && (interaction.particle2_name == "C"))) {
    	      if (l == 1 && r == 1) {
    	        ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 65.8583, -4.8063, 8.7735e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    	      } else
    	        if (l == 1 && r == 2) {
    		  ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 64.2555, -4.7476, 8.7729e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    		} else
    		  if (l == 1 && r == 3) {
    		    ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 63.0684, -4.7037, 8.7724e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    		  } else
    		    if (l == 1 && r == 4) {
    		      ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 62.1279, -4.6687, 8.7729e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    		    } else
    		      if (l == 1 && r == 5) {
    		        ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 61.3500, -4.6395, 8.7730e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    		      }
    		      tmp = sqrt(tmp * tmp + ch_exch * ch_exch);           

    	    } else if (((interaction.particle1_name == "CO") && (interaction.particle2_name == "CO+")) || ((interaction.particle1_name == "CO+") && (interaction.particle2_name == "CO"))) {
    	      if (l == 1 && r == 1) {
    	        ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 85.0889, -5.5980, 9.2122e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    	      } else
    		if (l == 1 && r == 2) {
    		  ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 83.2212, -5.5362, 9.2097e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    	 	} else
    	  	  if (l == 1 && r == 3) {
    		    ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 81.8376, -5.4902, 9.2105e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    		  } else
    		    if (l == 1 && r == 4) {
    		      ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 80.7381, -5.4530, 9.2078e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    		    } else
    		      if (l == 1 && r == 5) {
    		        ch_exch = kappa::Approximation::omega_integral_ESA_cn_corr( T, 79.8298, -5.4225, 9.2091e-2, interaction["beta Bruno"], interaction["re Bruno"]);
    		      }
    		      tmp = sqrt(tmp*tmp + ch_exch*ch_exch);           
    	    }
            if (dimensional == true) {
      	      return tmp * omega_integral_RS(T, l, r, interaction.collision_diameter, interaction.collision_mass);
      	    } else {
    	      return tmp;
      	    }
  	  } else {
    	    std::string error_string = "Incorrect values of (l, r) specified: (" + std::to_string(l) + ", " + std::to_string(r) + ");" +
    	         " can only be (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (2, 2), (2, 3), (2, 4), (3, 3), (4, 4) for the ESA potential";
            throw kappa::ModelParameterException(error_string.c_str());
  	  }
  	} else if  (interaction.interaction_type == interaction_types::interaction_neutral_electron) {

          if ((l == 1 && r == 1) || (l == 1 && r == 2) || (l == 1 && r == 3) || (l == 1 && r == 4) || (l == 1 && r == 5) || 
              (l == 2 && r == 2) || (l == 2 && r == 3) || (l == 2 && r == 4) || (l == 3 && r == 3)) {

            std::string lrc = "_Om" + std::to_string(l) + std::to_string(r) + "_";

            if (dimensional == true) {

              return omega_integral_ESA_ne(	T, 	interaction[lrc+"0"], 
							interaction[lrc+"1"], 
							interaction[lrc+"2"], 
							interaction[lrc+"3"], 
							interaction[lrc+"4"], 
							interaction[lrc+"5"],
                                           		interaction[lrc+"6"], 
							interaction[lrc+"7"], 
							interaction[lrc+"8"], 
							interaction[lrc+"9"],
                                           		interaction.collision_diameter) * 
							omega_integral_RS(T, l, r, interaction.collision_diameter, interaction.collision_mass);
            } else {
              return omega_integral_ESA_ne(	T, 
						interaction[lrc+"0"], 
						interaction[lrc+"1"], 
						interaction[lrc+"2"], 
						interaction[lrc+"3"], 
						interaction[lrc+"4"], 	
						interaction[lrc+"5"],
                                           	interaction[lrc+"6"], 
						interaction[lrc+"7"], 
						interaction[lrc+"8"], 
						interaction[lrc+"9"],
                                           	interaction.collision_diameter);
            }
          } else {
            std::string error_string = "Incorrect values of (l, r) specified: (" + std::to_string(l) + ", " + std::to_string(r) + ");" +
            " can only be (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (2, 2), (2, 3), (2, 4), (3, 3) for the ESA potential for neutral-electron interactions";
            throw kappa::ModelParameterException(error_string.c_str());
          }
        }
      	break;
      default:
	std::string error_string = "Unknown choice of interaction potential";
	throw kappa::ModelParameterException(error_string.c_str());
    }
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::omega_integral(	double T, 
						kappa::Interaction const &interaction, 
						int l, int r, 
						double debye_length, 
						kappa::models_omega model, 
						bool dimensional) {

    if (interaction.interaction_type != interaction_types::interaction_charged_charged) {
      return omega_integral(T, interaction, l, r, model, dimensional);
    } else {
      #ifdef KAPPA_STRICT_CHECKS
      if (T <= 0) {
        std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
        throw kappa::IncorrectValueException(error_string.c_str());
      }
      if (debye_length <= 0) {
        std::string error_string = "Non-positive Debye length: debye_length=" + std::to_string(debye_length);
        throw kappa::IncorrectValueException(error_string.c_str());
      }
      if ((l <= 0) || (r <= 0)) {
        std::string error_string = "Non-positive number density l or degree: (l,r)=(" + std::to_string(l) + "," + std::to_string(r) + ")";
        throw kappa::IncorrectValueException(error_string.c_str());
      }
      #endif
      
      if ((l == 1 && r == 1) || 
          (l == 1 && r == 2) || 
          (l == 1 && r == 3) || 
          (l == 1 && r == 4) || 
          (l == 1 && r == 5) || 
          (l == 2 && r == 2) || 
          (l == 2 && r == 3) || 
          (l == 2 && r == 4) || 
          (l == 3 && r == 3)) {

        if (dimensional) {
          return omega_integral_ESA_cc(T, l, r, interaction.charge1, interaction.charge2, debye_length) * 
                 omega_integral_RS(T, l, r, interaction.collision_diameter, interaction.collision_mass);
        } else {
          return omega_integral_ESA_cc(T, l, r, interaction.charge1, interaction.charge2, debye_length);
        }
      } else {
        std::string error_string = "Incorrect values of (l, r) specified: (" + std::to_string(l) + ", " + std::to_string(r) + ");" +
        " can only be (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (2, 2), (2, 3), (2, 4), (3, 3) for charged-charged interactions";
        throw kappa::ModelParameterException(error_string.c_str());
      }
    }
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::integral_VT(double T, int degree, kappa::Molecule const &molecule, kappa::Interaction const &interaction, int i, int delta_i, int e, kappa::models_cs_vt model) {

    #ifdef KAPPA_STRICT_CHECKS
    if (T <= 0) {
      std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (degree < 0) {
      std::string error_string = "Negative degree of integral: degree=" + std::to_string(degree);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e < 0) {
      std::string error_string = "Negative electronic level: e=" + std::to_string(e);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (delta_i == 0) {
      std::string error_string = "No change in vibrational levels of molecule after VT transitions";
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if ((i < 0) || (i + delta_i < 0)) {
      std::string error_string = "One of the vibrational levels before or after the VT transition is negative: i, i'=";
      error_string += std::to_string(i) + ", " + std::to_string(i + delta_i);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if ((interaction.particle1_name != molecule.name) && (interaction.particle2_name != molecule.name)) {
      std::string error_string = "Interaction is not for the molecule passed for the function, the interaction is for: " + interaction.particle1_name + " + " + interaction.particle2_name;
      error_string += "; the molecule passed to the function is " + molecule.name;
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    switch (model) {

      // Rigid Sphere (RS) + Forced Harmonic Oscillator (FHO)
      case models_cs_vt::model_cs_vt_rs_fho: {
        double omega = abs(molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i]) / (K_CONST_HBAR * abs(delta_i));
        double res = integral_VT_FHO_RS(T, degree, interaction.collision_mass, interaction.collision_diameter, molecule.reduced_osc_mass,
                                        molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i],
                                        i, delta_i, omega, molecule.mA_mAB, interaction["alpha_FHO"], interaction["E_FHO"], interaction["SVT_FHO"]);

        if (molecule.rot_symmetry == 2) {
          return res;
        } else {
          return 0.5 * (res + integral_VT_FHO_RS(T, degree, interaction.collision_mass, interaction.collision_diameter, molecule.reduced_osc_mass,
                        molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i],
                        i, delta_i, omega, molecule.mB_mAB, interaction["alpha_FHO"], interaction["E_FHO"], interaction["SVT_FHO"]));
        }
        break;
      }

      // VSS + FHO
      case models_cs_vt::model_cs_vt_vss_fho: {
        if (interaction.vss_data) {
          double omega = abs(molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i]) / (K_CONST_HBAR * abs(delta_i));
          double res = integral_VT_FHO_VSS(T, degree, interaction.collision_mass, interaction.vss_c_cs, interaction.vss_omega, molecule.reduced_osc_mass,
          				 molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i],
          				 i, delta_i, omega, molecule.mA_mAB, interaction["alpha_FHO"], interaction["E_FHO"], interaction["SVT_FHO"]);

          if (molecule.rot_symmetry == 2) {
            return res;
          } else {
            return 0.5 * (res + integral_VT_FHO_VSS(T, degree, interaction.collision_mass, interaction.vss_c_cs, interaction.vss_omega, molecule.reduced_osc_mass,
                          molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i],
                          i, delta_i, omega, molecule.mA_mAB, interaction["alpha_FHO"], interaction["E_FHO"], interaction["SVT_FHO"]));
          }
          break;
        } else {
          std::string error_string = "No VSS interaction data for " + interaction.particle1_name + " + " + interaction.particle2_name;
          throw kappa::ModelParameterException(error_string.c_str());
        }
        break;
      }
      default:
        std::string error_string = "Unknown choice of interaction potential";
        throw kappa::ModelParameterException(error_string.c_str());
    }
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::integral_VT(double T, int degree, kappa::Molecule const &molecule, kappa::Interaction const &interaction, int i, int delta_i, kappa::models_cs_vt model) {
    return integral_VT(T, degree, molecule, interaction, i, delta_i, 0, model);
   }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // main function for VT
   double kappa::Approximation::k_VT(	double T, 
					kappa::Molecule const &molecule, 
					kappa::Interaction const &interaction, 
					int i, 
					int delta_i, 
					int e, 
					kappa::models_k_vt model) {

    #ifdef KAPPA_STRICT_CHECKS
    if (T <= 0) {
      std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e < 0) {
      std::string error_string = "Negative electronic level: e=" + std::to_string(e);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (delta_i == 0) {
      std::string error_string = "No change in vibrational levels of molecule after VT transitions";
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if ((i < 0) || (i + delta_i < 0)) {
      std::string error_string = "One of the vibrational levels before or after the VT transition is negative: i, i'=";
      error_string += std::to_string(i) + ", " + std::to_string(i + delta_i);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if ((interaction.particle1_name != molecule.name) && (interaction.particle2_name != molecule.name)) {
      std::string error_string = "Interaction is not for the molecule passed for the function, the interaction is for: " 
                                 :+ interaction.particle1_name + " + " + interaction.particle2_name;
      error_string += "; the molecule passed to the function is " + molecule.name;
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    switch (model) {

      // rigid sphere + forced harmonic oscillator model
      case models_k_vt::model_k_vt_rs_fho: {

        if (delta_i < 0) {

          double omega = -(molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i]) / (K_CONST_HBAR * delta_i);
          double res = k_VT_FHO_RS(T, interaction.collision_mass, interaction.collision_diameter, molecule.reduced_osc_mass,
                                   molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i],
                                   i, delta_i, omega, molecule.mA_mAB, interaction["alpha_FHO"], interaction["E_FHO"], interaction["SVT_FHO"]);

          if (molecule.rot_symmetry == 2) {
            return res;
          } else {
            return 0.5 * (res + k_VT_FHO_RS(T, interaction.collision_mass, interaction.collision_diameter, molecule.reduced_osc_mass,
                   molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i],
                   i, delta_i, omega, molecule.mB_mAB, interaction["alpha_FHO"], interaction["E_FHO"], interaction["SVT_FHO"]));
          }
        } else { // compute the reverse coefficient (i->i-delta_i) and multiply by the factor obtained from detailed balance

          double omega = (molecule.vibr_energy[e][i + delta_i] - molecule.vibr_energy[e][i]) / (K_CONST_HBAR * delta_i);
          double res = k_VT_FHO_RS(T, interaction.collision_mass, interaction.collision_diameter, molecule.reduced_osc_mass,
                                   molecule.vibr_energy[e][i + delta_i] - molecule.vibr_energy[e][i],
                                   i + delta_i, -delta_i, omega, molecule.mA_mAB, interaction["alpha_FHO"], interaction["E_FHO"], interaction["SVT_FHO"]);

          if (molecule.rot_symmetry == 1) {
            res = 0.5 * (res + k_VT_FHO_RS(T, interaction.collision_mass, interaction.collision_diameter, molecule.reduced_osc_mass,
                  molecule.vibr_energy[e][i + delta_i] - molecule.vibr_energy[e][i],
                  i + delta_i, -delta_i, omega, molecule.mB_mAB, interaction["alpha_FHO"], interaction["E_FHO"], interaction["SVT_FHO"]));
          }
          return res * k_bf_VT(T, molecule, i + delta_i, -delta_i, e);
        }
        break;
      }

      // VSS + forced harmonic oscillator model
      case models_k_vt::model_k_vt_vss_fho: {

        if (!interaction.vss_data) {
          std::string error_string = "No VSS interaction data for " + interaction.particle1_name + " + " + interaction.particle2_name;
          throw kappa::ModelParameterException(error_string.c_str());
        }
        if (delta_i < 0) {

          double omega = -(molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i]) / (K_CONST_HBAR * delta_i);
          double res = k_VT_FHO_VSS(T, interaction.collision_mass, interaction.vss_c_cs, interaction.vss_omega, molecule.reduced_osc_mass,
                                    molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i],
                                    i, delta_i, omega, molecule.mA_mAB, interaction["alpha_FHO"], interaction["E_FHO"], interaction["SVT_FHO"]);

          if (molecule.rot_symmetry == 2) {
            return res;
          } else {
            return 0.5 * (res + k_VT_FHO_VSS(T, interaction.collision_mass, interaction.vss_c_cs, interaction.vss_omega, molecule.reduced_osc_mass,
                          molecule.vibr_energy[e][i] - molecule.vibr_energy[e][i + delta_i],
                          i, delta_i, omega, molecule.mB_mAB, interaction["alpha_FHO"], interaction["E_FHO"], interaction["SVT_FHO"]));
          }
        } else {

          double omega = (molecule.vibr_energy[e][i + delta_i] - molecule.vibr_energy[e][i]) / (K_CONST_HBAR * delta_i);
          double res = k_VT_FHO_VSS(T, interaction.collision_mass, interaction.vss_c_cs, interaction.vss_omega, molecule.reduced_osc_mass,
                                    molecule.vibr_energy[e][i + delta_i] - molecule.vibr_energy[e][i],
                                    i + delta_i, -delta_i, omega, molecule.mA_mAB, interaction["alpha_FHO"], interaction["E_FHO"], interaction["SVT_FHO"]);

          if (molecule.rot_symmetry == 1) {
            res = 0.5 * (res + k_VT_FHO_VSS(T, interaction.collision_mass, interaction.vss_c_cs, interaction.vss_omega, molecule.reduced_osc_mass,
                         molecule.vibr_energy[e][i + delta_i] - molecule.vibr_energy[e][i],
                         i + delta_i, -delta_i, omega, molecule.mB_mAB, interaction["alpha_FHO"], interaction["E_FHO"], interaction["SVT_FHO"]));
          }
          return res * k_bf_VT(T, molecule, i + delta_i, -delta_i, e);
        }
        break;
      }

      // SSH model
      case models_k_vt::model_k_vt_ssh: {

        if (delta_i == -1) {
          if (molecule.anharmonic_spectrum) {

            return k_VT_SSH(	T, 
      				i, 
      				interaction.collision_mass, 
      				interaction.collision_diameter, 
      				molecule.vibr_frequency[0], 
      				interaction.epsilon, 	
      				molecule.internuclear_distance,
                       		molecule.vibr_energy[0][i] - molecule.vibr_energy[0][i + delta_i], 
      				molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0]);
          } else {

            return k_VT_SSH(	T, 
      				i, 
      				interaction.collision_mass, 
      				interaction.collision_diameter, 
      				molecule.vibr_frequency[0], 
      				interaction.epsilon, 
      				molecule.internuclear_distance);
          }
        } else if (delta_i == 1) {
          if (molecule.anharmonic_spectrum) {

            return k_VT_SSH(	T, 
      				i + delta_i, 
      				interaction.collision_mass, 
      				interaction.collision_diameter, 	
      				molecule.vibr_frequency[0], 
      				interaction.epsilon, 
      				molecule.internuclear_distance,
               	        	molecule.vibr_energy[0][i + delta_i] - molecule.vibr_energy[0][i], 
      				molecule.vibr_energy[0][1] - molecule.vibr_energy[0][0]) * k_bf_VT(T, molecule, i + delta_i, -delta_i, 0);
          } else {

            return k_VT_SSH(	T, 	
      				i + delta_i, 
      				interaction.collision_mass, 
      				interaction.collision_diameter, 
      				molecule.vibr_frequency[0], 
      				interaction.epsilon, 
      				molecule.internuclear_distance) * k_bf_VT(T, molecule, i + delta_i, -delta_i, 0);
          } 
        } else {
          std::string error_string = "Only one-quantum transitions allowed for the SSH model";
          throw kappa::IncorrectValueException(error_string.c_str());
        }
        break;
      }

      // Billing
      case models_k_vt::model_k_vt_billing: {

        if (e != 0) {
         std::string error_string = "No Billing approximation for electronically excited molecules";
         throw kappa::ModelParameterException(error_string.c_str());
        }
        if ((interaction.particle1_name == "N2") && (interaction.particle2_name == "N2")) {
          if (delta_i == -1) {
            return k_VT_Billing_N2N2(T, i);
          } else if (delta_i == 1) {
            return k_VT_Billing_N2N2(T, i + delta_i) * k_bf_VT(T, molecule, i + delta_i, -delta_i, 0);
          } else {
            std::string error_string = "No Billing approximation for N2+N2 non-monoquantum VT transition";  // TODO: add reverse rate
            throw kappa::ModelParameterException(error_string.c_str());
          }
        } else if (((interaction.particle1_name == "N2") && (interaction.particle2_name == "N")) || 
                   ((interaction.particle1_name == "N") && (interaction.particle2_name == "N2"))) {
          if (delta_i < 0) {
            return k_VT_Billing_N2N(T, i, delta_i);
          } else {
            return k_VT_Billing_N2N(T, i + delta_i, -delta_i) * k_bf_VT(T, molecule, i + delta_i, -delta_i, 0);
          }
        } else  {
          std::string error_string = "No Billing approximation for VT rate for the following species: " + interaction.particle1_name + ", " 
          											      + interaction.particle2_name;
          throw kappa::ModelParameterException(error_string.c_str());
        }
      } 

      // ESA phys4entry 
      case models_k_vt::model_k_vt_phys4entry: {
        if (e != 0) {
          std::string error_string = "No Phys4Entry approximation for electronically excited molecules";
          throw kappa::ModelParameterException(error_string.c_str());
        }
        if (((interaction.particle1_name == "N2") && (interaction.particle2_name == "N")) || 
            ((interaction.particle1_name == "N") && (interaction.particle2_name == "N2"))) {
          if (delta_i < 0 ) {
            return k_VT_N2N_p4e(T, (int)convert_vibr_ladder_N2(molecule.vibr_energy[e][i]), delta_i);
          } else {
            return k_VT_N2N_p4e(T, (int)convert_vibr_ladder_N2(molecule.vibr_energy[e][i]), -delta_i) * k_bf_VT(T, molecule, i + delta_i, -delta_i, 0);
          }
          // else {
          //   std::string error_string = "No Billing approximation for N2+N2 non-monoquantum VT transition";  // TODO: add reverse rate
          //   throw kappa::ModelParameterException(error_string.c_str());
          // }
        } else if (((interaction.particle1_name == "O2") && (interaction.particle2_name == "O")) || 
                   ((interaction.particle1_name == "O") && (interaction.particle2_name == "O2"))) {
          if (delta_i < 0) {
            return k_VT_O2O_p4e(T, (int)convert_vibr_ladder_O2(molecule.vibr_energy[e][i]), delta_i);
          } else {
            return k_VT_O2O_p4e(T, (int)convert_vibr_ladder_O2(molecule.vibr_energy[e][i]), -delta_i) * k_bf_VT(T, molecule, i + delta_i, -delta_i, 0);
          }
        } else  {
          std::string error_string = "No Phys4Entry approximation for VT rate for the following species: " + interaction.particle1_name + ", " 
           											         + interaction.particle2_name;
          throw kappa::ModelParameterException(error_string.c_str());
        }
      }
      default:
        std::string error_string = "Unknown choice of interaction potential";
        throw kappa::ModelParameterException(error_string.c_str());
    }
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Simplified inteface (e=0) version for VT
  double kappa::Approximation::k_VT(	double T, 
					kappa::Molecule const &molecule, 
					kappa::Interaction const &interaction, 
					int i, 
					int delta_i, 
					kappa::models_k_vt model) {

    return k_VT(T, molecule, interaction, i, delta_i, 0, model);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::k_exch(	double T, 
					kappa::Molecule const &molecule, 
					kappa::Atom const &atom, 
					kappa::Interaction const &interaction,
                                    	int i, 
					int e, 
					int num_electron_levels, 
					kappa::models_k_exch model) {

    #ifdef KAPPA_STRICT_CHECKS
    if (T <= 0) {
      std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (i < 0) {
      std::string error_string = "Negative vibrational level: i=" + std::to_string(i);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if ((num_electron_levels < -1) || (num_electron_levels == 0)) {
      std::string error_string = "Negative or zero value of number of electronic levels (can only be -1 or positive): num_electron_levels=" 
                               + std::to_string(num_electron_levels);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (num_electron_levels > molecule.num_electron_levels) {
      std::string error_string = "Number of electronic levels passed to function larger than number of levels: num_electron_levels:"
                               + std::to_string(num_electron_levels) + ", number of levels" + std::to_string(molecule.num_electron_levels);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e < 0) {
      std::string error_string = "Negative electronic level: e=" + std::to_string(e);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (((interaction.particle1_name != molecule.name) || (interaction.particle2_name != atom.name)) && 
       ((interaction.particle1_name != atom.name) || (interaction.particle2_name != molecule.name))) {
      std::string error_string = "Interaction is for a different set of particles, the interaction is for: " 
      			       + interaction.particle1_name + " + " + interaction.particle2_name;
      error_string += "; the particles passed to the function are " + molecule.name + ", " + atom.name;
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    switch (model) {

      // model_k_exch_arrh_scanlon
      case models_k_exch::model_k_exch_arrh_scanlon:

        if (((molecule.name=="N2") && (atom.name=="O")) || ((molecule.name=="O2") && (atom.name=="N")) || 
	   ((molecule.name=="NO") && (atom.name=="O")) || ((molecule.name=="NO") && (atom.name=="N"))) {
          return k_Arrhenius(T, interaction["exch,Arrh_A,Scanlon"], interaction["exch,Arrh_n,Scanlon"], interaction["exch,Ea,Scanlon"]);
        } else {
          std::string error_string = "No Scanlon data for exchange reaction for the following species: " + molecule.name + ", " + atom.name;
          throw kappa::ModelParameterException(error_string.c_str());
        }
        break;

      // model_k_exch_arrh_park
      case models_k_exch::model_k_exch_arrh_park:

        if (((molecule.name=="N2") && (atom.name=="O")) || ((molecule.name=="NO") && (atom.name=="O"))) {
	  return k_Arrhenius(T, interaction["exch,Arrh_A,Park"], interaction["exch,Arrh_n,Park"], interaction["exch,Ed,Park"]);
        } else {
    	  std::string error_string = "No Park data for exchange reaction for the following species: " + molecule.name + ", " + atom.name;
          throw kappa::ModelParameterException(error_string.c_str());
        }
        break;

      // model_k_exch_warnatz
      case models_k_exch::model_k_exch_warnatz:

        if (((molecule.name=="N2") && (atom.name=="O")) || ((molecule.name=="O2") && (atom.name=="N"))) {
          return k_exch_WRFP(T, molecule.vibr_energy[e][i], interaction["exch,Ea,WRFP"], 1., 1., interaction["exch,A,Warnatz"] * (i+1), interaction["exch,n,Warnatz"]);
        } else {
          std::string error_string = "No Warnatz data for exchange reaction for the following species: " + molecule.name + ", " + atom.name;
          throw kappa::ModelParameterException(error_string.c_str());
        }
        break;

      // model_k_exch_rf (Rusanov-Friedman)
      case models_k_exch::model_k_exch_rf:

        if (((molecule.name=="N2") && (atom.name=="O")) || ((molecule.name=="O2") && (atom.name=="N"))) {
          return k_exch_WRFP(T, molecule.vibr_energy[e][i], interaction["exch,Ea,WRFP"], interaction["exch,alpha,RF"], 1., interaction["exch,A,RFP"], interaction["exch,n,RF"]);
        } else {
          std::string error_string = "No Rusanov-Friedman data for exchange reaction for the following species: " + molecule.name + ", " + atom.name;
          throw kappa::ModelParameterException(error_string.c_str());
        }
        break;

      // model_k_exch_polak
      case models_k_exch::model_k_exch_polak:

        if (((molecule.name=="N2") && (atom.name=="O")) || ((molecule.name=="O2") && (atom.name=="N"))) {
          return k_exch_WRFP(T, 
			     molecule.vibr_energy[e][i], 
			     interaction["exch,Ea,WRFP"], 
			     interaction["exch,alpha,Polak"], 
			     interaction["exch,beta,Polak"],
		   	     interaction["exch,A,RFP"], 
			     interaction["exch,n,Polak"]);
        } else {
          std::string error_string = "No Polak data for exchange reaction for the following species: " + molecule.name + ", " + atom.name;
          throw kappa::ModelParameterException(error_string.c_str());
        }
        break;
      default:
        std::string error_string = "Unknown choice of exchange reaction model";
        throw kappa::ModelParameterException(error_string.c_str());
    }
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // simplest interface (e=0 and only 1 electronic state) for the relaxation of chemical reactions
  double kappa::Approximation::k_exch(	double T, 	
					kappa::Molecule const &molecule, 
					kappa::Atom const &atom, 
					kappa::Interaction const &interaction, 
					int i, 
					kappa::models_k_exch model) {

    return k_exch(T, molecule, atom, interaction, i, 0, 1, model);
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::k_exch(	double T, 
					kappa::Molecule const &molecule, 
					kappa::Atom const &atom, 
					kappa::Molecule const &molecule_prod, 
					kappa::Interaction const &interaction,
                                    	int i, 
					int k, 
					int e, 
					int num_electron_levels, 
					kappa::models_k_exch model) {

    #ifdef KAPPA_STRICT_CHECKS
    if (T <= 0) {
      std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (i < 0) {
      std::string error_string = "Negative vibrational level: i=" + std::to_string(i);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (k < 0) {
      std::string error_string = "Negative vibrational level: k=" + std::to_string(k);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if ((num_electron_levels < -1) || (num_electron_levels == 0)) {
      std::string error_string = "Negative or zero value of number of electronic levels (can only be -1 or positive): num_electron_levels=" + std::to_string(num_electron_levels);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (num_electron_levels > molecule.num_electron_levels) {
      std::string error_string = "Number of electron levels passed to function larger than number of levels: num_electron_levels:"
                              + std::to_string(num_electron_levels) + ", number of levels" + std::to_string(molecule.num_electron_levels);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e < 0) {
      std::string error_string = "Negative electronic level: e=" + std::to_string(e);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (((interaction.particle1_name != molecule.name) || (interaction.particle2_name != atom.name)) && 
       ((interaction.particle1_name != atom.name) || (interaction.particle2_name != molecule.name))) {
      std::string error_string = "Interaction is for a different set of particles, the interaction is for: " 
			      + interaction.particle1_name + " + " + interaction.particle2_name;
      error_string += "; the particles passed to the function are " + molecule.name + ", " + atom.name;
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    if (num_electron_levels == -1) {
      num_electron_levels = molecule.num_electron_levels;
    }

    switch (model) {

      // model_k_exch_maliat_D6k_arrh_scanlon
      case models_k_exch::model_k_exch_maliat_D6k_arrh_scanlon:

        if (((molecule.name=="N2") && (atom.name=="O")) || ((molecule.name=="O2") && (atom.name=="N")) || 
	    ((molecule.name=="NO") && (atom.name=="O")) ||((molecule.name=="NO") && (atom.name=="N"))) {
          return k_Arrhenius(	T, 
				interaction["exch,Arrh_A,Scanlon"], 
				interaction["exch,Arrh_n,Scanlon"], 
				interaction["exch,Ea,Scanlon"]) *
              	 C_aliat(	T, 
				molecule.electron_energy, 
				molecule.statistical_weight, 
				num_electron_levels, 
				molecule.vibr_energy, 
               			molecule.num_vibr_levels, 
				molecule_prod.vibr_energy[0][k], 
				interaction["exch,Ea,Scanlon"], 
				i, 
				e, 
				molecule.diss_energy[e] / (6 * K_CONST_K));
        } else {
          std::string error_string = "No Scanlon data for exchange reaction for the following species: " + molecule.name + ", " + atom.name;
          throw kappa::ModelParameterException(error_string.c_str());
        }
        break;

      // model_k_exch_maliat_3T_arrh_scanlon
      case models_k_exch::model_k_exch_maliat_3T_arrh_scanlon:

        if (((molecule.name=="N2") && (atom.name=="O")) || ((molecule.name=="O2") && (atom.name=="N")) || 
	    ((molecule.name=="NO") && (atom.name=="O")) ||((molecule.name=="NO") && (atom.name=="N"))) {
          return k_Arrhenius(	T, 		
				interaction["exch,Arrh_A,Scanlon"], 
				interaction["exch,Arrh_n,Scanlon"], 
				interaction["exch,Ea,Scanlon"]) * 
               	 C_aliat(	T, 
				molecule.electron_energy, 
				molecule.statistical_weight, 
				num_electron_levels,
                         	molecule.vibr_energy, 
				molecule.num_vibr_levels, 
				molecule_prod.vibr_energy[0][k], 
				interaction["exch,Ea,Scanlon"], 
				i, 
				e, 
				3 * T);
        } else {
          std::string error_string = "No Scanlon data for exchange reaction for the following species: " + molecule.name + ", " + atom.name;
          throw kappa::ModelParameterException(error_string.c_str());
        }
        break;

      // model_k_exch_maliat_infty_arrh_scanlon
      case models_k_exch::model_k_exch_maliat_infty_arrh_scanlon:

        if (((molecule.name=="N2") && (atom.name=="O")) || ((molecule.name=="O2") && (atom.name=="N")) || 
	    ((molecule.name=="NO") && (atom.name=="O")) ||((molecule.name=="NO") && (atom.name=="N"))) {
          return k_Arrhenius(	T, 
				interaction["exch,Arrh_A,Scanlon"], 
				interaction["exch,Arrh_n,Scanlon"], 
				interaction["exch,Ea,Scanlon"]) * 
             	 C_aliat(	T, 
				molecule.electron_energy, 
				molecule.statistical_weight, 
				num_electron_levels,
                         	molecule.vibr_energy, 
				molecule.num_vibr_levels, 
				molecule_prod.vibr_energy[0][k], 
				interaction["exch,Ea,Scanlon"], 
				i, 
				e);
        } else {
          std::string error_string = "No Scanlon data for exchange reaction for the following species: " + molecule.name + ", " + atom.name;
          throw kappa::ModelParameterException(error_string.c_str());
        }
        break;

      // model_k_exch_maliat_D6k_arrh_park   
      case models_k_exch::model_k_exch_maliat_D6k_arrh_park:

        if (((molecule.name=="N2") && (atom.name=="O")) || ((molecule.name=="NO") && (atom.name=="O"))) {
          return k_Arrhenius(	T, 
				interaction["exch,Arrh_A,Park"], 	
				interaction["exch,Arrh_n,Park"], 
				interaction["exch,Ed,Park"]) * 
          	 C_aliat(	T, 
				molecule.electron_energy, 
				molecule.statistical_weight, 
				num_electron_levels,
                         	molecule.vibr_energy, 
				molecule.num_vibr_levels, 
				molecule_prod.vibr_energy[0][k], 
				interaction["exch,Ea,Scanlon"], 
				i, 
				e, 
				molecule.diss_energy[e] / (6 * K_CONST_K));
        } else {
          std::string error_string = "No Park data for exchange reaction for the following species: " + molecule.name + ", " + atom.name;
          throw kappa::ModelParameterException(error_string.c_str());
        }
        break;

      // model_k_exch_maliat_3T_arrh_park
      case models_k_exch::model_k_exch_maliat_3T_arrh_park:

        if (((molecule.name=="N2") && (atom.name=="O")) || ((molecule.name=="NO") && (atom.name=="O"))) {
          return k_Arrhenius(	T, 
				interaction["exch,Arrh_A,Park"], 
				interaction["exch,Arrh_n,Park"], 	
				interaction["exch,Ed,Park"]) * 
           	 C_aliat(	T, 
				molecule.electron_energy, 
				molecule.statistical_weight, 
				num_electron_levels,
                         	molecule.vibr_energy, 
				molecule.num_vibr_levels, 
				molecule_prod.vibr_energy[0][k], 
				interaction["exch,Ea,Scanlon"], 
				i, 
				e, 
				3 * T);
        } else {
          std::string error_string = "No Park data for exchange reaction for the following species: " + molecule.name + ", " + atom.name;
          throw kappa::ModelParameterException(error_string.c_str());
        }
        break;

      // model_k_exch_maliat_infty_arrh_park
      case models_k_exch::model_k_exch_maliat_infty_arrh_park:

        if (((molecule.name=="N2") && (atom.name=="O")) || ((molecule.name=="NO") && (atom.name=="O"))) {
          return k_Arrhenius(	T, 
				interaction["exch,Arrh_A,Park"], 
				interaction["exch,Arrh_n,Park"], 
				interaction["exch,Ed,Park"]) * 
              	 C_aliat(	T, 
				molecule.electron_energy, 
				molecule.statistical_weight, 
				num_electron_levels,
                         	molecule.vibr_energy, 
				molecule.num_vibr_levels, 
				molecule_prod.vibr_energy[0][k], 
				interaction["exch,Ea,Scanlon"], 
				i, 
				e);
        } else {
          std::string error_string = "No Park data for exchange reaction for the following species: " + molecule.name + ", " + atom.name;
          throw kappa::ModelParameterException(error_string.c_str());
        }
        break;

      default:
        return k_exch(T, molecule, atom, interaction, i, e, num_electron_levels, model);
    }
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Релаксационные члены, описывающие (бимолекулярные) обменные реакции
  double kappa::Approximation::k_exch(	double T, 
					kappa::Molecule const &molecule, 
					kappa::Atom const &atom, 
					kappa::Molecule const &molecule_prod, 
					kappa::Interaction const &interaction,
                                    	int i, 
					int k, 
					kappa::models_k_exch model) {

    return k_exch(T, molecule, atom, molecule_prod, interaction, i, k, 0, 1, model);
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::integral_diss(double T, int degree, kappa::Molecule const &molecule, kappa::Interaction const &interaction, int i, int e, kappa::models_cs_diss model) {

    #ifdef KAPPA_STRICT_CHECKS
    if (T <= 0) {
      std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (degree < 0) {
      std::string error_string = "Negative degree of integral: degree=" + std::to_string(degree);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e < 0) {
      std::string error_string = "Negative electronic level: e=" + std::to_string(e);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (i < 0) {
      std::string error_string = "Negative vibrational level: i=";
      error_string += std::to_string(i);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if ((interaction.particle1_name != molecule.name) && (interaction.particle2_name != molecule.name)) {
      std::string error_string = "Interaction is not for the molecule passed for the function, the interaction is for: " + interaction.particle1_name + " + " + interaction.particle2_name;
      error_string += "; the molecule passed to the function is " + molecule.name;
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    switch (model) {
      case models_cs_diss::model_cs_diss_rs_thresh_cmass_vibr:
	return integral_diss_RS(T, degree, interaction.collision_mass, interaction.collision_diameter, molecule.diss_energy[e], molecule.vibr_energy[e][i], true);
	break;
      case models_cs_diss::model_cs_diss_rs_thresh_vibr:
	return integral_diss_RS(T, degree, interaction.collision_mass, interaction.collision_diameter, molecule.diss_energy[e], molecule.vibr_energy[e][i], false);
	break;
      case models_cs_diss::model_cs_diss_rs_thresh_cmass:
	return integral_diss_RS(T, degree, interaction.collision_mass, interaction.collision_diameter, molecule.diss_energy[e], 0, true);
	break;
      case models_cs_diss::model_cs_diss_rs_thresh:
	return integral_diss_RS(T, degree, interaction.collision_mass, interaction.collision_diameter, molecule.diss_energy[e], 0, false);
	break;
      case models_cs_diss::model_cs_diss_vss_thresh_cmass_vibr:
        if (interaction.vss_data) {
          return integral_diss_VSS(T, degree, interaction.collision_mass, interaction.vss_c_cs, interaction.vss_omega, molecule.diss_energy[e], molecule.vibr_energy[e][i], true);
        } else {
          std::string error_string = "No VSS data found for " + interaction.particle1_name + " + " + interaction.particle2_name + " interaction";
          throw kappa::DataNotFoundException(error_string.c_str());
        }
	break;
      case models_cs_diss::model_cs_diss_vss_thresh_vibr:
        if (interaction.vss_data) {
          return integral_diss_VSS(T, degree, interaction.collision_mass, interaction.vss_c_cs, interaction.vss_omega, molecule.diss_energy[e], molecule.vibr_energy[e][i], false);
        } else {
          std::string error_string = "No VSS data found for " + interaction.particle1_name + " + " + interaction.particle2_name + " interaction";
          throw kappa::DataNotFoundException(error_string.c_str());
        }
	break;
      case models_cs_diss::model_cs_diss_vss_thresh_cmass:
        if (interaction.vss_data) {
          return integral_diss_VSS(T, degree, interaction.collision_mass, interaction.vss_c_cs, interaction.vss_omega, molecule.diss_energy[e], 0, true);
        } else {
          std::string error_string = "No VSS data found for " + interaction.particle1_name + " + " + interaction.particle2_name + " interaction";
          throw kappa::DataNotFoundException(error_string.c_str());
        }
	break;
      case models_cs_diss::model_cs_diss_vss_thresh:
        if (interaction.vss_data) {
          return integral_diss_VSS(T, degree, interaction.collision_mass, interaction.vss_c_cs, interaction.vss_omega, molecule.diss_energy[e], 0, false);
        } else {
          std::string error_string = "No VSS data found for " + interaction.particle1_name + " + " + interaction.particle2_name + " interaction";
          throw kappa::DataNotFoundException(error_string.c_str());
        }
	break;
      case models_cs_diss::model_cs_diss_ilt:
        if (e==0) {
          if (((interaction.particle1_name == "N2") && (interaction.particle2_name == "N")) || ((interaction.particle1_name == "N") && (interaction.particle2_name == "N2"))) {
            return integral_diss_ILT_N2N(T, degree, interaction.collision_mass, molecule.vibr_energy[e][i], convert_vibr_ladder_N2(molecule.vibr_energy[e][i]));
          } else if (((interaction.particle1_name == "O2") && (interaction.particle2_name == "O")) || ((interaction.particle1_name == "O") && (interaction.particle2_name == "O2"))) {
            return integral_diss_ILT_O2O(T, degree, interaction.collision_mass, molecule.vibr_energy[e][i], convert_vibr_ladder_O2(molecule.vibr_energy[e][i]));
          } else {
            std::string error_string = "No ILT model available for dissociation for the interaction: " + interaction.particle1_name + "+" + interaction.particle2_name;
            throw kappa::ModelParameterException(error_string.c_str());
          }
        } else {
          std::string error_string = "No ILT model available for dissociation from excited electronic states";
          throw kappa::ModelParameterException(error_string.c_str());
        }
        break;
      default:
        std::string error_string = "Unknown choice of dissociation model";
        throw kappa::ModelParameterException(error_string.c_str());
    }
  }
 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::Approximation::integral_diss(	double T, 
	     					int degree, 
     						kappa::Molecule const &molecule, 
     						kappa::Interaction const &interaction, 
     						int i, 
     						kappa::models_cs_diss model) {

    return integral_diss(T, degree, molecule, interaction, i, 0, model);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Computation of the dissociation rate for a selected "particle" interaction and dissociation model
  double kappa::Approximation::k_diss(	double T, 
					kappa::Molecule const &molecule, 
					kappa::Interaction const &interaction, 
					int i, 
					int e, 
					int num_electron_levels, 
					kappa::models_k_diss model) {

    #ifdef KAPPA_STRICT_CHECKS
    if (T <= 0) {
      std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e < 0) {
      std::string error_string = "Negative electronic level: e=" + std::to_string(e);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (i < 0) {
      std::string error_string = "Negative vibrational level: i=";
      error_string += std::to_string(i);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if ((interaction.particle1_name != molecule.name) && (interaction.particle2_name != molecule.name)) {
      std::string error_string = "Interaction is not for the molecule passed for the function, the interaction is for: " + interaction.particle1_name + 
														   " + " + interaction.particle2_name;
      error_string += "; the molecule passed to the function is " + molecule.name;
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif


    switch (model) {
      case models_k_diss::model_k_diss_rs_thresh_cmass_vibr:
        return k_diss_RS(T, interaction.collision_mass, interaction.collision_diameter, molecule.diss_energy[e], molecule.vibr_energy[e][i], true);
        break;
      case models_k_diss::model_k_diss_rs_thresh_vibr:
        return k_diss_RS(T, interaction.collision_mass, interaction.collision_diameter, molecule.diss_energy[e], molecule.vibr_energy[e][i], false);
        break;
      case models_k_diss::model_k_diss_rs_thresh_cmass:
        return k_diss_RS(T, interaction.collision_mass, interaction.collision_diameter, molecule.diss_energy[e], 0, true);
        break;
      case models_k_diss::model_k_diss_rs_thresh:
        return k_diss_RS(T, interaction.collision_mass, interaction.collision_diameter, molecule.diss_energy[e], 0, false);
        break;
      case models_k_diss::model_k_diss_vss_thresh_cmass_vibr:
        if (interaction.vss_data) {
          return k_diss_VSS(T, interaction.collision_mass, interaction.vss_c_cs, interaction.vss_omega, molecule.diss_energy[e], molecule.vibr_energy[e][i], true);
        } else {
          std::string error_string = "No VSS data found for " + interaction.particle1_name + " + " + interaction.particle2_name + " interaction";
          throw kappa::DataNotFoundException(error_string.c_str());
        }
        break;
      case models_k_diss::model_k_diss_vss_thresh_vibr:
        if (interaction.vss_data) {
          return k_diss_VSS(T, interaction.collision_mass, interaction.vss_c_cs, interaction.vss_omega, molecule.diss_energy[e], molecule.vibr_energy[e][i], false);
        } else {
          std::string error_string = "No VSS data found for " + interaction.particle1_name + " + " + interaction.particle2_name + " interaction";
          throw kappa::DataNotFoundException(error_string.c_str());
        }
        break;
      case models_k_diss::model_k_diss_vss_thresh_cmass:
        if (interaction.vss_data) {
          return k_diss_VSS(T, interaction.collision_mass, interaction.vss_c_cs, interaction.vss_omega, molecule.diss_energy[e], 0, true);
        } else {
          std::string error_string = "No VSS data found for " + interaction.particle1_name + " + " + interaction.particle2_name + " interaction";
          throw kappa::DataNotFoundException(error_string.c_str());
        }
        break;
      case models_k_diss::model_k_diss_vss_thresh:
        if (interaction.vss_data) {
          return k_diss_VSS(T, interaction.collision_mass, interaction.vss_c_cs, interaction.vss_omega, molecule.diss_energy[e], 0, false);
        } else {
          std::string error_string = "No VSS data found for " + interaction.particle1_name + " + " + interaction.particle2_name + " interaction";
          throw kappa::DataNotFoundException(error_string.c_str());
        }
        break;
      case models_k_diss::model_k_diss_ilt:
        if (e==0) {
          if (((interaction.particle1_name == "N2") && (interaction.particle2_name == "N")) || 
              ((interaction.particle1_name == "N" ) && (interaction.particle2_name == "N2"))) {

            return k_diss_ILT_N2N(T, interaction.collision_mass, molecule.vibr_energy[e][i], convert_vibr_ladder_N2(molecule.vibr_energy[e][i]));

          } else if (((interaction.particle1_name == "O2") && (interaction.particle2_name == "O")) || 
                     ((interaction.particle1_name == "O" ) && (interaction.particle2_name == "O2"))) {

            return k_diss_ILT_O2O(T, interaction.collision_mass, molecule.vibr_energy[e][i], convert_vibr_ladder_O2(molecule.vibr_energy[e][i]));
          } else {
            std::string error_string = "No ILT model available for dissociation for the interaction: " + interaction.particle1_name + "+" + interaction.particle2_name;
            throw kappa::ModelParameterException(error_string.c_str());
          }
        } else {
          std::string error_string = "No ILT model available for dissociation from excited electronic states";
          throw kappa::ModelParameterException(error_string.c_str());
        }
        break;
      case models_k_diss::model_k_diss_arrh_scanlon:
        return k_Arrhenius(	T, 
				interaction["diss," + molecule.name + ",Arrh_A,Scanlon"], 
				interaction["diss," + molecule.name + ",Arrh_n,Scanlon"], 
				interaction["diss," + molecule.name + ",Ea,Scanlon"]);
        break;
      case models_k_diss::model_k_diss_arrh_park:
        return k_Arrhenius(	T, 
				interaction["diss," + molecule.name + ",Arrh_A,Park"], 
				interaction["diss," + molecule.name + ",Arrh_n,Park"], 
				interaction["diss," + molecule.name + ",Ed,Park"]);
        break;
      case models_k_diss::model_k_diss_tm_D6k_arrh_scanlon:
        if (num_electron_levels == 1) {
          return k_Arrhenius(	T, 
				interaction["diss," + molecule.name + ",Arrh_A,Scanlon"], 
				interaction["diss," + molecule.name + ",Arrh_n,Scanlon"],
				interaction["diss," + molecule.name + ",Ea,Scanlon"]) * 
            		 	Z_diss(T, molecule.diss_energy[e] / (6 * K_CONST_K), molecule.vibr_energy[e], i);
        } else if (num_electron_levels == -1) {
          return k_Arrhenius(	T, 
				interaction["diss," + molecule.name + ",Arrh_A,Scanlon"], 
				interaction["diss," + molecule.name + ",Arrh_n,Scanlon"], 
				interaction["diss," + molecule.name + ",Ea,Scanlon"]) * 
            			Z_diss(	T, 
					molecule.diss_energy[e] / (6 * K_CONST_K), 
					molecule.electron_energy, 
					molecule.statistical_weight, 
					molecule.num_electron_levels, 
					molecule.vibr_energy, 
					i, e);
        } else {
          return k_Arrhenius(	T, 
				interaction["diss," + molecule.name + ",Arrh_A,Scanlon"], 
				interaction["diss," + molecule.name + ",Arrh_n,Scanlon"], 
				interaction["diss," + molecule.name + ",Ea,Scanlon"]) * 
 	    			Z_diss(	T, 
					molecule.diss_energy[e] / (6 * K_CONST_K), 
					molecule.electron_energy, 
					molecule.statistical_weight, 
					num_electron_levels, 
					molecule.vibr_energy, 	
					i, e);
        }
        break;
      case models_k_diss::model_k_diss_tm_3T_arrh_scanlon:
        if (num_electron_levels == 1) {
          return k_Arrhenius(	T, 
				interaction["diss," + molecule.name + ",Arrh_A,Scanlon"], 
				interaction["diss," + molecule.name + ",Arrh_n,Scanlon"], 
				interaction["diss," + molecule.name + ",Ea,Scanlon"]) * 
      				Z_diss(T, 3 * T, molecule.vibr_energy[e], i);

        } else if (num_electron_levels == -1) {
          return k_Arrhenius(	T, 
				interaction["diss," + molecule.name + ",Arrh_A,Scanlon"], 
				interaction["diss," + molecule.name + ",Arrh_n,Scanlon"], 
				interaction["diss," + molecule.name + ",Ea,Scanlon"]) * 
            			Z_diss(	T, 
					3 * T, 
					molecule.electron_energy,
					molecule.statistical_weight, 
					molecule.num_electron_levels, 
					molecule.vibr_energy, 
					i, e);
        } else {
          return k_Arrhenius(	T, 
				interaction["diss," + molecule.name + ",Arrh_A,Scanlon"], 
				interaction["diss," + molecule.name + ",Arrh_n,Scanlon"], 
				interaction["diss," + molecule.name + ",Ea,Scanlon"]) * 
	    			Z_diss(	T, 
					3 * T, 
					molecule.electron_energy, 
					molecule.statistical_weight, 
					num_electron_levels, molecule.vibr_energy, 
					i, e);
        }
        break;
      case models_k_diss::model_k_diss_tm_infty_arrh_scanlon:
        if (num_electron_levels == 1) {
          return k_Arrhenius(	T, 
				interaction["diss," + molecule.name + ",Arrh_A,Scanlon"], 
				interaction["diss," + molecule.name + ",Arrh_n,Scanlon"], 
				interaction["diss," + molecule.name + ",Ea,Scanlon"]) * 
            			Z_diss(T, molecule.vibr_energy[e], molecule.num_vibr_levels[e], i);

        } else if (num_electron_levels == -1) {
          return k_Arrhenius(	T, 
				interaction["diss," + molecule.name + ",Arrh_A,Scanlon"], 
				interaction["diss," + molecule.name + ",Arrh_n,Scanlon"], 
				interaction["diss," + molecule.name + ",Ea,Scanlon"]) * 
            			Z_diss(	T, 
					molecule.electron_energy, 
					molecule.statistical_weight, 
					molecule.num_electron_levels, 
					molecule.vibr_energy, 
					molecule.num_vibr_levels, 
					i, e);
        } else {
          return k_Arrhenius(	T, 
				interaction["diss," + molecule.name + ",Arrh_A,Scanlon"], 
				interaction["diss," + molecule.name + ",Arrh_n,Scanlon"], 
				interaction["diss," + molecule.name + ",Ea,Scanlon"]) * 
            			Z_diss(T, molecule.electron_energy, molecule.statistical_weight, num_electron_levels, molecule.vibr_energy, molecule.num_vibr_levels, i, e);
        } 
        break;
      case models_k_diss::model_k_diss_tm_D6k_arrh_park:
        if (num_electron_levels == 1) {
          return k_Arrhenius(	T, 
				interaction["diss," + molecule.name + ",Arrh_A,Park"], 
				interaction["diss," + molecule.name + ",Arrh_n,Park"], 
				interaction["diss," + molecule.name + ",Ed,Park"]) * 
            			Z_diss(T, molecule.diss_energy[e] / (6 * K_CONST_K), molecule.vibr_energy[e], i);

        } else if (num_electron_levels == -1) {
          return k_Arrhenius(	T, 
				interaction["diss," + molecule.name + ",Arrh_A,Park"], 
				interaction["diss," + molecule.name + ",Arrh_n,Park"], 
				interaction["diss," + molecule.name + ",Ed,Park"]) * 
            			Z_diss(	T, 
					molecule.diss_energy[e] / (6 * K_CONST_K), 
					molecule.electron_energy, 
					molecule.statistical_weight, 
					molecule.num_electron_levels, 
					molecule.vibr_energy, 
					i, e);
        } else {
          return k_Arrhenius(	T, 
				interaction["diss," + molecule.name + ",Arrh_A,Park"], 
				interaction["diss," + molecule.name + ",Arrh_n,Park"], 
				interaction["diss," + molecule.name + ",Ed,Park"]) * 
            			Z_diss(	T, 
					molecule.diss_energy[e] / (6 * K_CONST_K), 
					molecule.electron_energy, 
					molecule.statistical_weight, 
					num_electron_levels, 
					molecule.vibr_energy, 
					i, e);
        }
        break;
      case models_k_diss::model_k_diss_tm_3T_arrh_park:
        if (num_electron_levels == 1) {
          return k_Arrhenius(	T, 
				interaction["diss," + molecule.name + ",Arrh_A,Park"], 
				interaction["diss," + molecule.name + ",Arrh_n,Park"], 
				interaction["diss," + molecule.name + ",Ed,Park"]) * 
            			Z_diss(T, 3 * T, molecule.vibr_energy[e], i);

         } else if (num_electron_levels == -1) {
           return k_Arrhenius(	T, 
				interaction["diss," + molecule.name + ",Arrh_A,Park"], 
				interaction["diss," + molecule.name + ",Arrh_n,Park"], 
				interaction["diss," + molecule.name + ",Ed,Park"]) * 
            			Z_diss(T, 3 * T, molecule.electron_energy, molecule.statistical_weight, molecule.num_electron_levels, molecule.vibr_energy, i, e);
        } else {
          return k_Arrhenius( 	T, 
			    	interaction["diss," + molecule.name + ",Arrh_A,Park"], 
				interaction["diss," + molecule.name + ",Arrh_n,Park"], 
				interaction["diss," + molecule.name + ",Ed,Park"]) * 
     	    			Z_diss(T, 3 * T, molecule.electron_energy, molecule.statistical_weight, num_electron_levels, molecule.vibr_energy, i, e);
        }
        break;
      case models_k_diss::model_k_diss_tm_infty_arrh_park:
        if (num_electron_levels == 1) {
          return k_Arrhenius(	T, 	
				interaction["diss," + molecule.name + ",Arrh_A,Park"], 
				interaction["diss," + molecule.name + ",Arrh_n,Park"], 
				interaction["diss," + molecule.name + ",Ed,Park"]) * 
            			Z_diss(T, molecule.vibr_energy[e], molecule.num_vibr_levels[e], i);

        } else if (num_electron_levels == -1) {
          return k_Arrhenius(	T, 
				interaction["diss," + molecule.name + ",Arrh_A,Park"], 
				interaction["diss," + molecule.name + ",Arrh_n,Park"], 
				interaction["diss," + molecule.name + ",Ed,Park"]) * 
            			Z_diss(	T, 
					molecule.electron_energy, 
					molecule.statistical_weight, 
					molecule.num_electron_levels, 
					molecule.vibr_energy, 
					molecule.num_vibr_levels, 
					i, e);
        } else {
          return k_Arrhenius(	T, 
				interaction["diss," + molecule.name + ",Arrh_A,Park"], 
				interaction["diss," + molecule.name + ",Arrh_n,Park"], 
				interaction["diss," + molecule.name + ",Ed,Park"]) * 
  	    			Z_diss(	T, 
					molecule.electron_energy, 
					molecule.statistical_weight, 
					num_electron_levels, 
					molecule.vibr_energy, 
					molecule.num_vibr_levels, 
					i, e);
        }
        break;
      case models_k_diss::model_k_diss_phys4entry:
        if (e==0) {
          if (molecule.name == "N2" && (interaction.particle1_name == "N" || interaction.particle2_name == "N")) {
            return diss_rate_N2N_p4e(T, convert_vibr_ladder_N2(molecule.vibr_energy[e][i]));
          } else if (molecule.name == "O2" && (interaction.particle1_name == "O" || interaction.particle2_name == "O")) {
            return diss_rate_O2O_p4e(T, convert_vibr_ladder_O2(molecule.vibr_energy[e][i]));
          } else {
            std::string error_string = "No Phys4Entry data for dissociation of " + molecule.name + " in reaction: " + interaction.particle1_name + "+" 
														    + interaction.particle2_name;
            throw kappa::ModelParameterException(error_string.c_str());
          }
        } else {
           std::string error_string = "Cannot compute dissociation rate using Phys4Entry data for excited electronic states";
           throw kappa::ModelParameterException(error_string.c_str());
        }
        break;
      default:
        std::string error_string = "Unknown choice of dissociation model";
        throw kappa::ModelParameterException(error_string.c_str());
    }
  }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Simplified version of k_diss to be called
  double kappa::Approximation::k_diss(	double T, 
					kappa::Molecule const &molecule, 
					kappa::Interaction const &interaction, 
					int i, 
 					kappa::models_k_diss model) {

    return k_diss(T, molecule, interaction, i, 0, 1, model);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the backward and forward reaction rate coefficients for VT exchange (interface)
  double kappa::Approximation::k_bf_VT(double T, kappa::Molecule const &molecule, int i, int delta_i, int e) {

    return p_k_bf_VT(	T, 
 			molecule.vibr_energy[e][i], 
			molecule.vibr_energy[e][i + delta_i], 
			molecule.rot_energy[e][i], 
			molecule.num_rot_levels[e][i],
                 	molecule.rot_energy[e][i + delta_i], 
			molecule.num_rot_levels[e][i + delta_i], 
			molecule.rot_symmetry, 
			molecule.rigid_rotator);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the backward and forward reaction rate coefficients for VV exchange (interface)
  double kappa::Approximation::k_bf_VV(	double T, 
					kappa::Molecule const &molecule1, 
					kappa::Molecule const &molecule2, 
					int i, 
					int k, 
					int delta_i, 
					int e1, 
					int e2) {

    return p_k_bf_VV(	T, 
			molecule1.vibr_energy[e1][i], 
			molecule1.vibr_energy[e1][i + delta_i], 
			molecule2.vibr_energy[e2][k], 
			molecule2.vibr_energy[e2][k - delta_i],
                   	molecule1.rot_energy[e1][i], 
			molecule1.num_rot_levels[e1][i], 
			molecule1.rot_energy[e1][i + delta_i], 
			molecule1.num_rot_levels[e1][i + delta_i], 
			molecule1.rot_symmetry, 
                   	molecule2.rot_energy[e2][k], 
			molecule2.num_rot_levels[e2][k], 
			molecule1.rot_energy[e2][k - delta_i], 
			molecule2.num_rot_levels[e2][k - delta_i], 
			molecule2.rot_symmetry);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // compute the backward and forward reaction rate coefficients for chemical exchange (interface)
  double kappa::Approximation::k_bf_exch(double T, 
 					 kappa::Molecule const &molecule_before, 
					 kappa::Atom const &atom_before,
                             		 kappa::Molecule const &molecule_after, 
					 kappa::Atom const &atom_after,
                               		 kappa::Interaction const &interaction, 
					 int i_before, 
					 int i_after, 
					 int e_before, 
					 int e_after) {

    return p_k_bf_exch(	T, 
			molecule_before.mass, 
			atom_before.mass, 
			molecule_after.mass, 
			atom_after.mass,
                   	molecule_before.diss_energy[e_before], 
			molecule_after.diss_energy[e_after], 
			molecule_before.vibr_energy[e_before][i_before], 
			molecule_after.vibr_energy[e_after][i_after],
                   	molecule_before.rot_energy[e_before][i_before], 
			molecule_before.num_rot_levels[e_before][i_before], 
			molecule_before.rot_symmetry,
                   	molecule_after.rot_energy[e_after][i_after], 
			molecule_after.num_rot_levels[e_after][i_after], 
			molecule_after.rot_symmetry);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // compute the backward and forward reaction rate coefficients for dissociation-recombination (interface)
  // factor to convert k_diss -> k_rec for different atomic species
  double kappa::Approximation::k_bf_diss(double T, kappa::Molecule const &molecule, kappa::Atom const &atom1, kappa::Atom const &atom2, int i, int e) {

    return p_k_bf_diss(	T, 
			molecule.mass, 
			atom1.mass, 
			atom2.mass, 
			molecule.diss_energy[e], 
			molecule.vibr_energy[e][i], 
			molecule.rot_energy[e][i], 
			molecule.num_rot_levels[e][i], 
			molecule.rot_symmetry);
  } 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // factor to convert k_diss -> k_rec for the same atomic species
  double kappa::Approximation::k_bf_diss(double T, kappa::Molecule const &molecule, kappa::Atom const &atom, int i, int e) {

    return p_k_bf_diss(	T, 
			molecule.mass, 
			atom.mass, 
			atom.mass, 
			molecule.diss_energy[e],
			molecule.vibr_energy[e][i], 
			molecule.rot_energy[e][i],
			molecule.num_rot_levels[e][i], 
			molecule.rot_symmetry);
 } 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Compute the vibrational level populations according to the Boltzmann distribution at temperature T and total number density n
  arma::vec kappa::Approximation::Boltzmann_distribution(double T, double n, kappa::Molecule const &molecule, int e) {

    #ifdef KAPPA_STRICT_CHECKS
    if (T <= 0) {
      std::string error_string = "Non-positive temperature: T=" + std::to_string(T);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (n <= 0) {
      std::string error_string = "Non-positive number density: n=" + std::to_string(n);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    if (e < 0) {
      std::string error_string = "Negative electronic level: e=" + std::to_string(e);
      throw kappa::IncorrectValueException(error_string.c_str());
    }
    #endif

    return n * arma::exp(-molecule.vibr_energy[0] / (K_CONST_K * T)) / p_Z_vibr_eq(T, molecule.vibr_energy[e]);
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // LC -- FIXME: probably wrong 
  int kappa::Approximation::p_max_i(double T, double vibr_energy1, double vibr_frequency, double alpha) {

    int p_i = -1;
    p_i = (vibr_energy1 * T) / (2 * alpha * vibr_frequency * T * K_CONST_H) + 0.5;
    return p_i;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  const double kappa::Approximation::p4e_n2n_vt_bij1[5][5]={{-3.189200e+01, 4.915400e-01,3.303700e-03,-9.362800e-05,-3.517000e+00},
                  					    {-1.413600e+04, 4.269700e+02,-1.129500e+01, 8.816800e-02, 2.445900e+03},
                  					    {7.545700e+03,-2.282000e+03,7.615200e+01,-8.509700e-01, 6.103300e+03},
                 					    {-4.588900e+05, 1.564300e+05,-4.100500e+03,3.705600e+01,-4.393800e+05},
                					    {9.796900e-01,-5.611100e-02,-3.362100e-04,1.052400e-05,3.810500e-01}};

  const double kappa::Approximation::p4e_n2n_vt_bij2[5][5]={{-2.620600e+01,8.255100e-01,6.123000e-04,-9.805500e-05,-7.161100e+00},
  						            {-2.265800e+03,  1.615800e+03,-3.114100e+01, 2.170400e-01, -7.257200e+03},
						            {2.338600e+05,1.571600e+04,-1.909100e+02,6.338200e-01, -1.583300e+05},
  						            {-6.228100e+06,-3.163500e+05,2.992200e+03,-2.897700e+00, 3.819000e+06},
						            {3.333100e-01,-9.797200e-02, 3.756200e-06,1.106000e-05,8.232000e-01}};

  const double kappa::Approximation::p4e_n2n_vt_cij2[5][5]={{-5.611300e+00,  -4.912300e-01,6.910000e-03,-3.185400e-05,4.307800e+00},
  						            {-7.957600e+03,  -8.055500e+02, 1.344900e+01,-8.712100e-02, 6.528700e+03},
           						    {-1.589600e+05,-1.449600e+04,2.410800e+02,-1.532100e+00, 1.225500e+05},
           						    {3.886100e+06,3.431100e+05,-5.689700e+03, 3.633400e+01,-2.949500e+06},
           						    { 6.078700e-01,6.013200e-02,-8.229600e-04, 3.493500e-06,-5.128000e-01}};

  const double kappa::Approximation::p4e_n2n_vt_bij3[5][5]={{-3.371800e+02,  -2.602000e+01,4.249100e-01,-2.700200e-03,2.283700e+02},
           						    {-2.859800e+04,-2.867000e+02, -2.359700e+00,4.173600e-02, 1.094600e+04},
           						    {-2.011800e+06,-2.042300e+05,3.468900e+03,-2.291500e+01, 1.653800e+06},
           						    { 3.856700e+07,4.253000e+06,-7.432600e+04, 5.018800e+02,-3.309400e+07},
           						    { 3.598700e+01, 2.942500e+00,-4.767100e-02,3.007900e-046,-2.605400e+01}};

  const double kappa::Approximation::p4e_n2n_vt_cij3[5][5]={{3.525100e+01, 3.070400e+00,-4.988000e-02,3.179500e-04,-2.673400e+01},
           						    {6.733300e+02,-3.869300e+01, 1.178400e+00,-1.024700e-02,-9.191300e+01},
           						    { 2.211900e+05,2.316200e+04,-3.953300e+02,2.627600e+00,-1.849700e+05},
           						    {-4.030600e+06,-4.639000e+05, 8.147400e+03,-5.540000e+01, 3.549500e+06},
              						    {-4.193900e+00,-3.546000e-01,5.729100e-03,-3.640900e-05,3.126500e+00}};
 
  const double kappa::Approximation::p4e_n2n_vt_bij4[5][5]={{-5.144800e+02,-1.811900e+01, 2.071900e-01,-9.692600e-04, 2.553900e+02},
           						    {3.950300e+04,2.329000e+03,2.329000e+03,2.329000e+03,-2.567600e+04},
           						    {-5.504500e+06,-2.607000e+05,3.365900e+03,-1.814000e+01,3.174700e+06},
           					            {1.380000e+08,6.692000e+06,-8.798600e+04,4.828400e+02,-8.019000e+07},
           						    {4.876200e+01,1.742200e+00,-1.924000e-02,8.480100e-05,-2.511100e+01}};

  const double kappa::Approximation::p4e_n2n_vt_cij4[5][5]={{3.580700e+01,1.501200e+00, -1.845200e-02,9.490200e-05,-1.966800e+01},
           						    {-4.267800e+03,-2.086000e+02,2.816000e+00,-1.594900e-02,2.471900e+03},
           						    {3.556600e+05,1.801800e+04,-2.461400e+02,1.411600e+00,-2.094400e+05},
           						    {-3.901100e+00,-1.633400e-01,5.732900e+03,-3.354800e+01,4.690900e+06},
           						    {4.876200e+01,1.742200e+00,2.005500e-03,-1.027200e-05,2.139000e+00}};

  const double kappa::Approximation::p4e_n2n_vt_bij5[5][5]={{-1.913900e+03,-5.248100e+01,4.885700e-01,-1.879900e-03,9.005900e+02},
           						    {5.118300e+04,1.408200e+03,-9.481500e+00,3.505700e-03,-2.518000e+04},
           						    {2.867900e+07,8.497200e+05,-8.581400e+03,3.896300e+01,-1.395200e+07},
           						    {-1.343500e+09,-4.098400e+07,4.169100e+05,-1.873600e+03,6.606300e+08},
            						    {2.166000e+02,6.183200e+00,-5.853200e-02,2.295400e-04,-1.041800e+02}};

  const double kappa::Approximation::p4e_n2n_vt_cij5[5][5]={{1.118000e+02,3.463900e+00,-3.530000e-02,1.535600e-04,-5.523700e+01},
           						    {-7.693000e+03,-2.617800e+02, 2.849600e+00,-1.316300e-02,3.912200e+03},
           						    {-9.314000e+05,-2.472600e+04,2.224700e+02,-8.793600e-01,4.399400e+05},
           						    {-2.386600e+07,1.438500e+06,-1.402100e+04, 5.974700e+01,-2.386600e+07},
           						    {2.166000e+02,-2.386600e+07,4.526400e-03,-2.007900e-05,6.831900e+00}};

  const double kappa::Approximation::p4e_n2n_vt_bij6[5][5]={{7.730100e+01,7.730100e+01,7.730100e+01},
            						    {1.449000e+03,1.449000e+03,4.050100e-01},
            						    {2.492500e+08,-9.691600e+06,9.060600e+04},
            						    {-6.695300e+02,2.100000e+01,-1.808600e-01}};

  const double kappa::Approximation::p4e_n2n_vt_cij6[5][5]={{-2.929200e+00,9.561900e-02,-8.020000e-04},
	            					    {-1.393200e+02,4.419500e+00,-3.464200e-02},
            						    {-5.928900e+06,2.393900e+05,-2.354800e+03},
            						    {1.660300e+01,-5.307500e-01,4.228100e-03}};

  const double kappa::Approximation::p4e_o2o_vt_Cijk1[3][5][4]={{{-26.18993227,0,0,0},
                       						 {-1.69828917,0,0,0},
                       					         {3.349076e19,0,0,0},
                       						 {-3.946126e20,0,0,0},
                       						 {1.391056e19,0,0,0}},
                      						 {{7.833311e+00,0,0,0},
                       						 {3.712214e+00,0,0,0},
                       						 {3.573261e+20,0,0,0},
                       						 {6.433503e+20,0,0,0},
                       						 {-2.901352e+19,0,0,0}},
                      						 {{3.716395e-01,0,0,0},
                       						 {1.058709e-01,0,0,0},
                       						 {-5.312491e+19,0,0,0},
                       						 {3.754092e+19,0,0,0},
                       						 {-1.189832e+18,0,0,0}}};

  const double kappa::Approximation::p4e_o2o_vt_Cijk2[3][5][4]={{{-1.945676e+01,-3.380076e+00,8.985159e+01,5.853646e-02},
                                                                 {-1.487120e+00,-5.401834e-01,-5.333457e+01,-9.543789e-02},
                                                                 {1.505136e+21,1.621622e+21,-1.066183e+21,2.169160e+20},
                                                                 {-1.532916e+20,-4.105380e+19,1.185042e+21,-1.748646e+18},
                                                                 {-4.838473e+18,9.529399e+18,5.290114e+19,6.021495e+16}},
                                                                 {{1.447822e+01,-4.332225e+01,-3.481680e+02,-1.641860e-01},
                                                                 {3.266821e+01,-1.846219e-01,1.190313e+02,2.225885e-01},
                                                                 {1.522507e+21,2.654567e+21,-3.528669e+21,-2.861293e+20},
                                                                 {-1.533872e+21,-1.522587e+20,-3.124046e+21,2.322209e+19},
                                                                 {-3.762650e+19,1.955942e+19,8.847719e+18,-8.252347e+17}},
                                                                 {{8.865765e-01,-1.142588e+00,-4.530804e+00,4.732032e-03},
                                                                 {-1.856692e-01,-1.925642e-01,3.958191e+00,1.353007e-02},
                                                                 {-2.027215e+21,2.381051e+21,-1.248596e+20,-4.014395e+19},
                                                                 { 1.819921e+19,3.708631e+19,2.031805e+19,-6.021031e+17},
                                                                 { 1.530423e+18,-1.859900e+18,-8.207403e+18,9.708459e+15}}};

  const double kappa::Approximation::p4e_o2o_vt_Cijk3[3][5][4]={{{-3.993663e+00,2.684515e+00,-1.009927e+05,-2.836760e-01},
                                                               {-3.030965e+00,-4.594443e+00,3.590903e+04, 7.104764e-02},
                                                               {5.492061e+21,1.212196e+21,9.092488e+21,1.540038e+20},
                                                               {1.308503e+20,1.831856e+20,1.079540e+22,-5.608629e+18},
                                                               {2.160753e+19,-1.465914e+19,5.483520e+21,1.142128e+17}},
                                                               {{-7.575157e+01,-9.234026e+00,2.807916e+05,3.333970e-01},
                                                               {-7.713850e+00,2.545634e+01,-9.592245e+04,-1.792620e-01},
                                                               {4.002520e+21,-8.192010e+21,1.462011e+23,-1.821224e+20},
                                                               {-2.912948e+21,6.399791e+19,-3.531505e+22, 1.964744e+19},
                                                               {-7.070723e+19,4.805948e+19,5.201014e+21,-4.170396e+17}},
                                                               {{-1.271181e+00,3.178340e-01,5.830886e+03,7.186328e-03},
                                                               {-7.090280e-01, 4.753706e-02,-1.757570e+03,1.465161e-03},
                                                               {-8.686388e+20,9.428176e+20, 1.738719e+22,-1.251868e+19},
                                                               {1.877581e+19,1.097908e+19,-1.006633e+22,2.734279e+16},
                                                               {2.238756e+18,-7.688157e+17,-2.061752e+21,-3.427940e+15}}};

  const double kappa::Approximation::p4e_o2o_vt_Cijk4[3][5][4]={{{2.821759e+00,4.083138e+00,-8.809991e+01,2.369644e-02},
                                                                 {3.841080e+00,-9.986370e-01,-3.438479e+01,1.222381e-02},
                                                                 { 4.330479e+21,-1.677646e+22,-5.573334e+21,-3.089812e+19},
                                                                 {-1.194045e+20,-5.121704e+19,4.013656e+21,-1.052730e+18},
                                                                 {-9.939380e+18,4.180038e+18,-2.265448e+18,1.050411e+16}},
                                                                 {{-1.066105e+02,6.618737e+01,2.630406e+02,2.791153e-02},
                                                                 {-2.658559e+01,-9.167211e+00,1.678357e+02,-1.064740e-01},
                                                                 {1.312884e+22,-5.437653e+21, 5.735816e+21,1.568233e+20},
                                                                 {4.530951e+21,2.662341e+20,-2.932068e+22,6.788371e+18},
                                                                 {-3.473472e+19,-5.623449e+18,2.765213e+20,-6.030509e+16}},
                                                                 {{-3.476825e+00,-4.156750e-01,2.341590e+01,-1.760866e-03},
                                                                 {-5.624110e-01,1.663190e-01,2.356659e+00,-7.409818e-04},
                                                                 {1.908092e+21,1.107010e+21,-1.769244e+22,1.272578e+18},
                                                                 {-2.139241e+19,2.171483e+19,-2.478535e+19,7.612186e+16},
                                                                 {1.542813e+18,-7.531694e+17,-4.924709e+18,-5.546602e+14}}};

  const double kappa::Approximation::p4e_n2n_diss_bjk[4][5]={	{-4.102280e+01,-1.132030e+05,-5.081880e+05,6.584820e+07,2.009300e+00},
        							{-3.952490e-01,3.502120e+03,1.616070e+04,-3.360490e+06,6.080030e-02},
        							{2.810930e-02,-3.104770e+01,2.812760e+01,5.289150e+04,-3.362940e-03},
        							{-2.727510e-04,6.024030e-02,-2.413860e+00,-2.574970e+02,3.166530e-05}};

  const double kappa::Approximation::p4e_o2o_diss_bjk[3][8]={	{-27.0145,-3.19019,1.05845,-0.126656,0.00721311,-0.000211765,3.10206e-06,-1.79668e-08},
 			             				{-59620,2178.93,-11.2525,1.59843,-0.177455,0.00687121,-0.000119773,7.84053e-07},
              							{0.675455,0.376682,-0.121598,0.0144864,-0.00082352,2.41602e-05,-3.53967e-07,2.05237e-09}};

  const double kappa::Approximation::cc_ar[5] = {0., 1., 1.5, 1.83333333, 2.08333333333333};
  const double kappa::Approximation::cc_cl[5] = {0.5, 1., 1.1666666666666667, 1.333333333333333, 1.43333333333333};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
