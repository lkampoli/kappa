/*! 
    \file: mixture-sts-shear-keras.cpp 
    \brief Test for shear viscosity in binary mixtures
           computed with keras model regression
           imported by frugally-deep interface.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip> 

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#include "kappa.hpp"

#include <fdeep/fdeep.hpp>

using namespace kappa;

std::string GetCurrentWorkingDir( void ) {
  char buff[FILENAME_MAX];
  GetCurrentDir( buff, FILENAME_MAX );
  std::string current_working_dir(buff);
  return current_working_dir;
}

int main(int argc, char** argv) {

  std::cout << "Start test: computation of shear viscosity" << std::endl;

  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
  std::string particle_source    = m_source + "particles.yaml";
  std::string interaction_source = m_source + "interaction.yaml";
  std::string output_dir = GetCurrentWorkingDir();

  std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
  std::cout << "Current directory is: " << output_dir << std::endl;
  std::cout << "Particle dir: " << particle_source << '\n';
  std::cout << "Interaction dir: " << interaction_source << '\n';

  std::vector<kappa::Molecule> molecules;
  std::vector<kappa::Atom> atoms;

  std::cout << "Loading particles data" << std::endl;

  kappa::Molecule mol("N2", true, false, particle_source);
  kappa::Atom at("N", particle_source);

  molecules.push_back(mol);
  atoms.push_back(at);

  kappa::Mixture mixture(molecules, atoms, interaction_source, particle_source);
  std::cout << "Mixture created" << std::endl;

  //std::vector<double> T_vals = { 500., 1000., 5000., 10000., 20000., 40000. };
  std::vector<float> T_vals = { 5000. };
  float x_N2 = 0.9;
  float x_O2 = 0.0;
  float x_NO = 0.0;
  float x_N  = 0.1;
  float x_O  = 0.0;
  float pressure = 101325.;

  std::cout << T_vals[0] << " " << pressure << " " << x_N2 << " " << x_O2 << " " << x_NO << " " << x_N << " " << x_O << std::endl;

  // for (int i=0; i<80; i++) { //500-40000K
  //   T_vals.push_back(500 + i * 500);
  // }

  /////////////////   
  const auto model = fdeep::load_model("shear_fdeep_model.json");
  const auto result = model.predict(
          {fdeep::tensor(fdeep::tensor_shape(static_cast<std::size_t>(7)),
                  std::vector<float>{T_vals[0], pressure, x_N2, x_O2, x_NO, x_N, x_O})});
  std::cout << std::setw(20) << " shear viscosity KERAS = " << fdeep::show_tensors(result) << std::endl;
  /////////////////   
  
  std::vector<arma::vec> mol_ndens;
  mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], x_N2 * pressure / (K_CONST_K * T_vals[0]), mol)); // N2

  std::ofstream outf;
  outf.open(output_dir + "/shear_viscosity_" + mol.name + "_xat_" + std::to_string(x_N2) + ".txt");
  outf << std::setw(20) << "Temperature [K]";
  outf << std::setw(20) << "Eta";
  outf << std::endl;

  arma::vec atom_ndens(2);
  atom_ndens[0] = 0.;
  atom_ndens[1] = 0.;

  double tot_ndens;

//  std::vector<models_omega> omega_integral_models = {models_omega::model_omega_rs, models_omega::model_omega_vss,
//                                                     models_omega::model_omega_bornmayer, models_omega::model_omega_lennardjones,
//                                                     models_omega::model_omega_esa};

  for (auto T : T_vals) {
    tot_ndens =  pressure / (K_CONST_K * T);
    mol_ndens[0] = mixture.Boltzmann_distribution(T, x_N2 * pressure / (K_CONST_K * T), mol);
    atom_ndens[0] = x_N * pressure / (K_CONST_K * T);

    mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens, 0, models_omega::model_omega_rs, 0.);

    outf << std::setw(20) << T;
    outf << std::setw(20) << mixture.get_shear_viscosity();
    outf << std::endl;
  }

  std::cout << std::setw(20) << " shear viscosity = " << mixture.get_shear_viscosity() << std::endl;
  //std::cout << std::setw(20) << " bulk viscosity = " << mixture.get_bulk_viscosity() << std::endl;
  //std::cout << std::setw(20) << " thermal conductivity = " << mixture.get_thermal_conductivity() << std::endl;

  //arma::vec thd; 
  //thd = mixture.get_thermodiffusion();
  //std::cout << " thermal diffusion = " << mixture.get_thermodiffusion() << std::endl;

  //for (auto moles: molecules) {
  //  std::cout << std::setw(20) << mol_ndens[0];
  //  std::cout << std::endl;
  //}

  outf.close();
  return 0;
}
