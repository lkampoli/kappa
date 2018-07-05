
// Computation and comparison of all (available) transport coefficients

#include <iostream>
#include "kappa.hpp"
#include <fstream>
#include <iomanip> // std::setw

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

  std::string GetCurrentWorkingDir( void ) {
    char buff[FILENAME_MAX];
    GetCurrentDir( buff, FILENAME_MAX );
    std::string current_working_dir(buff);
    return current_working_dir;
  }

   int main(int argc, char** argv) {
   
   std::cout << "Start computation of transport coefficients" << std::endl;
   
   std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
   std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
   std::string particle_source    = m_source + "particles.yaml";
   std::string interaction_source = m_source + "interaction.yaml";
   std::string output_dir = GetCurrentWorkingDir();
   std::cout << "Current directory is: " << output_dir << std::endl;

   std::cout << "Loading particles data" << std::endl;

   // N2 molecule
   //kappa::Molecule mol("N2", true, true, particle_source); // rigid rotator
   kappa::Molecule mol("N2", true, false, particle_source); // non-rigid rotator
  
   // N atom
   kappa::Atom at("N", particle_source);

   std::cout << "Finished loading particles data" << std::endl;

   std::cout << "Molecule's name " << mol.name << std::endl;
   std::cout << "Atom's name " << at.name << std::endl;
   std::cout << "Molecule vibrational levels " << mol.num_vibr_levels[0] << std::endl;

   // N2/N binary mixture creation
   kappa::Mixture mixture(mol, at, interaction_source, particle_source);

   // some check print
   std::cout << "particles: " << mixture.get_n_particles() << std::endl;
   std::cout << "names: " << mixture.get_names() << std::endl;

   // spectroscopic constant for diatomic molecules in ground electronic state
   int vibr_l = 0; // for N2 d0
   double Re = 1.097E-10; // N2
   double be = 2.25E-10;
   double omega_e = 235860; // m^-1 (N2)
   double mu = 0.028; // kg/mol (N2)
   double l_alpha = sqrt( 16.863/(omega_e * mu) );
   double beta = 2.6986E+10; // N2
   double d0 = Re + be + (9./2.)*beta*l_alpha*l_alpha*exp( 2*sqrt(beta*l_alpha)*(vibr_l - 1) );
    
   // set a range for temperature
   std::vector<double> T_vals;
   std::vector<double> X_vals;
   double tot_ndens;

   // vibrational levels
   int i;

   for (i = 0; i < 80; i++) { //500-40000K
      T_vals.push_back(500 + i * 500);
   }

   for (i = 0; i <= 20; i++) { //0.-1.
      X_vals.push_back(i/20.);
   }

   // arma vector for atom number density
   arma::vec atom_ndens(1);

   // vector of arma vector for molecular number density
   std::vector<arma::vec> mol_ndens;

   // assume an initial boltzmann ditribution at the lowest temperature FIXME why it is necessary?
   mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (kappa::K_CONST_K * T_vals[0]), mol));

   // prepare the output files
   std::ofstream outf, outf2, outf3, outf4;
   
   //outf.open(output_dir + "/TCs_with_T" + mol.name + "_" + at.name + "_Xat" + std::to_string(x_atom_perc) + ".txt");
   outf.open(output_dir + "/TCs_with_T_" + mol.name + "_" + at.name + ".txt");

   outf << std::setw(20) << "X"; 
   outf << std::setw(20) << "Temperature [K]"; 
   outf << std::setw(20) << "Thermal conduction."; 
   outf << std::setw(20) << "Shear viscosity.";
   outf << std::setw(20) << "Bulk viscosity.";
   outf << std::endl;

   // transport coefficients
   double th_c; 	// thermal conductivity
   double sh_v; 	// shear viscosity
   double bk_v; 	// bulk viscosity
   double bk_o_sh;      // ratio bulk/shear viscosity
   arma::vec th_d;	// thermal diffusion
   arma::mat diff;	// diffusion coeffs.
///////////////////////////////////////////////////////////////////////////////
// Now we compute the transport coefficients varying the temperature T

   // main loop on temperatures
   for (auto X : X_vals) {
    for (auto T : T_vals) {

       // total number density
       tot_ndens =  101325.0 / (kappa::K_CONST_K * T);

       // vibrational level population distribution
       // mol_ndens[0] = mixture.Boltzmann_distribution(T, (1 - x_atom) * tot_ndens, mol);
       mol_ndens[0] = mixture.Boltzmann_distribution(T, (1 - X) * tot_ndens, mol);
       //mol_ndens.push_back(mixture.Boltzmann_distribution(T, (1 - X) * tot_ndens, mol));
       // atom_ndens[0] = x_atom * tot_ndens;
       atom_ndens[0] = X * tot_ndens;

       std::cout << std::setw(20) << T << std::endl;

       // computation of transport coefficients
       mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens);

       // retrieve thermal conductivity coefficients
       th_c = mixture.get_thermal_conductivity();

       // retrieve shear viscosity coefficients
       sh_v = mixture.get_shear_viscosity();

       // retrieve bulk viscosity coefficients
       bk_v = mixture.get_bulk_viscosity();

       outf << std::setw(20) << X; 
       outf << std::setw(20) << T; 
       outf << std::setw(20) << th_c; 
       outf << std::setw(20) << sh_v;
       outf << std::setw(20) << bk_v;
       outf << std::endl;
    }
   }

/////////////////////////////////////////////////////////////////////////////////////////
// Now we compute the same quantities but varying the molar mass fractions X
   outf4.open(output_dir + "/TCs_with_X_" + mol.name + "_" + at.name + ".txt");

   outf4 << std::setw(20) << "Temperature [K]"; 
   outf4 << std::setw(20) << "Mass fraction []"; 
   outf4 << std::setw(20) << "Thermal conduction."; 
   outf4 << std::setw(20) << "Shear viscosity.";
   outf4 << std::setw(20) << "Bulk viscosity.";
   outf4 << std::endl;

   // main loop on mass fractions
   for (auto T : T_vals) {
   std::cout << std::setw(20) << T; 
    for (auto X : X_vals) {
    std::cout << std::setw(20) << X; 

      // total number density
      tot_ndens =  101325.0 / (kappa::K_CONST_K * T);

      // vibrational level population distribution
      mol_ndens[0] = mixture.Boltzmann_distribution(T, (1 - X) * tot_ndens, mol); // FIXME but why?
      //mol_ndens.push_back(mixture.Boltzmann_distribution(T, (1 - X) * tot_ndens, mol));
      atom_ndens[0] = X * tot_ndens;

      // computation of transport coefficients
      mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens);

      // retrieve thermal conductivity coefficients
      th_c = mixture.get_thermal_conductivity();

      // retrieve shear viscosity coefficients
      sh_v = mixture.get_shear_viscosity();

      // retrieve bulk viscosity coefficients
      bk_v = mixture.get_bulk_viscosity();

      outf4 << std::setw(20) << T; 
      outf4 << std::setw(20) << X; 
      outf4 << std::setw(20) << th_c; 
      outf4 << std::setw(20) << sh_v;
      outf4 << std::setw(20) << bk_v;
      outf4 << std::endl;
    }
   }

/////////////////////////////////////////////////////////////////////////////////////////
   // compute thermo-diffusion coefficients

   //outf2.open(output_dir + "/thermo_diffusion_" + mol.name + "_" + at.name + "_Xat" + std::to_string(x_atom_perc) + ".txt");
   outf2.open(output_dir + "/thermo_diffusion_" + mol.name + "_" + at.name + ".txt");
  
   //double x_atom = 50.; // just N2 mix
   double D0;

   for (auto X : X_vals) {
    for (auto T : T_vals) {

       tot_ndens =  101325.0 / (kappa::K_CONST_K * T);
 
       mol_ndens[0] = mixture.Boltzmann_distribution(T, (1 - X) * tot_ndens, mol);
       atom_ndens[0] = X * tot_ndens;
 
       mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens);
 
       // retrieve thermal-diffusion coefficients
       th_d = mixture.get_thermodiffusion();
       //D0 = (3./(8.* mixture.compute_n(mol_ndens[0] * d0 * d0) )) * sqrt( (kappa::K_CONST_K * T)/(kappa::K_CONST_PI*mol.mass) );
       D0 = (3./(8.* tot_ndens * d0 * d0) ) * sqrt( (kappa::K_CONST_K * T)/(kappa::K_CONST_PI*mol.mass) );
       th_d /= D0;
 
       for (i=0; i<th_d.n_elem; i++) {

          outf2 << std::setw(20) << X; 
          outf2 << std::setw(20) << T; 
          outf2 << std::setw(20) << i; 
          outf2 << std::setw(25) << std::setprecision(18) << th_d[i]; 
          outf2 << std::endl;
       }
    }
   }
/////////////////////////////////////////////////////////////////////////////////////////  
   // compute diffusion coefficients (from simplified binary mixture formulas)

   outf3.open(output_dir + "/diffusion_" + mol.name + "_" + at.name + ".txt");

   double x_atom = 0.; // just N2 mix
   for (auto T : T_vals) {

      double tot_ndens;
      tot_ndens =  101325.0 / (kappa::K_CONST_K * T);

      mol_ndens[0] = mixture.Boltzmann_distribution(T, (1 - x_atom) * tot_ndens, mol);
      atom_ndens[0] = x_atom * tot_ndens;

      std::cout << std::setw(20) << T << std::endl;

      mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens);

      // retrieve thermal-diffusion coefficients
      diff = mixture.get_diffusion();
      diff /= D0;
  
//      outf3 << std::setw(20) << T; 
//      outf3 << std::setw(20) << i;
//      outf3 << std::setw(20) << diff;
//      outf3 << std::endl;

      for (i=0; i<diff.n_rows; i++) {
         outf3 << std::setw(20) << T; 
         outf3 << std::setw(20) << i; 
         outf3 << std::setw(20) << diff(i,0);
         outf3 << std::endl;
      }
   }
      
   outf.close();
   outf2.close();
   outf3.close();
   outf4.close();
   return 0;
}
