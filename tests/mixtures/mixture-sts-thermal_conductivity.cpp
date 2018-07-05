
// Computation of thermal conductivity (with Maria's data)

#include <iostream>
#include "kappa.hpp"
#include <fstream>
#include <iomanip> // std::setw
#include <omp.h> // openmp
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
   kappa::Molecule mol("N2", true, false, particle_source); // anharmonic, non-rigid rotator
  
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

   arma::vec atom_ndens(1);
   std::vector<arma::vec> mol_ndens;

   double t=1000.;
   double tmax=40500.; //25500.
   double p=101325.0;

   int i;
   int x_atom_perc = 20.;
   double x_atom = x_atom_perc / 100.;
   double th_c; //thermal conductivity

   std::ofstream outf;
   
   outf.open(output_dir + "/TRANSPORT_COEFFICIENTS/thermal_conductivity/" + mol.name + "_" + at.name + "_xat" + std::to_string(x_atom_perc) + ".txt");

//   outf << std::setw(20) << "X"; 
//   outf << std::setw(20) << "Temperature [K]"; 
//   outf << std::setw(20) << "Thermal conductivity"; 
//   outf << std::endl;
     
   while (t<tmax) {

   double tot_ndens =  p / (kappa::K_CONST_K * t);
   atom_ndens[0] = x_atom * tot_ndens;

   mol_ndens.push_back(mixture.Boltzmann_distribution(t, p / (kappa::K_CONST_K * t), mol));
   mol_ndens[0] = mixture.Boltzmann_distribution(t, (1 - x_atom) * tot_ndens, mol);

   //std::cout << std::setw(20) << t << std::endl;c

//#pragma omp parallel

   // computation of transport coefficients
   mixture.compute_transport_coefficients(t, mol_ndens, atom_ndens);

   // retrieve thermal conductivity coefficients
   th_c = mixture.get_thermal_conductivity();

// std::cout << std::setw(20) << th_c << std::endl;
   outf << std::setw(20) << t; 
   outf << std::setw(20) << th_c; 
   outf << std::endl;
 
   t += 500; 
   }

   outf.close();
   return 0;
}
