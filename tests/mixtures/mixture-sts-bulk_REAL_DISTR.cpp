#include <iostream>
#include "kappa.hpp"
#include <fstream>
#include <string>
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
    std::cout << "Start test: computation of bulk viscosity" << std::endl;

    std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
    std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
    std::string particle_source    = m_source + "particles.yaml";
    std::string interaction_source = m_source + "interaction.yaml";
    std::cout << "particle_source: " << particle_source << std::endl;
    std::cout << "interaction_source: " << interaction_source << std::endl;
    std::string output_dir = GetCurrentWorkingDir();
    std::cout << "Current directory is: " << output_dir << std::endl;    

    std::vector<kappa::Molecule> molecules;
    std::vector<kappa::Atom> atoms;

    std::cout << "Loading particles data" << std::endl;

    kappa::Molecule mol("N2", true, false, particle_source);
    kappa::Atom at("N", 		   particle_source);
    //kappa::Molecule mol("O2", true, false, particle_source);
    //kappa::Atom at("O", 		   particle_source);

    std::cout << "Finished loading particles data" << std::endl;

    kappa::Mixture mixture(mol, at, interaction_source, particle_source);
    std::cout << "Mixture created!" << std::endl;

    std::vector<double> T_vals;
    int i;
    for (i = 0; i < 80; i++) {
     T_vals.push_back(500 + i * 500); // 500-40000 K
    }

    double x_atom_perc = 0.; // 0 - 20 - 50
    double x_atom = x_atom_perc / 100.0;

    std::vector<arma::vec> mol_ndens;
    mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (kappa::K_CONST_K * T_vals[0]), mol));
    arma::vec atom_ndens(1);

    std::ofstream outf;
    outf.open(output_dir + "/TRANSPORT_COEFFICIENTS/bulk_viscosity/" +  mol.name + "_" + at.name + "_xat" + std::to_string(x_atom_perc) + ".txt");

    std::cout << "Output file opened!" << std::endl;
    outf << std::setw(20) << "Temperature [K]";
    outf << std::setw(20) << "zeta";
    //outf << std::setw(20) << "zeta / eta";
    outf << std::endl;

    double tot_ndens;

    for (auto T : T_vals) {

     tot_ndens =  101325.0 / (kappa::K_CONST_K * T); // atmospheric pressure
     mol_ndens[0] = mixture.Boltzmann_distribution(T, (1 - x_atom) * tot_ndens, mol);
     atom_ndens[0] = x_atom * tot_ndens;

     mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens);

     outf << std::setw(20) << T;
     outf << std::setw(20) << mixture.get_bulk_viscosity();
     //outf << std::setw(20) << mixture.get_bulk_viscosity() / mixture.get_shear_viscosity();
     outf << std::setw(20) << mixture.get_bulk_viscosity() / tot_ndens;
     outf << std::endl;
      
    }

    outf.close();
    return 0;
}
