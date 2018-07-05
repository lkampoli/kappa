/* 
 * Test for caching in mixture transport coefficients computation
 */

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

    std::cout << "Start test: computation of shear viscosity" << std::endl;
 
    std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
    std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << std::endl;
    std::string particle_source    = m_source + "particles.yaml";
    std::string interaction_source = m_source + "interaction.yaml";
    std::string output_dir = GetCurrentWorkingDir();
    std::cout << "Current directory is: " << output_dir << std::endl;
   
    std::vector<kappa::Molecule> molecules;
    std::vector<kappa::Atom> atoms;

    std::cout << "Loading particles data" << std::endl;

    kappa::Molecule mol("N2", true, false, particle_source);
    kappa::Atom at(     "N",              particle_source);

    molecules.push_back(mol);
    atoms.push_back(at);

    std::cout << "Finished loading particles data" << std::endl;

    // create N2-N binary mixture
    kappa::Mixture mixture(molecules, atoms, interaction_source, particle_source);
    // create a pure N2 mixture
    kappa::Mixture mixture_pure(molecules,   interaction_source, particle_source);

    std::vector<double> T_vals = { 300., 500., 1000., 2000., 10000., 40000. };
    double x_atom = 0.5;
   
    std::vector<arma::vec> mol_ndens;
    mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (kappa::K_CONST_K * T_vals[0]), mol));
    arma::vec atom_ndens(1);

    std::ofstream outf;
    std::ofstream outf2;
    //outf.open(output_dir + "/nocache_" + mol.name + "_" + at.name + "_xat50.txt");
    outf2.open(output_dir + "/cache_" + mol.name + "_" + at.name + "_xat50.txt");
 
    outf2 << std::setw(20) << "Temperature [K]";
    outf2 << std::setw(20) << "Lambda";
    outf2 << std::setw(20) << "Eta";
    outf2 << std::endl;

    double tot_ndens;
    double sh_v, th_c;

    std::cout << std::setw(20) << "Temperature [K]" << std::endl;
    for (auto T : T_vals) {
       tot_ndens =  101325.0 / (kappa::K_CONST_K * T);
       mol_ndens[0] = mixture.Boltzmann_distribution(T, (1 - x_atom) * tot_ndens, mol);
       atom_ndens[0] =  x_atom * tot_ndens;

       mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens);
//     sh_v = mixture.get_thermal_conductivity();
//     th_c = mixture.get_shear_viscosity();

//       outf2 << std::setw(20) << T;
//       outf2 << std::setw(20) << sh_v;
//       outf2 << std::setw(20) << th_c;
//       outf2 << std::endl;
    }

    //outf.close();
    outf2.close();

    return 0;
}
