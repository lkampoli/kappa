
#include <iostream>
#include "kappa.hpp"
//#include "utilities.hpp"
#include <fstream>
#include <string>
#include <iomanip> // std::setw

using namespace kappa;

int main(int argc, char** argv) {

    std::cout << "Start Test state-to-state mixture" << std::endl;

//  kappa::ReadDataMixture();

    std::string m_libPath = std::getenv("KAPPA_DATA_DIRECTORY");
//  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
//  std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
    std::string particle_source    = m_libPath + "particles.yaml";
    std::string interaction_source = m_libPath + "interaction.yaml";
//  std::string output_dir = GetCurrentWorkingDir();
//  std::cout << "Current directory is: " << output_dir << std::endl;


    // Instantiates molecules and atoms kappa vectors
    std::vector<kappa::Molecule> molecules;
    std::vector<kappa::Atom> atoms;

    std::cout << "Loading particles data" << std::endl;

    // Chooses some particles from file
    kappa::Molecule N2("N2", true, true, particle_source);
    kappa::Molecule O2("O2", true, true, particle_source);
    kappa::Molecule NO("NO", true, true, particle_source);
    kappa::Atom N("N",                   particle_source);
    kappa::Atom O("O",                   particle_source);

    // Populate the vectors with such particles
    molecules.push_back(N2);
    molecules.push_back(O2);
    molecules.push_back(NO);
    atoms.push_back(N);
    atoms.push_back(O);
   
    std::cout << "Finished loading particles data" << std::endl;

    std::cout << "Size of molecules: " << molecules.size() << std::endl;
    std::cout << "Size of atoms: " << atoms.size() << std::endl;

    // Instantiates the mixture with those molecules, atoms
    // and defines an interaction from file
    kappa::Mixture mixture(molecules, atoms, interaction_source, particle_source);

    for (auto at: atoms) {
        for (auto mo: molecules) {
            std::cout << at.name << "+" << mo.name << ", interaction=" 
                      << mixture.interaction(at, mo).particle1_name << "+" 
   		      << mixture.interaction(at, mo).particle2_name << std::endl;
        }
    }

    for (auto at: atoms) {
        for (auto at2: atoms) {
            std::cout << at.name << "+" << at2.name << ", interaction=" 
		      << mixture.interaction(at, at2).particle1_name << " " 
		      << mixture.interaction(at, at2).particle2_name << std::endl;
        }
    }

    for (auto mo: molecules) {
        for (auto mo2: molecules) {
            std::cout << mo.name << "+" << mo2.name << ", interaction=" 
		      << mixture.interaction(mo, mo2).particle1_name << " " 
		      << mixture.interaction(mo, mo2).particle2_name << std::endl;
        }
    }

     std::vector<double> T_vals = { 500., 1000., 5000., 10000., 20000., 40000. };
     std::vector<arma::vec> mol_ndens;
     mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (K_CONST_K * T_vals[0]), N2));
     mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (K_CONST_K * T_vals[0]), O2));
     mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (K_CONST_K * T_vals[0]), NO));

     std::cout << std::setw(20) << " Temp. [K] ";
     std::cout << std::setw(20) << " N2 ";
     std::cout << std::setw(20) << " O2 ";
     std::cout << std::setw(20) << " NO ";
     std::cout << std::endl;
     for (auto T : T_vals) {
        mol_ndens[0] = mixture.Boltzmann_distribution(T, 101325.0 / (K_CONST_K * T), N2);
        mol_ndens[1] = mixture.Boltzmann_distribution(T, 101325.0 / (K_CONST_K * T), O2);
        mol_ndens[2] = mixture.Boltzmann_distribution(T, 101325.0 / (K_CONST_K * T), NO);
        std::cout << std::setw(20) << T;
        std::cout << std::setw(20) << mol_ndens[0][0] / (101325.0 / (K_CONST_K * T));
        std::cout << std::setw(20) << mol_ndens[1][0] / (101325.0 / (K_CONST_K * T));
        std::cout << std::setw(20) << mol_ndens[2][0] / (101325.0 / (K_CONST_K * T));
        //std::cout << std::setw(20) << mol_ndens[0] / (101325.0 / (K_CONST_K * T));
        //std::cout << std::setw(20) << mol_ndens[1] / (101325.0 / (K_CONST_K * T));
        //std::cout << std::setw(20) << mol_ndens[2] / (101325.0 / (K_CONST_K * T));
        std::cout << std::endl;
     }
     arma::vec atom_ndens(2);

//     double x_N = 0.30; 
//     double x_O = 0.70;
 
     arma::vec x_molar_fractions = arma::zeros(molecules.size()+atoms.size());
     x_molar_fractions[0] = 0.1; // N2
     x_molar_fractions[1] = 0.2; // O2
     x_molar_fractions[2] = 0.3; // NO
     x_molar_fractions[3] = 0.2; // N
     x_molar_fractions[4] = 0.2; // O

//     atom_ndens[1] = 0.0;
//     atom_ndens[2] = 0.0;
//     atom_ndens[3] = x_N * 101325.0 / (K_CONST_K * 1000);
//     atom_ndens[4] = x_O * 101325.0 / (K_CONST_K * 1000);
//
//     std::cout << std::setw(20) << "N2 density: " <<  mixture.compute_density(mol_ndens[0]) << "\n";
//     std::cout << std::setw(20) << "O2 density: " <<  mixture.compute_density(mol_ndens[1]) << "\n";
//     std::cout << std::setw(20) << "NO density: " <<  mixture.compute_density(mol_ndens[2]) << "\n";
//
//     std::cout << std::setw(20) << "compute atoms' density: " <<  mixture.compute_density(atom_ndens) << "\n";
//     std::cout << std::setw(20) << "compute mols' density: " <<  mixture.compute_density(mol_ndens) << "\n";
//     std::cout << std::setw(20) << "compute density (at + mol): " <<  mixture.compute_density(mol_ndens,atom_ndens,0.0) << "\n";
//
//     std::cout << std::setw(20) << "compute n (at + mol): " <<  mixture.compute_n(mol_ndens, atom_ndens, 0.0) << "\n";
//     std::cout << std::setw(20) << "compute n (mol): " <<  mixture.compute_n(mol_ndens, 0.0) << "\n";
//     std::cout << std::setw(20) << "compute n mol: " <<  mixture.compute_n_molecule(mol_ndens) << "\n";
//     std::cout << std::setw(20) << "compute n (at): " <<  mixture.compute_n(atom_ndens) << "\n";
     
     arma::vec x_mass_fractions = arma::zeros(molecules.size()+atoms.size());
     x_mass_fractions = mixture.convert_molar_frac_to_mass(x_molar_fractions);
     //std::cout << std::setw(20) << " x_mass_fractions: " <<  x_mass_fractions << "\n";
     std::cout << std::setw(20) << x_mass_fractions << "\n";
     std::cout << std::setw(20) << x_mass_fractions.size() << "\n";
     std::cout << std::setw(20) << x_mass_fractions.n_rows << "\n";
     std::cout << std::setw(20) << x_mass_fractions.n_cols << "\n";

     std::vector<double> mass_fractions;
     //mass_fractions(x_mass_fractions.n_rows);
     for ( int i = 0; i < x_mass_fractions.n_rows; ++i) {
        mass_fractions.push_back(x_mass_fractions(i));
        std::cout << std::setw(20) << mass_fractions[i] << "\n";
    };

    return 0;
}
