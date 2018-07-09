/*! 
    \file mixture-sts-multi-component-diffusion_air5.cpp 
    \brief Test for multi-component diffusion coefficients' computation for air5 mixture in the StS approach.
 */

#include <iostream>
#include <fstream>
#include <iomanip> 

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#include "kappa.hpp"

std::string GetCurrentWorkingDir( void ) {
  char buff[FILENAME_MAX];
  GetCurrentDir( buff, FILENAME_MAX );
  std::string current_working_dir(buff);
  return current_working_dir;
}

int main(int argc, char** argv) {
  
  std::cout << "Start test: computation of diffusion coefficients" << std::endl;

  std::string line;
  int linecount=0;
  int _NS=0;
  int _nMol=0;
  int _nAt=0;

  std::vector<std::string> _speciesNames; 
  std::vector<std::pair<std::string,double>> composition;

  // mixture selection
  std::string _mixtureName = "air5";   
  std::cout << "Mixture: " << _mixtureName << "\n";   

  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
  std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
  std::string particle_source    = m_source + "particles.yaml";
  std::string interaction_source = m_source + "interaction.yaml";
  std::string output_dir = GetCurrentWorkingDir();
  std::cout << "Current directory is: " << output_dir << std::endl;

  std::ifstream mixfile((m_source+"mixtures/"+_mixtureName+".mix").c_str());
   
  // open and read input file
  if (mixfile.is_open()) {
    while ( mixfile.good() ) {
      getline(mixfile,line);

      if (line[0] != '/') {
        linecount+=1;
        if (linecount==1) {
          _NS= static_cast<int> (atoi(line.c_str()));
          _speciesNames.resize(_NS);
          // checks the species vector size
          // std::cout << "check size of _speciesNames: " << _speciesNames.size() << std::endl;
        }
        if (linecount>1 && linecount<_NS+2)  _speciesNames[linecount-2]=line;
      }
    }
    mixfile.close();
  } else {
    std::cout << "Unable to open mixture file: "
              << m_libPath+"mixtures/"+_mixtureName+".mix" << "\n";
              abort();
  }
  std::cout << "Number of species: " << _NS << "\n";
  std::cout << "Loading particles data ..." << std::endl;

  // instantiates molecules and atoms kappa vectors
  std::vector<kappa::Molecule> molecules;
  std::vector<kappa::Atom> atoms;

  // instantiates temperature vector
  std::vector<double> T_vals = { 500., 1000., 5000., 10000., 20000., 40000. };
  int j = 1;

/*  
    We will assume that element names are no more than 2 characters
    long and that the species name is in "standard" form which is AaBbCc+
    where the upper case letters represent element names, lower case letters
    are the optional stoichiometric numbers for the element, and the '+'
    represents either a '+' or '-' character which represents minus or plus
    one electron respectively. If the name does not fit this form, then we
    forego the check.
*/

  // checks the species vector and create molecules and atoms
  for (int i=0; i<_NS; ++i) {
    // std::cout << _speciesNames[i] << " " <<  _speciesNames[i].length() << std::endl;

    // Step 1: Check to see if the name is in standard form
    bool standard = isalpha(_speciesNames[i].at(0));
    // std::cout << "_speciesNames[i].length()-1: " << _speciesNames[i].length()-1 << std::endl; 
    if (_speciesNames[i].length() > 1) {
      for (j=0; j<_speciesNames[i].length()-1; ++j)
        standard = standard && isalnum(_speciesNames[i].at(j));
        standard = standard && (isalnum(_speciesNames[i].at(j)) ||
                                        _speciesNames[i].at(j) == '+' ||
                                        _speciesNames[i].at(j) == '-');
      } // /=1
      if (!standard) {
        std::cout << standard << " ouch! " << std::endl;
        return 0;
        // std::cout << "_speciesNames[i].at(j) : " << _speciesNames[i].at(j) << std::endl; }
      }

      int length = _speciesNames[i].length();
      std::cout << length << std::endl;

      // Step 2: Now we know that the name is in standard form, parse it to get the stoichiometry
      // Handle the special case of the electron
      if (_speciesNames[i] == "e-") {
        length = 0;
      } else if (_speciesNames[i].at(length-1) == '+') {
        length--;
      } else if (_speciesNames[i].at(length-1) == '-') {
        length--;
      }
      // std::cout << _speciesNames[i] << " " <<  length << std::endl;

      // At this point we know that the name must be AaBbCc, use finite-state machine to parse
      if (length==0 || length == 1) { // electrons or atoms
        kappa::Atom at(_speciesNames[i], particle_source);
        atoms.push_back(at);
      } else if ( length > 1 ) {
        kappa::Molecule mol(_speciesNames[i], true, false, particle_source);
        molecules.push_back(mol);
      }
  }

  std::cout << " Loading particles data completed!" << std::endl;

  _nMol= molecules.size();
  _nAt = atoms.size();
  std::cout << " Size of molecules: " << _nMol << std::endl;
  std::cout << " Size of atoms: " << _nAt << std::endl;

  // instantiates the mixture object and some print info
  kappa::Mixture mixture(molecules, atoms, interaction_source, particle_source);
  std::cout << " The mixture has been created! " << std::endl;
  std::string MixtureComponents = mixture.get_names();
  std::cout << " Mixture's components: " << MixtureComponents << std::endl;
  std::cout << " Number of component species: " << mixture.get_n_particles() << std::endl;

  std::cout << "Finished loading particles data" << std::endl;


  int vibr_l = 0.; // for N2 d0
  double Re = 1.097E-10; // N2
  double be = 2.25E-10;
  double omega_e = 235860.; // m^-1 (N2)
  double mu = 0.028; // kg/mol (N2)
  double l_alpha = sqrt(16.863/(omega_e * mu));
  double beta = 2.6986E+10; // N2
  double d0 = Re + be + (9./2.)*beta*l_alpha*l_alpha*exp( 2.*sqrt(beta*l_alpha)*(vibr_l - 1.) );  

  // set a range for temperature
  std::vector<double> T_vals;

  double p = 101325.;

  // vibrational levels
  int i;
   for (i = 0; i < 80; i++) {
    T_vals.push_back(500 + i * 500); // 500-40000 K
   }

  // set an arbitrary atom mass fraction (0.20)
  int x_atom_perc = 25.;
  double x_atom = x_atom_perc / 100.;

  // arma vector for atom number density
  arma::vec atom_ndens(2);

  // vector of arma vector for molecular number density
  std::vector<arma::vec> mol_ndens;

  // mixture number density
  double tot_ndens;
  tot_ndens =  p / (kappa::K_CONST_K * T_vals[0]);
  std::cout << " tot_ndens " << tot_ndens << std::endl;

  for (auto mol : molecules) {
    mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], (1. - 2.*x_atom) * tot_ndens, mol, 0));
  }

  for (auto i=mol_ndens.begin(); i!=mol_ndens.end(); ++i) {
    std::cout << (*i) << std::endl;
  }

  // atom number density
  // atom_ndens = x_atom * tot_ndens;
  for (int at=0; at<atom_ndens.size(); at++) {
    atom_ndens[at] = x_atom * tot_ndens;
    std::cout << " atom_ndens " << atom_ndens.at(at) << std::endl;
  }

  std::ofstream outf;
  outf.open(output_dir + "/TRANSPORT_COEFFICIENTS/bulk_viscosity/" + "air5_xat" + std::to_string(x_atom_perc) + ".txt");
  std::cout << "Output file opened!" << std::endl;
  outf << std::setw(20) << "Temperature [K]";
  outf << std::setw(20) << "zeta";
  outf << std::setw(20) << "zeta / n_tot";
  outf << std::endl;

  // main loop on temperatures
  for (auto T : T_vals) {

    tot_ndens =  101325.0 / (kappa::K_CONST_K * T);

    for (auto mol : molecules) {
      mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], (1. - 2.*x_atom) * tot_ndens, mol, 0));
    }

    // for (auto i=mol_ndens.begin(); i!=mol_ndens.end(); ++i) {
    //   std::cout << (*i) << std::endl;
    // }

    // for (int at=0; at<atom_ndens.size(); at++) {
    //   atom_ndens[at] = x_atom * tot_ndens;
    //   std::cout << " atom_ndens " << atom_ndens.at(at) << std::endl;
    // }

    // compute transport coefficients 
    mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens);
   
    // retrieve bulk viscosity
    outf << std::setw(20) << T;
    outf << std::setw(20) << mixture.get_bulk_viscosity();
    outf << std::setw(20) << mixture.get_bulk_viscosity() / tot_ndens;
    outf << std::endl;  
   
  }
  outf.close();
  return 0;
}
