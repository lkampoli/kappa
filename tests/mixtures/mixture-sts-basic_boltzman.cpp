
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip> // std::setw
#include <cstdlib>
#include <vector>
#include "kappa.hpp"


int main(int argc, char** argv) {

    std::cout << "Start Test state-to-state mixture" << std::endl;

    std::string line;
    int linecount=0;

    int _NS=0;
    int _nMol = 0;
    int _nAt = 0;

    std::vector<std::string> _speciesNames; // just empty declaration
    std::vector<std::pair<std::string,double>> composition;


//  mixture selection
//  std::string _mixtureName = "air5";		//OK!
//  std::string _mixtureName = "air11"; 	//NO+ not in DB
    std::string _mixtureName = "N2";		//OK!
//  std::string _mixtureName = "N2ions";	//OK!
//  std::string _mixtureName = "O2";		//OK!
//  std::string _mixtureName = "argon";		// bug 
    std::cout << "Mixture: " << _mixtureName << "\n";

//  some path check
    std::string m_libPath = std::getenv("KAPPA_DATA_DIRECTORY");
    std::cout << "KAPPA_DATA_DIRECTORY is: " << m_libPath << '\n';
    std::string particle_source    = m_libPath + "particles.yaml";
    std::string interaction_source = m_libPath + "interaction.yaml";
    std::cout << "particle_source is: " << particle_source << '\n';
    std::cout << "interaction_source is: " << interaction_source << '\n';

    std::ifstream mixfile((m_libPath+"mixtures/"+_mixtureName+".mix").c_str());

//  open and read input file
     if (mixfile.is_open()) {

      while ( mixfile.good() )
      {
       getline(mixfile,line);

       if (line[0] != '/') // skip lines starting with /
        {
         linecount+=1;
         if (linecount==1)
         {
          _NS= static_cast<int> (atoi(line.c_str()));
          _speciesNames.resize(_NS);

//        checks the species vector size
//        std::cout << "check size of _speciesNames: " << _speciesNames.size() << std::endl;

         }
          if (linecount>1 && linecount<_NS+2)  _speciesNames[linecount-2]=line;
// 	read composition from file
// 	else composition.push_back(std::make_pair(? , ?));
        }
      }
      mixfile.close();
   }
    else {
      std::cout << "Unable to open mixture file: "
                << m_libPath+"mixtures/"+_mixtureName+".mix" << "\n";
                abort();
    }
    std::cout << "Number of species: " << _NS << "\n";
    std::cout << "Loading particles data ..." << std::endl;

//    instantiates molecules and atoms kappa vectors
      std::vector<kappa::Molecule> molecules;
      std::vector<kappa::Atom> atoms;

 
//    instantiates temperature vector
      std::vector<double> T_vals = { 500., 1000., 5000., 10000., 20000., 40000. };
      int j = 1;

/*  
    We will assume that element names are no more than 2 characters
    long and that the species name is in "standard" form which is AaBbCc+
    where the upper case letters represent element names, lower case letters
    are the optional stoichiometric numbers for the element, and the '+'
    represents either a '+' or '-' character which represents minus or plus
    one electron respectively.  If the name does not fit this form, then we
    forego the check.
*/

//    checks the species vector and create molecules and atoms
      for ( int i = 0; i < _NS; ++i ) {
//      std::cout << _speciesNames[i] << " " <<  _speciesNames[i].length() << std::endl;

// 	Step 1: Check to see if the name is in standard form
	bool standard = isalpha(_speciesNames[i].at(0));
//	std::cout << "_speciesNames[i].length()-1: " << _speciesNames[i].length()-1 << std::endl; 
  	if (_speciesNames[i].length() > 1) {
	 for ( j=0; j < _speciesNames[i].length()-1; ++j) 
           standard = standard && isalnum(_speciesNames[i].at(j));
	   standard = standard && (isalnum(_speciesNames[i].at(j)) || 
				 	   _speciesNames[i].at(j) == '+' || 
				 	   _speciesNames[i].at(j) == '-');
	} // /=1
  	if (!standard) {
	 std::cout << standard << " ouch! " << std::endl;
	 return 0;
//	 std::cout << "_speciesNames[i].at(j) : " << _speciesNames[i].at(j) << std::endl; }
	}

 	int length = _speciesNames[i].length();
        std::cout << length << std::endl;

	// Step 2: Now we know that the name is in standard form, parse it to get the stoichiometry
	// Handle the special case of the electron
    	if (_speciesNames[i] == "e-") {
         length = 0;
    	} 
	else if (_speciesNames[i].at(length-1) == '+') {
         length--;
    	} 
	else if (_speciesNames[i].at(length-1) == '-') {
         length--;
    	}
//        std::cout << _speciesNames[i] << " " <<  length << std::endl;

 	// At this point we know that the name must be AaBbCc, use finite-state machine to parse
        if ( length == 0 || length == 1 ) { // electrons or atoms
	 kappa::Atom at(_speciesNames[i], particle_source);
         atoms.push_back(at);
        } 
        else if ( length > 1 ) { 
         kappa::Molecule mol(_speciesNames[i], true, false, particle_source);
         molecules.push_back(mol);
        }
       }

       std::cout << " Loading particles data completed!" << std::endl;

    _nMol= molecules.size();
    _nAt = atoms.size();
    std::cout << " Size of molecules: " << _nMol << std::endl;
    std::cout << " Size of atoms: " << _nAt << std::endl;

    kappa::Mixture mixture(molecules, atoms, interaction_source, particle_source);
    std::cout << " The mixture has been created! " << std::endl;
    std::string MixtureComponents = mixture.get_names();
    std::cout << " Mixture's components: " << MixtureComponents << std::endl;
    std::cout << " Number of component species: " << mixture.get_n_particles() << std::endl;

    double molmass=0.; // just for 1 molecule mixture FIXME: generalize
    for (auto mol : molecules) {
      std::cout << "Molecule name: " << mol.name << std::endl;
      std::cout << "Molecule mass: " << mol.mass << std::endl;
      molmass = mol.mass;
    }
    std::cout << "Molecule mass: " << molmass << std::endl;

    double T0 = 1833.; //42.6; //233.61;
    //double T0 = 42.6; //233.61;
    double p0 = 2908.8; //2.23; //668.1;
    //double p0 = 2.23; //668.1;
    double tot_ndens =  p0 / (kappa::K_CONST_K * T0);

   std::vector<arma::vec> mol_ndens;
   arma::vec atom_ndens(1);

   double x_N=0.;
   atom_ndens[0] = x_N * (p0 / (kappa::K_CONST_K * T0));
   std::cout << "atom_ndens[0]: " << atom_ndens[0] << std::endl;

   for (auto mol : molecules) {
    mol_ndens.push_back(mixture.Boltzmann_distribution(T0, (1.-x_N)*tot_ndens, mol, 0));
   } 

   for ( auto i = mol_ndens.begin(); i != mol_ndens.end(); i++ ) {
    std::cout << (*i) << std::endl;
   }

// mol_ndens2m_x
  arma::vec arma_m_x_tot = arma::zeros(48); 
  std::vector<double> tmp_m_x_tot;
  for (auto i=mol_ndens.begin(); i!=mol_ndens.end(); i++) {
   for (int j=0; j<(*i).size(); ++j) {
    tmp_m_x_tot.push_back((*i).at(j));
   }
  }
  j=0;
  for (auto i=tmp_m_x_tot.begin(); i!=tmp_m_x_tot.end(); i++) {
   arma_m_x_tot[j] = (*i);
   j++;
  }
  for (j=0; j<arma_m_x_tot.size(); ++j) {
   std::cout << " in arma::vec mx " <<  arma_m_x_tot.at(j) << std::endl;
  }

  arma::vec arma_m_y_tot = arma::zeros(arma_m_x_tot.size());
  std::cout << " in arma::vec my " <<  arma_m_y_tot.size() << std::endl;
  std::cout << " NS " <<  _NS << std::endl;

  //arma_m_y_tot = mixture.convert_molar_frac_to_mass(arma_m_x_tot); FIXME doens't work in this way!
  // Molar mass
  double sum=0.;
  for (int i=0; i<arma_m_y_tot.size(); ++i) {
   sum += arma_m_x_tot[i] * molmass;
  }

// Conversion molar -> mass fractions
  for (int i=0; i<arma_m_y_tot.size(); ++i) {
   arma_m_y_tot[i] = arma_m_x_tot[i] * molmass / sum;
  }
 
  for (j=0; j<arma_m_y_tot.size(); ++j) {
   std::cout << " in arma::vec my " <<  arma_m_y_tot.at(j) << std::endl;
  }
  
  double density = mixture.compute_density(mol_ndens, atom_ndens, 0);

  std::cout << " density = " << density << " at temperature = " << T0 << " and pressure = " << p0 << std::endl;

  for (j=0; j<arma_m_y_tot.size(); ++j) {
   std::cout << " rho_ci =  " <<  arma_m_y_tot.at(j)*density << std::endl;
  }

  return 0;
  }
