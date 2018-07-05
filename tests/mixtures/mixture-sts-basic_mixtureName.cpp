
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
    std::string _mixtureName = "air5";		//OK!
//  std::string _mixtureName = "air11"; 	//NO+ not in DB
//  std::string _mixtureName = "N2";		//OK!
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

//    instantiates molecular density vector
      std::vector<arma::vec> mol_ndens;
 
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
         kappa::Molecule mol(_speciesNames[i], true, true, particle_source);
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

    for (auto at: atoms) {
        for (auto mo: molecules) {
            std::cout << at.name << "+" << mo.name << ", interaction=" 
                      << mixture.interaction(at, mo).particle1_name << "+" 
	   	        << mixture.interaction(at, mo).particle2_name << std::endl;
        }
    }
    std::cout << " " << std::endl;

    for (auto at: atoms) {
        for (auto at2: atoms) {
            std::cout << at.name << "+" << at2.name << ", interaction=" 
		      << mixture.interaction(at, at2).particle1_name << " " 
   		      << mixture.interaction(at, at2).particle2_name << std::endl;
        }
    }
    std::cout << " " << std::endl;

    for (auto mo: molecules) {
        for (auto mo2: molecules) {
            std::cout << mo.name << "+" << mo2.name << ", interaction=" 
		      << mixture.interaction(mo, mo2).particle1_name << " " 
		      << mixture.interaction(mo, mo2).particle2_name << std::endl;
        }
    }
    std::cout << " " << std::endl;

    std::cout << std::setw(20) << "Temperature [K]";
    std::cout << std::setw(20) << "Boltzmann distr.";
    std::cout << std::endl;

   for (auto mol : molecules) {
     mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], 101325.0 / (kappa::K_CONST_K * T_vals[0]), mol, 0));
//    std::cout << std::setw(20) << mol_ndens[0][0] / (101325.0 / (kappa::K_CONST_K * T_vals[0])) << std::setw(20) <<  mol.name ;
//    std::cout << std::endl;
   }

   for ( int i = 0; i < mol_ndens.size(); ++i ) {
    std::cout << std::setw(20) << T_vals[0];
    std::cout << std::setw(20) << mol_ndens[i][0] / (101325.0 / (kappa::K_CONST_K * T_vals[0]));
    std::cout << std::endl;
   }

    std::cout << " " << std::endl;
    std::cout << std::setw(20) << "Temperature [K]";
    std::cout << std::setw(20) << "Boltzmann distr.";
    std::cout << std::endl;

    for (auto mol : molecules) {
      std::cout << " " << std::endl;
      std::cout << "Molecule: " << mol.name << std::endl;
     for (auto T : T_vals) {
      mol_ndens[0] = mixture.Boltzmann_distribution(T, 101325.0 / (kappa::K_CONST_K * T), mol, 0);
      std::cout << std::setw(20) << T;
      std::cout << std::setw(20) << mol_ndens[0][0] / (101325.0 / (kappa::K_CONST_K * T));
      std::cout << std::endl;
     }
    }

//    std::cout << " " << std::endl;
//    for (auto T : T_vals) {
//     std::cout << std::setw(20) << T;
//     std::cout << std::setw(20) << mol_ndens[0][0] / (101325.0 / (kappa::K_CONST_K * T));
//     std::cout << std::setw(20) << mol_ndens[1][0] / (101325.0 / (kappa::K_CONST_K * T)); 
//     std::cout << std::setw(20) << mol_ndens[2][0] / (101325.0 / (kappa::K_CONST_K * T)); 
//     std::cout << std::endl;
//    }
    return 0;
  }
