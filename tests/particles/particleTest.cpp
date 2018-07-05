/* 
 * \file particleTest.cpp
 * \author aspera
 *
 * Created on 4 марта 2016 г., 12:36
 */

#include <cstdlib>
#include <iostream>

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

using namespace std;
using namespace kappa;

int main(int argc, char** argv) {

  cout << "Start Test for Particle classes" << endl;

  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
  std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
  std::string particle_source    = m_source + "particles.yaml";
  std::string interaction_source = m_source + "interaction.yaml";

  try {
    Particle Ar_badfilename("Ar","articles.yaml");  
    // test that trying to read a non-existent database file throws the correct error
  } catch (const UnopenedFileException& e) {
    std::cout << e.what() << endl;
  }
   
  try {
    // Particle Bad("B");  // Bug !
    Particle Bad("B", particle_source);  
    // test that trying to load particles that are not in the database throw the correct error
  } catch (const DataNotFoundException& e) {
    std::cout << e.what() << endl;
  }

  Particle Ar("Ar", particle_source);
  Particle C("C", particle_source);
  Particle e("e-", particle_source);
  Particle N("N", particle_source);
  Particle Nplus("N+", particle_source);
  Particle O("O", particle_source);
  Particle Oplus("O+", particle_source);
  Particle Ominus("O-", particle_source);
   
  Atom Ar_atom("Ar", particle_source);
  Atom C_atom("C", particle_source);
  Atom N_atom("N", particle_source);
  Atom Nplus_atom("N+", particle_source);
  Atom O_atom("O", particle_source);
  Atom Oplus_atom("O+", particle_source);
  Atom Ominus_atom("O-", particle_source);
    
  // Molecule C2("C2","anfalse",false, particle_source);
  Molecule C2("C2",true,false, particle_source);
  Molecule CO("CO",true,false, particle_source);
  Molecule N2("N2",true,false, particle_source);
  Molecule N2plus("N2+", true,false, particle_source);
  Molecule NO("NO", true,false, particle_source);
  Molecule O2("O2", true,false, particle_source);
  Molecule O2plus("O2+", true,false, particle_source);

  std::cout << "Number of electron levels in CO: " << CO.num_electron_levels << endl;  // should be 28
  //std::cout << "Vibrational spectrum of CO: " << CO.vibr_spectrum << endl;  // should be "anfalse"
  std::cout << "Vibrational spectrum of CO: " << CO.anharmonic_spectrum << endl;  // should be true
  std::cout << "Rigid rotator model used for CO molecule: " << CO.rigid_rotator << endl; // should be false (0)
  std::cout << "N2 number of vibrational levels in ground state: " << N2.num_vibr_levels[0] << endl;  // should be 48

  // try {
  //   Molecule CObad1("CO","nhrmonic",false,"particles.yaml");  // test that giving an incorrect value of the spectrum parameter throws the correct error
  // } catch (const ModelParameterException& e) {
  //   std::cout << e.what() << endl;
  // }

  try {
    Molecule CObad2("CO",true,false,"articles.yaml");  // test that trying to read a non-existent database file throws the correct error
  } catch (const UnopenedFileException& e) {
    std::cout << e.what() << endl;
  }

  Molecule C2_h("C2","false",false, particle_source);
  Molecule CO_h("CO","false",false, particle_source);
  Molecule N2_h("N2","false", false, particle_source);
  Molecule N2plus_h("N2+","false", false, particle_source);
  Molecule NO_h("NO","false", false, particle_source);
  Molecule O2_h("O2","false", false, particle_source);
  Molecule O2plus_h("O2+","false", false, particle_source);

  std::cout << "N2 number of vibrational levels in ground state: " << N2_h.num_vibr_levels[0] << endl;  // should be 33
  
  Atom wrong_at("N2", particle_source);
  std::cout << "loaded N2 as atom";

  // Molecule wrong_mol("N", true, false, particle_source);
  // std::cout << "loaded N as molecule";

  return 0;
}
