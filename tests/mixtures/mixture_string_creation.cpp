
#include "kappa.hpp"
#include <iostream>
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

    std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
    std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
    std::string particle_source    = m_source + "particles.yaml";
    std::string interaction_source = m_source + "interaction.yaml";
    std::string output_dir = GetCurrentWorkingDir();
    std::cout << "Current directory is: " << output_dir << std::endl;

//  kappa::Mixture N2N("N2, N", particle_source, interaction_source);
//  kappa::Mixture N2N("N2, N", interaction_source, particle_source);
    kappa::Mixture N2N("N2, N", interaction_source, particle_source);

    std::cout << N2N.get_names();

}
