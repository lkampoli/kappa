/* 
 * Test for multi-component diffusion coefficients' computation in the StS approach.
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
   
   std::cout << "Start test: computation of thermal diffusion coefficients" << std::endl;
   
   std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
   std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
   std::string particle_source    = m_source + "particles.yaml";
   std::string interaction_source = m_source + "interaction.yaml";
   std::string output_dir = GetCurrentWorkingDir();
   std::cout << "Current directory is: " << output_dir << std::endl;

   std::cout << "Loading particles data" << std::endl;

   // N2 molecule (non-rigid model)
   kappa::Molecule mol("N2", true, false, particle_source);
  
   // N atom
   kappa::Atom at("N", particle_source);

   std::cout << "Finished loading particles data" << std::endl;

   std::cout << "Molecule's name " << mol.name << std::endl;
   std::cout << "Atom's name " << at.name << std::endl;
   std::cout << "Molecule vibrational levels " << mol.num_vibr_levels[0] << std::endl;

   // N2/N binary mixture creation
   kappa::Mixture mixture(mol, at, interaction_source, particle_source);
   //kappa::Mixture mixture(mol, interaction_source, particle_source);

   // some check print
   std::cout << "particles: " << mixture.get_n_particles() << std::endl;
   std::cout << "names: " << mixture.get_names() << std::endl;

   // set a range for temperature
   //std::vector<double> T_vals = {2500.0, 5000.0, 20000.0, 50000.0};
   std::vector<double> T_vals = {10000.0};
   //std::vector<double> T_vals = {15000.0};
   //std::vector<double> T_vals = {30000.0};

   //std::vector<double> T_vals = {1.1986210569919167e+04};
   //double p = 2.6298517223967945e+04; //101325.;
   double p = 101325.;

   // vibrational levels
   int i;

//   for (i = 0; i < 80; i++) { // assume a max num. of vibr. levels a priori
//    T_vals.push_back(500 + i * 500);
//   }

   // set an arbitrary atom mass fraction (0.20)
   int x_atom_perc = 0.;
   double x_atom = x_atom_perc / 100.;

   // arma vector for atom number density
   arma::vec atom_ndens(1);

   // vector of arma vector for molecular number density
   std::vector<arma::vec> mol_ndens;

   double tot_ndens;
   tot_ndens =  p / (kappa::K_CONST_K * T_vals[0]);
   std::cout << " tot_ndens " << tot_ndens << std::endl;

   mol_ndens.push_back(mixture.Boltzmann_distribution(T_vals[0], (1. - x_atom) * tot_ndens, mol));

   //arma::vec arma_tmp = arma::zeros(46);
   //arma_tmp[0]=1.6903014643887344e+22;
   //arma_tmp[1]=2.1296668779717096e+21;
   //arma_tmp[2]=2.7931992082665574e+20;
   //arma_tmp[3]=3.8171814412708643e+19;
   //arma_tmp[4]=5.4404393924037304e+18;
   //arma_tmp[5]=8.0936179383704602e+17;
   //arma_tmp[6]=1.2577471220841520e+17;
   //arma_tmp[7]=2.0429620765342120e+16;
   //arma_tmp[8]=3.4702684112498175e+15;
   //arma_tmp[9]=6.1668407878582462e+14;
   //arma_tmp[10]=1.1467277652510819e+14;
   //arma_tmp[11]=2.2314976367988434e+13;
   //arma_tmp[12]=4.5441448115914492e+12;
   //arma_tmp[13]=9.6816910169120422e+11;
   //arma_tmp[14]=2.1575597049658203e+11;
   //arma_tmp[15]=5.0269596200232300e+10;
   //arma_tmp[16]=1.2239137258029718e+10;
   //arma_tmp[17]=3.1119258248120232e+09;
   //arma_tmp[18]=8.2572508940577209e+08;
   //arma_tmp[19]=2.2847565747377831e+08;
   //arma_tmp[20]=6.5874389233444184e+07;
   //arma_tmp[21]=1.9780155711133137e+07;
   //arma_tmp[22]=6.1874493264987366e+06;
   //arma_tmp[23]=2.0226649733369094e+06;
   //arma_tmp[24]=6.9905592968323291e+05;
   //arma_tmp[25]=2.6413774918950797e+05;
   //arma_tmp[26]=1.1746206619806583e+05;
   //arma_tmp[27]=6.7615803677532967e+04;
   //arma_tmp[28]=5.1384869660211079e+04;
   //arma_tmp[29]=4.7140368605032192e+04;
   //arma_tmp[30]=4.7201724077435814e+04;
   //arma_tmp[31]=4.8722924260764979e+04;
   //arma_tmp[32]=5.0590309726684478e+04;
   //arma_tmp[33]=5.2334346783359055e+04;
   //arma_tmp[34]=5.3737047904416693e+04;
   //arma_tmp[35]=5.4686660586866346e+04;
   //arma_tmp[36]=5.5122481433952518e+04;
   //arma_tmp[37]=5.5013096318090240e+04;
   //arma_tmp[38]=5.4347018745491623e+04;
   //arma_tmp[39]=5.3127741558997302e+04;
   //arma_tmp[40]=5.1370128760917847e+04;
   //arma_tmp[41]=4.9096978320252056e+04;
   //arma_tmp[42]=4.6335341881050022e+04;
   //arma_tmp[43]=4.3112484070257960e+04;
   //arma_tmp[44]=3.9451465101371847e+04;
   //arma_tmp[45]=3.5366337174221044e+04;

//   for (int mo=0; mo<arma_tmp.size(); mo++) {
//    std::cout << " mo_ndens " << arma_tmp.at(mo) << std::endl;
//   }
//   for (auto i=mol_ndens.begin(); i!=mol_ndens.end(); ++i) {
//    std::cout << " mol_ndens " << (*i) << std::endl;
//   }

//   mol_ndens.push_back(arma_tmp);

// atom number density
   atom_ndens = x_atom * tot_ndens;
//   atom_ndens[0] = 4.0571848823129030e+17;


   for (int at=0; at<atom_ndens.size(); at++) {
    std::cout << " atom_ndens " << atom_ndens.at(at) << std::endl;
   }

   std::ofstream outf;
   outf.open(output_dir + "/diff_" + mol.name + "_" + at.name + "_xat" + std::to_string(x_atom_perc) + ".txt");

   // output files' header
   outf << std::setw(20) << "Temperature [K]"; 
   outf << std::setw(20) << "Diff. coeffs."; 
   outf << std::endl;


   std::cout << std::setw(20) << "Temperature [K]" << std::endl;

// main loop on temperatures
//   for (auto T : T_vals) {

    //std::cout << std::setw(20) << T_vals[0] << std::endl;
    //std::cout << std::setw(20) << T << std::endl;
    //outf << std::setw(20) << T_vals[0] << std::endl;
    //outf << std::setw(20) << T << std::endl;

//   tot_ndens =  101325.0 / (kappa::K_CONST_K * T);
//   mol_ndens[0] = mixture.Boltzmann_distribution(T, (1 - x_atom) * tot_ndens, mol);
//   atom_ndens[0] = x_atom * tot_ndens;


    mixture.compute_transport_coefficients(T_vals[0], mol_ndens, atom_ndens);
    //mixture.compute_transport_coefficients(T, mol_ndens, atom_ndens);
    mixture.get_diffusion();
    //mixture.get_lite_diffusion();
    
    // simplified binary diffusion coefficients (for checking)
    // mixture.binary_diffusion(T_vals[0]);
    //mixture.binary_diffusion(T);
    
 //  }

   //outf.close();

   return 0;
}
