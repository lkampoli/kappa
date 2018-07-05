/*!
 *  \file atom.hpp
 *  \brief The Atom class loads and stores spectroscopic data for atoms
 */

#ifndef atom_hpp
#define atom_hpp

#include <map>
#include <string>
#include <vector>
#include "constants.h"
#include "particle.hpp"
#include <armadillo>

namespace kappa {

class Atom : public Particle {

 public:

  //! Creates an Atom-type object by loading the atom data from the database
  //! @param name atom's name
  //! @param filename path to the database file
  Atom(const std::string &name, const std::string &filename = "particles.yaml");

 private:
        
  void readData(const std::string &name, const std::string &filename);

};
}
#endif 
