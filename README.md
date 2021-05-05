# KAPPA
Kinetic Approach to Physical Processes in Atmospheres is an open-source library 
written in C++ and developed at the Department of Hydroaeromechanics of the 
Saint Petersburg State University (SPBSU), designed to be coupled with 
conventional CFD codes to provide thermodynamic, transport, chemistry, and 
energy transfer properties associated with non-equilibrium reacting flows.

[![GitHub Workflow Status](https://img.shields.io/github/workflow/status/lkampoli/kappa/CMake?style=plastic)](https://github.com/lkampoli/kappa/actions)
[![license: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![GitHub forks](https://img.shields.io/github/forks/lkampoli/kappa?style=plastic)](https://github.com/lkampoli/kappa/watchers)

The main features of _KAPPA_ code are the following:
* Kinetic theory approximations implemented:
    - Vibrational state-to-state (STS);
    - 1-temperature;
    - multi-temperature.
* The tree of directories reflects the class structure in /src:
    - particles;
    - mixtures:
    - approximations;
    - interactions;
    - numerics.
* Database files for particle and interaction are in /data:
    - particles.yaml;
    - interaction.yaml.
* Example testcases are located in /tests.

## Installation
Detailed information about installation can be found at:
https://github.com/lkampoli/kappa/wiki/HOWTO.

## Copyrights

KAPPA is an open source project, it is distributed under the 
[LGPL v3](https://www.gnu.org/licenses/lgpl-3.0.en.html). Anyone interested in 
using, developing or contributing to KAPPA is welcome. Take a look at the 
[contributing guidelines](CONTRIBUTING.md) to start to contribute to the 
project.

## Contributors

A list of contributors can found in [CONTRIBUTORS.md](CONTRIBUTORS.md).

## Source code structure

* Files in /src/approximations implement the various (state-to-state, multi- and one-temperature) kinetic theory approximations:
    - approximation.hpp, approximation_multit.hpp and approximation_onet.hpp are the header files for the state-to-state, multi-temperature and one-temperatue approximations, correspdondingly;
    - approximations.hpp is a include file containing includes for approximation.hpp, approximation_multit.hpp and approximation_onet.hpp
    - exceptions.hpp defines exceptions which can be thrown by KAPPA
    - kappa.hpp is the include file intented for use in new code and includes headers providing all functionaly of KAPPA
    - models.h defines the strongly typed enums of various models used by KAPPA (for calculation of rates, cross-sections, etc.)
    - .cpp files are the implementation files for the header files listed above
* Files in /src/interactions implement the Interaction class:
    - interactions.hpp is the header file for the Interaction class
    - interactions.cpp files is the implementation file for the Interaction class
* Files in /src/mixtures implement the Mixture class (used for state-to-state computations of thermodynamic properties and transport coefficients in mixtures):
    - mixture.hpp is the header file for the Mixture class
    - mixture.cpp files is the implementation file for the Mixture class
* Files in /src/numerics implement various numerical routines (integration over finite and infinite intervals, factorial calculation, etc):
    - numeric.hpp is the header file containing definitions of various numerical methods
    - numeric.cpp is the implementation file of the various numerical routines
    - constants.h is a header file containing definitions of various mathematical and physical constants
* Files in /src/particles implement the Particle class and its derivative classes (Atom, Molecule):
    - particle.hpp, atom.hpp and molecule.hpp are the header files for the Particle, Atom and Molecule classes
    - particle.cpp, atom.cpp and molecule.cpp are the implementation files for the Particle, Atom and Molecule classes
    
## Reference 
https://www.researchgate.net/publication/328706345_KAPPA_Kinetic_approach_to_physical_processes_in_atmospheres_library_in_C
