# KAPPA
Kinetic Approach to Physical Processes in Atmospheres is an open-source library 
written in C++ and developed at the Department of Hydroaeromechanics of the 
Saint Petersburg State University (SPBSU), designed to be coupled with 
conventional CFD codes to provide thermodynamic, transport, chemistry, and 
energy transfer properties associated with non-equilibrium reacting flows.

[![GitHub Workflow Status](https://img.shields.io/github/workflow/status/lkampoli/kappa/CMake?style=plastic)](https://github.com/lkampoli/kappa/actions)
[![license: LGPL v3](https://img.shields.io/github/license/lkampoli/kappa?color=orange&style=plastic)](https://www.gnu.org/licenses/lgpl-3.0)
[![GitHub forks](https://img.shields.io/github/forks/lkampoli/kappa?style=plastic)](https://github.com/lkampoli/kappa/network/members)
[![GitHub Repo stars](https://img.shields.io/github/stars/lkampoli/kappa?color=yellow&style=plastic)](https://github.com/lkampoli/kappa/stargazers)
[![GitHub watchers](https://img.shields.io/github/watchers/lkampoli/kappa?color=green&style=plastic)](https://github.com/lkampoli/kappa/watchers)
[![GitHub commit activity](https://img.shields.io/github/commit-activity/m/lkampoli/kappa?color=red&style=plastic)](https://github.com/lkampoli/kappa/graphs/commit-activity)
[![GitHub size](https://img.shields.io/github/languages/code-size/lkampoli/kappa?color=violet&style=plastic)]()
[![GitHub lines](https://img.shields.io/tokei/lines/github/lkampoli/kappa?color=pink&style=plastic)]()
<!-- [![GitHub language count](https://img.shields.io/github/languages/count/lkampoli/kappa?color=cyan&style=plastic)]() -->

The main features of KAPPA code are the following:
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
* Example testcases are located in /examples.

## Installation
Detailed information about installation can be found at:
https://github.com/lkampoli/kappa/wiki/HOWTO.

or briefly, 

`sudo apt-get install libopenblas-dev libarmadillo-dev libyaml-cpp-dev catch`

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

## Copyrights

KAPPA is an open source project, it is distributed under the 
[LGPL v3](https://www.gnu.org/licenses/lgpl-3.0.en.html). Anyone interested in 
using, developing or contributing to KAPPA is welcome. Take a look at the 
[contributing guidelines](CONTRIBUTING.md) to start to contribute to the 
project.

## Contributors

A list of contributors can found in [CONTRIBUTORS.md](CONTRIBUTORS.md).
    
## Citing KAPPA
Please cite the following article when mentioning KAPPA in your own papers.

* Campoli, L., Oblapenko, G. P., & Kustova, E. V. (2019). [Kappa: Kinetic approach to physical processes in atmospheres library in C++.](https://doi.org/10.1016/j.cpc.2018.10.016), Computer Physics Communications, 236, 244-267.

**Bibtex**
```bibtex
@article{campoli2019kappa,
  title={Kappa: Kinetic approach to physical processes in atmospheres library in C++},
  author={Campoli, Lorenzo and Oblapenko, GP and Kustova, Elena V},
  journal={Computer Physics Communications},
  volume={236},
  pages={244--267},
  year={2019},
  publisher={Elsevier}
}
