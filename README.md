# KAPPA
Kinetic Approach to Physical Processes in Atmospheres is an open-source library written in C++ and developed at the Department of Hydroaeromechanics of the Saint Petersburg State University (SPBSU), designed to be coupled with conventional CFD codes to provide thermodynamic, transport, chemistry, and energy transfer properties associated with non-equilibrium reacting flows.

The main features of _KAPPA_ code are the following:
* Models implemented:
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

## Copyrights

KAPPA is an open source project, it is distributed under the [LGPL v3](https://www.gnu.org/licenses/lgpl-3.0.en.html). Anyone interested in using, developing or contributing to KAPPA is welcome. Take a look at the [contributing guidelines](CONTRIBUTING.md) to start to contribute to the project.

## Contributors

A list of contributors to the project can found in [CONTRIBUTORS.md](CONTRIBUTORS.md)
