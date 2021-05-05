Dependencies 
************

### **KAPPA depends on the following packages:**

[YAML](https://github.com/jbeder/yaml-cpp)

[OpenBLAS](https://github.com/xianyi/OpenBLAS)

[Armadillo](https://gitlab.com/conradsnicta/armadillo-code.git)

---

KAPPA can be compiled
* on Linux using GNU or Intel compilers;

It used to compile on Windows too, but needs to be verified ...

---

**NOTE:** if you don't want to manually install all or some of the dependencies, a script is available for automatically install them (install-kappa-deps.pl).

In order to use the installation script, the user can go into the folder scripts/ and run:

`./install-kappa-deps.pl`

This will install all default dependencies (yaml, amradillo, openblas) in the default location $HOME/kappa/kappa-deps/$ARCH (inside bin, lib, include subfolders), after having compiled all packages inside $HOME/kappa/kappa_tmp.

---

**NOTE:** once all KAPPA dependencies have been installed, the user must set PATH and LD_LIBRARY_PATH variables to point to the right locations.

---

### UPDATE -- 14/05/2020
The manual installation of each third-party dependency through the script `./install-kappa-deps.pl` may be by-passed, at least on Linux machines, by directly installing the packages from their repositores with the following line:

`sudo apt-get install libopenblas-dev libarmadillo-dev libyaml-cpp-dev`
