/**
@page installation Installation

@tableofcontents

@section prerequisites Prerequisites

All that is required is a C++ compiler, a Python interpreter and/or Perl (shipped with every UNIX operating system installation) to use
the automated install process for the third-party dependencies and a working installation of CMake (make sure that command line support is enabled if your computer is running macOS).

Before compiling KAPPA, users need to have installed all dependencies: YAML, Armadillo (which by itself relies on BLAS, LAPACK, or for best performances OpenBLAS, MKL).

However, if you would like to build this reference manual, then you will need the
[Doxygen](http://www.stack.nl/~dimitri/doxygen/download.html) tool installed.  See
the section on [building the documentation](@ref build_docs) below for more instructions.

@section basic_install Install on Unix/Linux or Mac OS X

First, download the library and unpack the source to a directory where you have write access.

Then, you must add the following environment
variables to your .bashrc (Linux) or .bash_profile (Mac) file located in your
home directory

__Linux__

     export KPP_DIRECTORY="path_to_kappa_directory"
     export KPP_DATA_DIRECTORY=$KPP_DIRECTORY/data
     export PATH=$MPP_DIRECTORY/install/bin:$PATH
     export LD_LIBRARY_PATH=$KPP_DIRECTORY/install/lib:$LD_LIBRARY_PATH

__Mac__

     export KPP_DIRECTORY="path_to_kapp_directory"
     export KPP_DATA_DIRECTORY=$KPP_DIRECTORY/data
     export PATH=$KPP_DIRECTORY/install/bin:$PATH
     export DYLD_LIBRARY_PATH=$KPP_DIRECTORY/install/lib:$DYLD_LIBRARY_PATH

Note that you must change the `"path_to_kappa_directory"` to the correct path
on your machine.

Once the environment variables are set, you must install third-party dependencies.
This can be done manually or by invoking the install script from /scritps in the
KAPPA source directory:

    cd $KPP_DIRECTORY
    ./install-kappa-deps.pl

Once the dependencies are installed then use the following commands from the root
of the KAPPA repository to install.

     mkdir build
     cd build
     cmake ..
     make -j N install

where `N` is the number of CPU units to use for the build process (e.g. 4).

@note
If MKL is already installed on the system, it may conflict with OpenBLAS. It this case, it is 
sufficient to manually modify the CMakeCache.txt file in the build directory 
leaving blank where the MKL (libmkl_rt.so) is assumed and execute `make`.

*/
