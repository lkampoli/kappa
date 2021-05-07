/*! 
    \file constants.h
    \brief Various physical constants
 */

#ifndef kappa_constants_h
#define kappa_constants_h

#include <cmath>

namespace kappa {

/// \defgroup constants
/// @{

/// \defgroup math Mathematical constants
/// @{
const double PI    = 4.0 * std::atan(1.0); //!< \f$\pi\f$
const double TWOPI = 2.0 * PI;             //!< \f$2\pi\f$
const double SQRT2  = std::sqrt(2.0);      //!< \f$\sqrt{2}\f$
/// @}

/// \defgroup physical Physical constants (all units in SI)
/// @{
const double NA     = 6.0221415E23;    //!< Avagadro's number (molecule/mol)
const double KB     = 1.3806503E-23;   //!< Boltzmann's constant (J/molecule-K)
const double RU     = NA * KB;         //!< Universal Gas constant (J/mole-K)
const double HP     = 6.626068E-34;    //!< Planck's constant (J-s)
const double MU0    = PI * 4.0E-7;     //!< Magnetic constant (H/m)
const double C0     = 299792458.0;     //!< Speed of light in vacuum (m/s)
const double EPS0   = 1.0/(MU0*C0*C0); //!< Vacuum permittivity (F/m)
const double QE     = 1.602176565E-19; //!< Elementary positive charge (C)
const double ONEATM = 101325.0;        //!< 1 atm in Pa

const double K_CONST_K                 = 1.3806504e-23;                   // Boltzmann constant
const double K_CONST_C                 = 2.99792458e8;                    // Speed of light
const double K_CONST_PI                = 3.14159265358979323846264338328; // Pi
const double K_CONST_H                 = 6.62606957e-34;                  // Planck constant
const double K_CONST_HBAR              = 1.05457162825e-34;               // Planck constant, hbar
const double K_CONST_NA                = 6.02214179e23;                   // Avogadro number
const double K_CONST_ELEMENTARY_CHARGE = 1.60217646e-19;                  // Elementary charge in Coulombs
const double K_CONST_R                 = 8.3144598;                       // universal gas constant [J·mole−1·K−1]
const double K_CONST_SVV_FHO           = 0.04;                            // FHO model steric factor for VV-exchanges
const double K_CONST_EV                = 1.602176565e-19;                 // 1 electron-volt in Joules
const double K_CONST_E0                = 8.854187817e-12;                 // vacuum permittivity
const double K_CONST_EULER             = 0.577216;                        // Euler's constant

const int K_CONST_SUBDIVISIONS = 5;       // amount of subdivisions for numerical integration over semi-infinite intervals
const int K_CONST_OMEGA_D_STEP_SIZE = 10; // temperature step size for evaluating derivatives when calculating omega integrals of higher orders
/// @}
/// @}

} // namespace kappa

#endif /* constants_h */
