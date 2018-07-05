// \file numeric.cpp

#include "numeric.hpp"
#include "constants.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // integrate the function f on the semi-infinite interval [a, +infinity)
  // it performs a change of variables: x = a + (1-t)/t and then integrates over the interval [0,1]
  // the interval is split into equal-sized subintervals and the integral is computed over each
  // subinterval using the 15-point Kronrod rule
  double kappa::integrate_semi_inf(std::function<double(double) > f, double a, int subdivisions, double* error_estimate) {

    double result = 0., tmp_error;
    double step_size = 1. / subdivisions;

    auto f_transform = [f, a](double t) {
      return f(a + (1 - t) / t) / (t * t);
    };

    for (int i = 0; i < subdivisions; i++) {
      if (error_estimate == nullptr) {
        result += kappa::integrate_interval(f_transform, i*step_size, (i + 1) * step_size);
      } else {
        result += kappa::integrate_interval(f_transform, i*step_size, (i + 1) * step_size, &tmp_error);
        *error_estimate += tmp_error;
      }
    }
    return result;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // integrate the function f
  // on the finite interval [a, b], where b > a, using a 15-point Kronrod rule
  double kappa::integrate_interval(std::function<double (double) > f, double a, double b, double* error_estimate) {

    double result;
    double bma = (b - a) / 2., bpa = (a + b) / 2.;

    result  = (f(-0.991455371120813 * bma + bpa) + f(0.991455371120813 * bma + bpa)) * 0.022935322010529;
    result += (f(-0.949107912342759 * bma + bpa) + f(0.949107912342759 * bma + bpa)) * 0.063092092629979;
    result += (f(-0.864864423359769 * bma + bpa) + f(0.864864423359769 * bma + bpa)) * 0.104790010322250;
    result += (f(-0.741531185599394 * bma + bpa) + f(0.741531185599394 * bma + bpa)) * 0.140653259715525;
    result += (f(-0.586087235467691 * bma + bpa) + f(0.586087235467691 * bma + bpa)) * 0.169004726639267;
    result += (f(-0.405845151377397 * bma + bpa) + f(0.405845151377397 * bma + bpa)) * 0.190350578064785;
    result += (f(-0.207784955007898 * bma + bpa) + f(0.207784955007898 * bma + bpa)) * 0.204432940075298;
    result += f(0. + bpa) * 0.209482141084728;

    result *= bma;

    if (error_estimate != nullptr) {
      *error_estimate  = (f(-0.949107912342759 * bma + bpa) + f(0.949107912342759 * bma + bpa)) * 0.129484966168870;
      *error_estimate += (f(-0.741531185599394 * bma + bpa) + f(0.741531185599394 * bma + bpa)) * 0.279705391489277;
      *error_estimate += (f(-0.405845151377397 * bma + bpa) + f(0.405845151377397 * bma + bpa)) * 0.381830050505119;
      *error_estimate += f(0. + bpa) * 0.417959183673469;
      *error_estimate *= bma;
      *error_estimate = pow(200. * fabs(*error_estimate - result), 1.5);
    }
    return result;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int kappa::find_max_value(std::function<double(int) > f, double max_value, int start) {

    for (int i = start; i < INT_MAX-1; i++) {
      if ((f(i) < max_value) && (f(i + 1) > max_value)) {
	return i;
      }
    }
    return -1;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::factorial(int n) {
  
    double res = 1.0;
    for (int i = 2; i <= n; i++) {
      res *= i;
    }
    return res;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::fact_div_fact(int start, int end) {

    double res = 1.0;
    for (int i = start + 1; i <= end; i++) {
      res *= i;
    }
    return res;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double kappa::convert_cm_to_Joule(double x) { 
    return K_CONST_H * K_CONST_C * 100 * x;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  const extern arma::vec::fixed<70> kappa::factorial_table = kappa::compute_factorial_table();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  arma::vec::fixed<70> kappa::compute_factorial_table(void) {

    arma::vec::fixed<70> result;
    for (int i = 0; i < 70; i++) {
      result[i] = kappa::factorial(i);
    }
    return result;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
