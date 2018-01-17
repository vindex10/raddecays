#include <complex>
#include <cmath>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gamma.h>

using namespace std::complex_literals;
using std::complex;

complex<double> sphHarm(int l, int m, double theta, double phi) {
    return gsl_sf_legendre_sphPlm(l, (m<0?-m:m), cos(theta))*exp(1i*(double)m*phi);
}
