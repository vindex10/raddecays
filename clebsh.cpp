#include <complex>
#include <cmath>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gamma.h>

using namespace std::complex_literals;
using std::complex;

complex<double> sphHarm(int l, int m, double theta, double phi) {
    int sign = m<0 ? -1 : 1;
    double prefactor;

    if (m < 0) {
        prefactor = (m%2 ? 1 : -1)*(double)gsl_sf_fact(l-m)/gsl_sf_fact(l+m);
    } else {
        prefactor = 1;
    }

    return prefactor*gsl_sf_legendre_sphPlm(l, sign*m, cos(theta))*exp(1i*(double)m*phi);
}
