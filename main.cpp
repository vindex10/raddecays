#include <cmath>
#include <complex>
#include <array>
#include <boost/numeric/odeint.hpp>
#include <Eigen/Dense>
#include <iostream>

#include "constants.h"

using namespace std;
using namespace std::complex_literals;
using namespace boost::numeric::odeint;
using namespace Eigen;

typedef complex<double> fldvar;
typedef Array<fldvar, 2, 1> fldarr;

namespace boost { namespace numeric { namespace odeint {
template<>
struct vector_space_norm_inf<fldarr>
{
    typedef double result_type;
    double operator()( const fldarr &p ) const
    {
        return real((p.conjugate()*p).sum());
    }
};
} } }

typedef bulirsch_stoer_dense_out<fldarr
                         , double
                         , fldarr
                         , double
                         , vector_space_algebra
                         > stepper; 


class SmearedDeltaFunc {
    double sigma;

    double operator() (double r) {
        return pow(sigma/sqrt(PI), 3)*exp(-sigma*sigma*r*r);
    }
};

class VlinFunc {
public:
    double b, alphaS;

    fldvar operator() (const double r) const {
        return -4/3*alphaS/r + b*r;
    }
};

class VscrFunc {
public:
    double b, mu, alphaS;

    fldvar operator() (const double r) const {
        return -4/3*alphaS + b/mu*(1 - exp(-mu*r));
    }
};

class VssFunc {
    double mC, alphaS;
    SmearedDeltaFunc smearedDelta;

    double operator() (double r) {
        return 32*PI*alphaS/9/mC/mC*smearedDelta(r)*(xS*xS - xS1*xS1 - xS2*xS2 + 1)/8;
    }
};

class eq_radial_shrod {
public:
    double rC, muR, E, L;
    VlinFunc Vlin;
    VscrFunc Vscr;

    void operator() (const fldarr &u, fldarr &dudr, const double r) {
        double ruse = r>rC ? r : rC;

        dudr[0] = u[1];
        dudr[1] = -2*muR*(E - Vlin(ruse) - Vscr(ruse) - L*(L+1)/2/muR/ruse/ruse)*u[0];
    }
};

template <class Eq>
class testEnerg {
public:
    double threshold, cutScale;
    Eq eq;
    fldarr u;
    stepper stpr;

    bool operator() (double E) {
        eq.E = E;
        u[0] = 0.2;
        u[1] = 0.5;

        integrate_adaptive(stpr, eq, u, 0., cutScale, 0.01);
        return (abs(u[0]) < threshold);
    }
};

int main() {

    testEnerg<eq_radial_shrod> test;

    test.eq.muR = 1;
    test.eq.L = 2;
    test.eq.rC = 0.001;
    test.eq.Vlin.alphaS = 0.5;
    test.eq.Vlin.b = 4;
    test.eq.Vscr.alphaS = 0.5;
    test.eq.Vscr.b = 4;
    test.eq.Vscr.mu = 4;
    test.threshold = 1;
    test.cutScale = 50;

    cout << test(10) << endl;

    return 0;
}
