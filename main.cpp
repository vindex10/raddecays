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

class Vlin {
public:
    double b, alphaS;

    fldvar operator() (const fldarr &u, const double r) {
        return -4/3*alphaS/r + b*r;
    }
};

class Vscr {
public:
    double b, mu, alphaS;

    fldvar operator() (const fldarr &u, const double r) {
        return -4/3*alphaS + b/mu*(1 - exp(-mu*r));
    }
};

template <class Vtype>
class eq_radial_shrod {
public:
    double rC, muR, E, L;
    Vtype V;

    void operator() (const fldarr &u, fldarr &dudr, const double r) {
        double ruse = r>rC ? r : rC;

        dudr[0] = u[1];
        dudr[1] = -2*muR*(E - V(u, ruse) - L*(L+1)/2/muR/ruse/ruse)*u[0];
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

    testEnerg<eq_radial_shrod<Vlin>> test;

    test.eq.muR = 1;
    test.eq.L = 2;
    test.eq.rC = 0.001;
    test.eq.V.alphaS = 0.5;
    test.eq.V.b = 4;
    test.threshold = 1;
    test.cutScale = 50;

    cout << test(10) << endl;

    return 0;
}
