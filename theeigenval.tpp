#include <cmath>
#include <boost/numeric/odeint.hpp>


template <typename Eq>
double TheEigenVal<Eq>::f() {
    u[0] = 0;
    u[1] = eq.xL < 1.2 ? 1 : 0; // xL is double, so == is bad idea

    boost::numeric::odeint::integrate_adaptive(stpr, eq, u, 0., cutscale, intstep);

    return std::abs(std::real(u[0])*std::exp(2*eq.env.muR*eq.E*cutscale));
}

template <typename Eq>
void TheEigenVal<Eq>::findmin() {
    double prevF = f();
    double inc = etol + 1;
    while (abs(inc) > etol) {
        eq.E += estep;
        double inc = prevF*estep/(f() - prevF);
        eq.E -= estep + inc;
    }
}
