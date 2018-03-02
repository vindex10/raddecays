#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/util/odeint_error.hpp>
#include <limits>
#include <iostream>

template <typename Eq>
double TheEigenVal<Eq>::f() {
    u[0] = 0;
    u[1] = eq.xL < 1.2 ? 1 : 0; // xL is double, so == is bad idea

    try {
        boost::numeric::odeint::integrate_adaptive(stpr, eq, u, 0., cutscale, intstep);
    } catch(const boost::numeric::odeint::step_adjustment_error& e) {
        u[0] = std::numeric_limits<double>::infinity();
    }

    return std::abs(std::real(u[0])*std::exp(2*eq.env.muR*eq.E*cutscale));
}

template <typename Eq>
void TheEigenVal<Eq>::findmin() {
    double prevF = f();
    double inc = etol + 1;
    while (abs(inc) > etol &&  eq.E < 0) {
        eq.E += estep;
        //inc = prevF*estep/(f() - prevF);
        inc = (f()-prevF)/estep*1E-20;
        eq.E -= estep + inc;
        std::cout << eq.E << " ± " << inc << std::endl;
    }
}
