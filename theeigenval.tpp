#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/util/odeint_error.hpp>
#include <limits>
#include <iostream>
#include <nlopt.hpp>
#include <vector>

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
double TheEigenVal<Eq>::f_wrap(unsigned n, const double* x, double* grad, void* data) {
    TheEigenVal<Eq>* self = static_cast<TheEigenVal<Eq>*>(data);
    self->eq.E = x[0];
    double res = self->f();

    if (grad) {
        self->eq.E += self->estep;
        grad[0] = (self->f()-res)/self->estep;
        self->eq.E -= self->estep;
    }

    return res;
}

template <typename Eq>
double TheEigenVal<Eq>::findmin() {
    nlopt::opt opt(nlopt::LN_COBYLA, 1);
    std::vector<double> lb{-ewindow+eq.E};
    std::vector<double> ub{ewindow+eq.E};
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
    opt.set_min_objective(TheEigenVal<Eq>::f_wrap, static_cast<void*>(this));
    opt.set_stopval(etol);
    std::vector<double> x{eq.E};
    double minF;

    nlopt::result result = opt.optimize(x, minF);
    eq.E = x[0];
    return minF;
}
