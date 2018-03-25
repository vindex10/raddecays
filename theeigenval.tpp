#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/util/odeint_error.hpp>
#include <limits>
#include <iostream>
#include <nlopt.hpp>
#include <vector>
#include "json_types.hpp"

template <typename Eq>
double TheEigenVal<Eq>::f() {
    eq.initU(u, intstep);
    eq.initTu(stpr.prevdu, intstep);

    try {
        boost::numeric::odeint::integrate_adaptive(stpr, eq, u, intstep, cutscale, intstep);
    } catch(const boost::numeric::odeint::step_adjustment_error& e) {
        u[0] = std::numeric_limits<double>::infinity();
    }

    return std::abs(std::real(u[0])*std::exp(-std::sqrt(2*eq.env.muR*std::abs(eq.E))*cutscale));
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

template <typename Eq>
void from_json(const json &j, TheEigenVal<Eq>& p) {
    try {
        p.intstep = j.at("intstep").get<double>();
    } catch(json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }

    try {
        p.etol = j.at("etol").get<double>();
    } catch(json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }

    try {
        p.estep = j.at("estep").get<double>();
    } catch(json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }
    
    try {
        p.ewindow = j.at("ewindow").get<double>();
    } catch(json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }
    
    try {
        p.eq = j.at("eq");
    } catch(json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }
}
