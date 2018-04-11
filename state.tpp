#include <vector>
#include <cmath>
#include <gsl/gsl_spline.h>
#include "cubature.h"
#include "json_types.hpp"

template <class Eq>
int State<Eq>::norm_f(unsigned ndim
                     ,const double *x
                     ,void *fdata
                     ,unsigned fdim
                     ,double *fval) {
    *fval = (*static_cast<State<Eq>*>(fdata))(*x);
    *fval *= *fval;
    return 0;
}

template <class Eq>
double State<Eq>::norm() {
    double res,err;
    double minR = 0.;
    hcubature(1, norm_f, this, 1, &minR, &maxR, 0, 1E-5, 0, ERROR_INDIVIDUAL, &res, &err);
    return std::sqrt(res);
}

template <class Eq>
void State<Eq>::renorm() {
    double norm = this->norm();
    for (double& v: data[1]) {
        v /= norm;
    }
}

template <class Eq>
double State<Eq>::operator()(double r) {
    return r > maxR ? 0 : gsl_spline_eval(spline, r, accel);
}

template <class Eq>
void State<Eq>::init_spline() {
    gsl_spline_free(spline);
    spline = gsl_spline_alloc(gsl_interp_cspline, data[0].size());
    gsl_spline_init(spline, data[0].data(), data[1].data(), data[0].size());
    gsl_interp_accel_reset(accel);
}

template <class Eq>
State<Eq>::State() {
    spline = nullptr;
    accel = gsl_interp_accel_alloc();
}

template <class Eq>
State<Eq>::State(State<Eq>&& rhs) {
    gsl_spline_free(spline);
    spline = rhs.spline;
    rhs.spline = nullptr;

    gsl_interp_accel_free(accel);
    accel = rhs.accel;
    rhs.accel = nullptr;

    eq = rhs.eq;
    maxR = rhs.maxR;
}

template <class Eq>
State<Eq>& State<Eq>::operator=(State<Eq>&& rhs) {
    if (this == &rhs) {
        return *this;
    }

    gsl_spline_free(spline);
    spline = rhs.spline;
    rhs.spline = nullptr;

    gsl_interp_accel_free(accel);
    accel = rhs.accel;
    rhs.accel = nullptr;

    eq = rhs.eq;
    maxR = rhs.maxR;

    return *this;
}

template <class Eq>
State<Eq>::~State() {
    gsl_spline_free(spline);
    gsl_interp_accel_free(accel);
}

template <class Eq>
void from_json(const json &j, State<Eq>& p) {
    try {
        p.eq = j.at("eq");
    } catch(json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }

    try {
        p.maxR = j.at("rMax");
    } catch(json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }
}
