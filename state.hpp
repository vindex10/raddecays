#ifndef STATE_HPP
#define STATE_HPP

#include <vector>
#include <gsl/gsl_spline.h>

template <class Eq>
class State {
public:
    Eq eq;
    std::vector<std::vector<double> > data;
    double maxR;

    gsl_spline* spline;
    gsl_interp_accel* accel;

    static int norm_f(unsigned ndim
                 ,const double *x
                 ,void *fdata
                 ,unsigned fdim
                 ,double *fval);

    double norm();
    void renorm();

    double operator()(double r);

    void init_spline();
    State();
    State(State<Eq>&& rhs);
    State<Eq>& operator=(State<Eq>&& rhs);
    ~State();
};

template <class Eq>
void from_json(const json &j, State<Eq>& p);

#include "state.tpp"

#endif
