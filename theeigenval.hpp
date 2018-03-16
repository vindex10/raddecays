#ifndef THEEIGENVAL_HPP
#define THEEIGENVAL_HPP

#include <nlohmann/json.hpp>
#include "types.hpp"

template <typename Eq>
class TheEigenVal {
public:
    Eq eq;
    fldarr u;
    double cutscale, intstep;
    double etol, estep, ewindow;
    stepper stpr;

    // coefficient near positive exponent on u-asymptotics
    double f();
    static double f_wrap(unsigned n, const double* x, double* grad, void* data);
 
    // find minimum of coefficient with gradient descent
    double findmin();
};

template <typename Eq>
void from_json(const nlohmann::json &j, TheEigenVal<Eq>& p);

#include "theeigenval.tpp"

#endif
