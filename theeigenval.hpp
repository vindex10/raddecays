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

#include "theeigenval.tpp"
