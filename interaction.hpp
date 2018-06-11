#ifndef INTERACTION_HPP
#define INTERACTION_HPP

#include <cmath>
#include <complex>
#include <gsl/gsl_sf_coupling.h>
#include "cubature.h"
#include "state.hpp"

template <class Eq>
class Interaction {
public:
    State<Eq> instate, outstate;
    double alphaEM;

    double coefQ(double xL, double xJ, double xlam, double xjf, double xji);
    double coefC(double xL, double xJ, double xlam, double xjf, double xji);

    struct melParamBundle {
        Interaction<Eq>* obj;
        double coefC1, coefC2, coefQ1, coefQ2;
        double xJ, xlam;
    };

    static int melMxJ_f(unsigned ndim
                 ,const double *x
                 ,void *fdata
                 ,unsigned fdim
                 ,double *fval);
    std::complex<double> melMxJ(double xJ, double xlam, double xjf, double xji);
    static int melExJ_f(unsigned ndim
                 ,const double *x
                 ,void *fdata
                 ,unsigned fdim
                 ,double *fval);
    std::complex<double> melExJ(double xJ, double xlam, double xjf, double xji);
    static int melMLW_f(unsigned ndim
                 ,const double *x
                 ,void *fdata
                 ,unsigned fdim
                 ,double *fval);
    std::complex<double> melMLW(double xlam, double xjf, double xji);
    static int melELW_f(unsigned ndim
                 ,const double *x
                 ,void *fdata
                 ,unsigned fdim
                 ,double *fval);
    std::complex<double> melELW(double xlam, double xjf, double xji);


    double widthMel(std::complex<double> mel);

    double widthMxJ(double xJ, double xlam, double xjf, double xji);
    double widthMxJ(double xJ, double xlam, double xj, bool subthr=false);
    double widthMxJ(double xJ, double xlam);
    double widthMxJ(double xJ);
    double widthMxJHel(double xJ, double xH, bool subthr=false);

    double widthExJ(double xJ, double xlam, double xjf, double xji);
    double widthExJ(double xJ, double xlam, double xj, bool subthr=false);
    double widthExJ(double xJ, double xlam);
    double widthExJ(double xJ);
    double widthExJHel(double xJ, double xH, bool subthr=false);

    double widthELW(double xlam, double xjf, double xji);
    double widthELW(double xlam, double xj, bool subthr=false);
    double widthELW(double xlam);
    double widthELW();
    double widthELWHel(double xH, bool subthr=false);

    double width(double xJ, double xlam, double xjf, double xji);
    double width(double xJ, double xlam, double xj, bool subthr=false);
    double width(double xJ, double xlam);
    double width(double xJ);
    double widthHel(double xJ, double xH, bool subthr=false);

    double reduceWidth(double width);
};

#include "interaction.tpp"

#endif
