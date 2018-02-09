#include <cmath>
#include "env_deng2016lin.hpp"

double EnvLin::smearedDelta (double r) {
    return std::pow(sigma/std::sqrt(M_PI), 3)*std::exp(-sigma*sigma*r*r);
}

double EnvLin::Vv (double r) {
    return -4./3*alphaS/r;
}

double EnvLin::Vlin(double r) {
    return b*r;
}

double EnvLin::dVlin(double r) {
    return b;
}

double EnvLin::Vss(double r, double xS, double xS1, double xS2) {
    return 32*M_PI*alphaS/9/mC/mC*smearedDelta(r)*((xS*xS-1.) - (xS1*xS1-1.) - (xS2*xS2-1.))/8;
}

double EnvLin::Vsl(double r, double xJ, double xL, double xS)  {
    double ruse = r ? r > rC : rC;
    return 1/2/mC/mC*(4*alphaS/ruse/ruse/ruse - b/r)*((xJ*xJ-1.) - (xL*xL-1.) - (xS*xS - 1.))/8;
}
