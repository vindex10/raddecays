#include <cmath>
#include "env_deng2016scr.hpp"

double EnvScr::smearedDelta (double r) {
    return std::pow(sigma/std::sqrt(M_PI), 3)*std::exp(-sigma*sigma*r*r);
}

double EnvScr::Vv (double r) {
    return -4./3*alphaS/r;
}

double EnvScr::Vscr(double r) {
    return b/mu*(1 - std::exp(-mu*r));
}

double EnvScr::dVscr(double r) {
    return b*std::exp(-mu*r);
}

double EnvScr::Vss(double r, double xS, double xS1, double xS2) {
    return 32*M_PI*alphaS/9/mC/mC*smearedDelta(r)*((xS*xS-1.) - (xS1*xS1-1.) - (xS2*xS2-1.))/8;
}

double EnvScr::Vsl(double r, double xJ, double xL, double xS)  {
    double ruse = r ? r > rC : rC;
    return 1/2/mC/mC*(4*alphaS/ruse/ruse/ruse - b/r)*((xJ*xJ-1.) - (xL*xL-1.) - (xS*xS - 1.))/8;
}
