#include <cmath>
#include <complex>
#include <functional>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_bessel.h>
#include "cubature.h"
#include "utils.hpp"
#include <iostream>

using namespace std;

template <class Eq>
double Interaction<Eq>::widthMel(std::complex<double> mel) {
    double Mf = 2*outstate.eq.env.mC + outstate.eq.E;
    double Mi = 2*instate.eq.env.mC + instate.eq.E;
    /*cout << "widthMel: " << mel << ", " << alphaEM/2.*(Mi*Mi*Mi*Mi - Mf*Mf*Mf*Mf)/(Mi*Mi*Mi) << ", " << (Mi - Mf) << endl;*/
    return alphaEM/2.*(Mi*Mi*Mi*Mi - Mf*Mf*Mf*Mf)/(Mi*Mi*Mi)*std::norm(mel);
}

template <class Eq>
int Interaction<Eq>::melMxJ_f(unsigned ndim
                     ,const double *x
                     ,void *fdata
                     ,unsigned fdim
                     ,double *fval) {
    struct melParamBundle* params = static_cast<struct melParamBundle*>(fdata);
    State<Eq>* in = &(params->obj->instate);
    State<Eq>* out = &(params->obj->outstate);
    /*cout << "in MxJ" << endl;*/

    double k = 2.*out->eq.env.mC + out->eq.E;
    k = -k + std::sqrt(k*k + 2*k*(in->eq.E - out->eq.E));

    fval[0] = (*in)(*x)*(*out)(*x)*gsl_sf_bessel_jl(std::lround((params->xJ-3.)/2.), (*x)*k/2.);
    fval[1] = (*in)(*x)*(*out)(*x)*gsl_sf_bessel_jl(std::lround((params->xJ+1)/2.), (*x)*k/2.);

    return 0;
}

template <class Eq>
std::complex<double> Interaction<Eq>::melMxJ(double xJ, double xlam, double xjf, double xji) {
    if (std::lround((instate.eq.xS + outstate.eq.xS + xJ - 3.)/2.) % 2 != 0) {
        return 0.;
    }

    double* jmels = new double[2];
    double err;
    struct melParamBundle params;
    /*cout << "MxJ.  xJ: " << xJ << ", xlam: " << xlam << ", xjf: " << xjf << ", xji: " << xji << endl;*/
    params.obj = this;
    params.xJ = xJ;

    double minR = 0.;
    double maxR = std::min(instate.maxR, outstate.maxR);
    hcubature(2, melMxJ_f, &params, 1, &minR, &maxR, 0, 1E-5, 0, ERROR_INDIVIDUAL, jmels, &err);

    double k = 2.*outstate.eq.env.mC + outstate.eq.E;
    k = -k + std::sqrt(k*k + 2*k*(instate.eq.E - outstate.eq.E));

    std::complex<double> res = 0.;
    res += gsl_sf_coupling_9j(outstate.eq.xJ-1., instate.eq.xJ-1., xJ-1.
                             ,outstate.eq.xS-1., instate.eq.xS-1., 2.
                             ,outstate.eq.xL-1., instate.eq.xL-1., xJ-3.)
        *  clebsch(xJ-2., 1., instate.eq.xL, 1., outstate.eq.xL, 1.)
        *  std::sqrt((xJ+1.)/2.*(xJ-2.))
        *  jmels[0];
    res -= gsl_sf_coupling_9j(outstate.eq.xJ-1., instate.eq.xJ-1., xJ-1.
                             ,outstate.eq.xS-1., instate.eq.xS-1., 2.
                             ,outstate.eq.xL-1., instate.eq.xL-1., xJ+1.)
        *  clebsch(xJ+2., 1., instate.eq.xL, 1., outstate.eq.xL, 1.)
        *  std::sqrt((xJ-1.)/2.*(xJ+2.))
        *  jmels[1];
    res *= std::sqrt(3.)/2.*(xlam-1)/2.*k/instate.eq.env.muR
        *  clebsch(xJ, 2.-xlam, instate.eq.xJ, xji, outstate.eq.xJ, xjf)
        *  gsl_sf_coupling_6j(outstate.eq.xS-1., 2., instate.eq.xS-1.
                                          , 1., 1., 1.)
        *  std::sqrt(instate.eq.xL*instate.eq.xS*outstate.eq.xS*instate.eq.xJ*xJ);

    if (std::lround((xJ-1.)/2.+outstate.eq.xS)%2 == 0) {
        res *= (std::lround((xJ-1.)/2.+outstate.eq.xS)%4 == 0 ? 1. : -1.);
    } else {
        res *= (std::lround((xJ-1.)/2.+outstate.eq.xS - 1)%4 == 0 ? 1. : -1.)
            *  std::complex<double>(0., 1.);
    }

    return res;
}

template <class Eq>
double Interaction<Eq>::widthMxJ(double xJ, double xlam, double xjf, double xji) {
    std::complex<double> mel = 0.;
    for (int i=1; i <= std::lround((xJ-1.)/2.); ++i) {
        mel += melMxJ(2.*i+1., xlam, xjf, xji);
    }
    return widthMel(mel);
}

template <class Eq>
double Interaction<Eq>::widthMxJ(double xJ, double xlam, double xj, bool subthr) {
    double res = 0.;

    if (!subthr) {
        for (int jf_cnt = 0; jf_cnt < std::lround(outstate.eq.xJ); ++jf_cnt) {
            res += widthMxJ(xJ, xlam, -outstate.eq.xJ+2.+2.*(double)jf_cnt, xj);
        }
        return res;
    } else {
        for (int ji_cnt = 0; ji_cnt < std::lround(instate.eq.xJ); ++ji_cnt) {
            res += widthMxJ(xJ, xlam, xj, -instate.eq.xJ+2.+2.*(double)ji_cnt);
        }
        return res/instate.eq.xJ;
    }
}

template <class Eq>
double Interaction<Eq>::widthMxJ(double xJ, double xlam) {
    double res = 0.;
    for (int ji_cnt = 0; ji_cnt < std::lround(instate.eq.xJ); ++ji_cnt) {
        res += widthMxJ(xJ, xlam, -instate.eq.xJ+2.+2.*(double)ji_cnt);
    }
    return res/instate.eq.xJ;
}

template <class Eq>
double Interaction<Eq>::widthMxJ(double xJ) {
    return widthMxJ(xJ, -1.) + widthMxJ(xJ, 3.);
}

template <class Eq>
double Interaction<Eq>::widthMxJHel(double xJ, double xH, bool subthr) {
    return (widthMxJ(xJ, -1., xH, subthr) + widthMxJ(xJ, 3., xH, subthr) + (std::lround(xH) != 1 ? widthMxJ(xJ, -1., -xH+2., subthr) + widthMxJ(xJ, 3., -xH+2., subthr) : 0.))/(subthr ? 1. : instate.eq.xJ);
}

template <class Eq>
int Interaction<Eq>::melExJ_f(unsigned ndim
                     ,const double *x
                     ,void *fdata
                     ,unsigned fdim
                     ,double *fval) {
    struct melParamBundle* params = static_cast<struct melParamBundle*>(fdata);
    State<Eq>* in = &(params->obj->instate);
    State<Eq>* out = &(params->obj->outstate);

    double k = 2*out->eq.env.mC + out->eq.E;
    k = -k + std::sqrt(k*k + 2*k*(in->eq.E - out->eq.E));

    fval[0] = (*in)(*x)*(*out)(*x)*gsl_sf_bessel_jl(std::lround((params->xJ-1.)/2.), (*x)*k/2.);
    /*cout << "in ExJ" << endl;*/

    return 0;
}

template <class Eq>
std::complex<double> Interaction<Eq>::melExJ(double xJ, double xlam, double xjf, double xji) {
    double* jmel = new double[1];
    double err;
    struct melParamBundle params;
    params.obj = this;
    params.xJ = xJ;

    /*cout << "melExJ. xJ: " << xJ << ", xlam: " << xlam << ", xjf: " << xjf << ", xji: " << xji << endl;*/

    double minR = 0.;
    double maxR = std::min(instate.maxR, outstate.maxR);
    hcubature(1, melExJ_f, &params, 1, &minR, &maxR, 0, 1E-5, 0, ERROR_INDIVIDUAL, jmel, &err);

    double k = 2*outstate.eq.env.mC + outstate.eq.E;
    k = -k + std::sqrt(k*k + 2*k*(instate.eq.E - outstate.eq.E));

    std::complex<double> res = 0.;

    if (std::lround((xJ-1.)/2.)%2 != 0
            &&
        std::lround(instate.eq.xS-outstate.eq.xS) == 0) {
        res += -(std::lround((outstate.eq.xJ+instate.eq.xL-2.)/2.)%2 == 0 ? 1. : -1.)
            *  std::sqrt(2*(xJ-1.)/2.*(xJ+1.)/2)*(instate.eq.E - outstate.eq.E)
            *  gsl_sf_coupling_6j(outstate.eq.xJ-1., xJ-1., instate.eq.xJ-1.
                                ,instate.eq.xL-1., instate.eq.xS-1., outstate.eq.xL-1.);
    }

    if (std::lround((instate.eq.xS+outstate.eq.xS+xJ-3.)/2.)%2 != 0) {
        res += std::sqrt(3*xJ*instate.eq.xS*outstate.eq.xS)/2./instate.eq.env.muR
            *  k*k
            *  gsl_sf_coupling_6j(outstate.eq.xS-1., 2., instate.eq.xS-1.
                                 ,1. ,1., 1.)
            *  gsl_sf_coupling_9j(outstate.eq.xJ-1., instate.eq.xJ-1, xJ-1.
                                 ,outstate.eq.xS-1., instate.eq.xS-1., 2.
                                 ,outstate.eq.xL-1., instate.eq.xL-1., xJ-1.);
    }

    res *= jmel[0]/k*std::sqrt(instate.eq.xJ*instate.eq.xL)*xJ
        *  clebsch(xJ, 2.-xlam, instate.eq.xJ, xji, outstate.eq.xJ, xjf)
        * clebsch(xJ, 1., instate.eq.xL, 1., outstate.eq.xL, 1.);

    if (std::lround((xJ-1.)/2. + outstate.eq.xS + 1.)%2 == 0) {
        res *= (std::lround((xJ-1.)/2. + outstate.eq.xS + 1.)%4 == 0 ? 1. : -1.);
    } else {
        res *= (std::lround((xJ-1.)/2. + outstate.eq.xS)%4 == 0 ? 1. : -1.)
            *  std::complex<double>(0., 1.);
    }

    return res;
}

template <class Eq>
int Interaction<Eq>::melELW_f(unsigned ndim
                     ,const double *x
                     ,void *fdata
                     ,unsigned fdim
                     ,double *fval) {
    struct melParamBundle* params = static_cast<struct melParamBundle*>(fdata);
    State<Eq>* in = &(params->obj->instate);
    State<Eq>* out = &(params->obj->outstate);

    double k = 2*out->eq.env.mC + out->eq.E;
    k = -k + std::sqrt(k*k + 2*k*(in->eq.E - out->eq.E));

    fval[0] = (*x)*(*in)(*x)*(*out)(*x);

    return 0;
}

template <class Eq>
std::complex<double> Interaction<Eq>::melELW(double xlam, double xjf, double xji) {
    if (std::lround(instate.eq.xS-outstate.eq.xS) != 0) {
        return 0.;
    }

    double* rmel = new double[1];
    double err;
    struct melParamBundle params;
    params.obj = this;
    /*cout << "xjf: " << xjf << ", xji: " << xji << ", xlam: " << xlam << " :: " << params.coefQ << endl;*/
    /*cout << "Ef: " << outstate.eq.E << ", Ei: " << instate.eq.E << endl;*/

    double minR = 0.;
    double maxR = std::min(instate.maxR, outstate.maxR);
    hcubature(1, melELW_f, &params, 1, &minR, &maxR, 0, 1E-5, 0, ERROR_INDIVIDUAL, rmel, &err);
    /*cout << "Mel: " <<  res[0] << " + i" << res[1] << endl;*/

    std::complex<double> res = -std::complex<double>(0., 1.)*(instate.eq.E - outstate.eq.E)*std::sqrt(instate.eq.xL*instate.eq.xJ)
                             * clebsch(3., 2.-xlam, instate.eq.xJ, xji, outstate.eq.xJ, xjf)
                             * clebsch(3., 1., instate.eq.xL, 1., outstate.eq.xL, 1.)
                             * gsl_sf_coupling_6j(outstate.eq.xJ-1., 2., instate.eq.xJ-1.
                                                 ,instate.eq.xL-1., instate.eq.xS-1., outstate.eq.xL-1.)
                             * (std::lround((outstate.eq.xS + outstate.eq.xJ + instate.eq.xL -3.)/2. + 1.)%2 ? 1. : -1.)
                             * rmel[0];
    return res;
}

template <class Eq>
double Interaction<Eq>::widthExJ(double xJ, double xlam, double xjf, double xji) {
    std::complex<double> mel = 0.;
    for (int i=1; i <= std::lround((xJ-1.)/2.); ++i) {
        mel += melExJ(2.*i+1., xlam, xjf, xji);
    }
    return widthMel(mel);
}

template <class Eq>
double Interaction<Eq>::widthExJ(double xJ, double xlam, double xj, bool subthr) {
    double res = 0.;

    if (!subthr) {
        for (int jf_cnt = 0; jf_cnt < std::lround(outstate.eq.xJ); ++jf_cnt) {
            res += widthExJ(xJ, xlam, -outstate.eq.xJ+2.+2.*(double)jf_cnt, xj);
        }
        return res;
    } else {
        for (int ji_cnt = 0; ji_cnt < std::lround(instate.eq.xJ); ++ji_cnt) {
            res += widthExJ(xJ, xlam, xj, -instate.eq.xJ+2.+2.*(double)ji_cnt);
        }
        return res/instate.eq.xJ;
    }
}

template <class Eq>
double Interaction<Eq>::widthExJ(double xJ, double xlam) {
    double res = 0.;
    for (int ji_cnt = 0; ji_cnt < std::lround(instate.eq.xJ); ++ji_cnt) {
        res += widthExJ(xJ, xlam, -instate.eq.xJ+2.+2.*(double)ji_cnt);
    }
    return res/instate.eq.xJ;
}

template <class Eq>
double Interaction<Eq>::widthExJ(double xJ) {
    return widthExJ(xJ, -1.) + widthExJ(xJ, 3.);
}

template <class Eq>
double Interaction<Eq>::widthExJHel(double xJ, double xH, bool subthr) {
    return (widthExJ(xJ, -1., xH, subthr) + widthExJ(xJ, 3., xH, subthr) + (std::lround(xH) != 1 ? widthExJ(xJ, -1., -xH+2., subthr) + widthExJ(xJ, 3., -xH+2., subthr) : 0.))/(subthr ? 1. : instate.eq.xJ);
}

template <class Eq>
double Interaction<Eq>::widthELW(double xlam, double xjf, double xji) {
    std::complex<double> mel = 0.;
    mel = melELW(xlam, xjf, xji);
    return widthMel(mel);
}

template <class Eq>
double Interaction<Eq>::widthELW(double xlam, double xj, bool subthr) {
    double res = 0.;

    if (!subthr) {
        for (int jf_cnt = 0; jf_cnt < std::lround(outstate.eq.xJ); ++jf_cnt) {
            res += widthELW(xlam, -outstate.eq.xJ+2.+2.*(double)jf_cnt, xj);
        }
        return res;
    } else {
        for (int ji_cnt = 0; ji_cnt < std::lround(instate.eq.xJ); ++ji_cnt) {
            res += widthELW(xlam, xj, -instate.eq.xJ+2.+2.*(double)ji_cnt);
        }
        return res/instate.eq.xJ;
    }
}

template <class Eq>
double Interaction<Eq>::widthELW(double xlam) {
    double res = 0.;
    for (int ji_cnt = 0; ji_cnt < std::lround(instate.eq.xJ); ++ji_cnt) {
        res += widthELW(xlam, -instate.eq.xJ+2.+2.*(double)ji_cnt);
    }
    return res/instate.eq.xJ;
}

template <class Eq>
double Interaction<Eq>::widthELW() {
    return widthELW(-1.) + widthELW(3.);
}

template <class Eq>
double Interaction<Eq>::widthELWHel(double xH, bool subthr) {
    return (widthELW(-1., xH, subthr) + widthELW(3., xH, subthr) + (std::lround(xH) != 1 ? widthELW(-1., -xH+2., subthr) + widthELW(3., -xH+2., subthr) : 0.))/(subthr ? 1. : instate.eq.xJ);
}

template <class Eq>
double Interaction<Eq>::width(double xJ, double xlam, double xjf, double xji) {
    std::complex<double> mel = 0.;
    std::complex<double> buf;
    for (int i=1; i <= std::lround((xJ-1.)/2.); ++i) {
        mel += melExJ(2.*i+1., xlam, xjf, xji) + melMxJ(2.*i+1., xlam, xjf, xji);
        /*std::cout  << "ExJ::> "<< 2*i+1 << "," << xlam << "," << xjf << "," << xji << ": " << melExJ(2.*i+1., xlam, xjf, xji) << std::endl;*/
        /*std::cout << "MxJ::> " << 2*i+1 << "," << xlam << "," << xjf << "," << xji << ": " << melMxJ(2.*i+1., xlam, xjf, xji) << std::endl;*/
    }
    return widthMel(mel);
}

template <class Eq>
double Interaction<Eq>::width(double xJ, double xlam, double xj, bool subthr) {
    double res = 0.;

    if (!subthr) {
        for (int jf_cnt = 0; jf_cnt < std::lround(outstate.eq.xJ); ++jf_cnt) {
            res += width(xJ, xlam, -outstate.eq.xJ+2.+2.*(double)jf_cnt, xj);
        }
        return res;
    } else {
        for (int ji_cnt = 0; ji_cnt < std::lround(instate.eq.xJ); ++ji_cnt) {
            res += width(xJ, xlam, xj, -instate.eq.xJ+2.+2.*(double)ji_cnt);
        }
        return res/instate.eq.xJ;
    }

}

template <class Eq>
double Interaction<Eq>::width(double xJ, double xlam) {
    double res = 0.;
    for (int ji_cnt = 0; ji_cnt < std::lround(instate.eq.xJ); ++ji_cnt) {
        res += width(xJ, xlam, -instate.eq.xJ+2.+2.*(double)ji_cnt);
    }
    return res/instate.eq.xJ;
}

template <class Eq>
double Interaction<Eq>::width(double xJ) {
    return width(xJ, -1.) + width(xJ, 3.);
}

template <class Eq>
double Interaction<Eq>::widthHel(double xJ, double xH, bool subthr) {
    return (width(xJ, -1., xH, subthr) + width(xJ, 3., xH, subthr) + (std::lround(xH) != 1 ? width(xJ, -1., -xH+2., subthr) + width(xJ, 3., -xH+2., subthr) : 0.))/(subthr ? 1. : instate.eq.xJ);
}

template <class Eq>
double  Interaction<Eq>::reduceWidth(double width) {
    double Mf = 2*outstate.eq.env.mC + outstate.eq.E;
    double Mi = 2*instate.eq.env.mC + instate.eq.E;
    /*cout << "reduce: " << width << ", " << alphaEM/2.*(Mi*Mi*Mi*Mi - Mf*Mf*Mf*Mf)/(Mi*Mi*Mi) << ", " << (Mi - Mf) << endl;*/
    return 1./(alphaEM/2.*(Mi*Mi*Mi*Mi - Mf*Mf*Mf*Mf)/(Mi*Mi*Mi)/width);
}

