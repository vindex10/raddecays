#include <cmath>
#include <iostream>
#include <gsl/gsl_sf_coupling.h>
#include "env_deng2016scr.hpp"
#include "json_types.hpp"

double EnvScr::smearedDelta (double r) {
    return std::pow(sigma/std::sqrt(M_PI), 3)*std::exp(-sigma*sigma*r*r);
}

double EnvScr::Vv (double r) {
    return -4./3.*alphaS/r;
}

double EnvScr::dVv(double r) {
    return 4./3.*alphaS/r/r;
}

double EnvScr::ddVv(double r) {
    return -8./3.*alphaS/r/r/r;
}

double EnvScr::Vs(double r) {
    return b/mu*(1.-std::exp(-mu*r));
}

double EnvScr::dVs(double r) {
    return b*exp(-mu*r);
}

double EnvScr::Vss(double r, double xS, double xS1, double xS2) {
    if (xS1 < 1.2 || xS2 < 1.2) {
        return 0.;
    } else {
        return 32.*M_PI*alphaS/9./mC/mC*smearedDelta(r)*((xS*xS-1.) - (xS1*xS1-1.) - (xS2*xS2-1.))/8.;
    }
}

double EnvScr::Vsl(double r, double xJ, double xL, double xS)  {
    if (xL < 1.2 || xS < 1.2) {
        return 0;
    } else {
        return 1./2./mC/mC*(3.*dVv(r) - dVs(r))/r*((xJ*xJ-1.) - (xL*xL-1.) - (xS*xS - 1.))/8.;
    }
}

double EnvScr::St(double xJ, double xL, double xS) {
    return 2.*(fmod((xL-1.+xS-1.-(xJ-1.))/2 + 1., 2.) < 0.5 ? 1. : -1.)*std::sqrt(xL*(xL+1.)*(xL-1.)/4./(xL-2.)/(xL+2))*std::sqrt(xS*(xS+1.)*(xS-1.)*(xS-2.)*(xS+2.)/4.)*gsl_sf_coupling_6j(std::round(xL-1.), std::round(xJ-1.), std::round(xS-1.), std::round(xS-1.), 4, std::round(xL-1.));
}

double EnvScr::Vt(double r, double xJ, double xL, double xS) {
    if (xL < 1.2 || xS < 1.2) {
        return 0.;
    } else {
        return 1./12./mC/mC*(1./r*dVv(r) - ddVv(r))*St(xJ, xL, xS);
    }
}

void from_json(const json& j, EnvScr& p) {
    try {
        p.alphaS = j.at("alphaS").get<double>();
    } catch(json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }
    
    try {
        p.b = j.at("b").get<double>();
    } catch(json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }
    
    try {
        p.mu = j.at("mu").get<double>();
    } catch(json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }
    
    try {
        p.mC = j.at("mC").get<double>();
    } catch(json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }
    
    p.muR = p.mC/2.;
    
    try {
        p.sigma = j.at("sigma").get<double>();
    } catch(json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }

    try {
        p.rC = j.at("rC").get<double>();
    } catch(json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }
}
