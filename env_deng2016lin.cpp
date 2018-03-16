#include <cmath>
#include <iostream>
#include <nlohmann/json.hpp>
#include <gsl/gsl_sf_coupling.h>
#include "env_deng2016lin.hpp"

double EnvLin::smearedDelta (double r) {
    return std::pow(sigma/std::sqrt(M_PI), 3)*std::exp(-sigma*sigma*r*r);
}

double EnvLin::Vv (double r) {
    return -4./3.*alphaS/r;
}

double EnvLin::dVv(double r) {
    return 4./3.*alphaS/r/r;
}

double EnvLin::ddVv(double r) {
    return -8./3.*alphaS/r/r/r;
}

double EnvLin::Vs(double r) {
    return b*r;
}

double EnvLin::dVs(double r) {
    return b;
}

double EnvLin::Vss(double r, double xS, double xS1, double xS2) {
    return 32.*M_PI*alphaS/9./mC/mC*smearedDelta(r)*((xS*xS-1.) - (xS1*xS1-1.) - (xS2*xS2-1.))/8.;
}

double EnvLin::Vsl(double r, double xJ, double xL, double xS)  {
    return 1./2./mC/mC*(3.*dVv(r) - dVs(r))/r*((xJ*xJ-1.) - (xL*xL-1.) - (xS*xS - 1.))/8.;
}

double EnvLin::St(double xJ, double xL, double xS) {
    return 2.*(fmod((xL-1.+xS-1.-(xJ-1.))/2 + 1., 2.) < 0.5 ? 1. : -1.)*std::sqrt(xL*(xL+1.)*(xL-1.)/4./(xL-2.)/(xL+2))*std::sqrt(xS*(xS+1.)*(xS-1.)*(xS-2.)*(xS+2.)/4.)*gsl_sf_coupling_6j(std::round(xL-1.), std::round(xJ-1.), std::round(xS-1.), std::round(xS-1.), 4, std::round(xL-1.));
}

double EnvLin::Vt(double r, double xJ, double xL, double xS) {
    return 1./12./mC/mC*(1./r*dVv(r) - ddVv(r))*St(xJ, xL, xS);
}

void from_json(const nlohmann::json& j, EnvLin& p) {
    try {
        p.alphaS = j.at("alphaS").get<double>();
    } catch(nlohmann::json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(nlohmann::json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }
    
    try {
        p.b = j.at("b").get<double>();
    } catch(nlohmann::json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(nlohmann::json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }
    
    try {
        p.mC = j.at("mC").get<double>();
    } catch(nlohmann::json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(nlohmann::json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }
    
    try {
        p.muR = p.mC/2.;
    } catch(nlohmann::json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(nlohmann::json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }
    
    try {
        p.sigma = j.at("sigma").get<double>();
    } catch(nlohmann::json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(nlohmann::json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }

    try {
        p.rC = j.at("rC").get<double>();
    } catch(nlohmann::json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(nlohmann::json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }
}
