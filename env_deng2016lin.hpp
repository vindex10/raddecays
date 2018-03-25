#ifndef ENV_DENG2016LIN_HPP
#define ENV_DENG2016LIN_HPP

#include "json_types.hpp"

class EnvLin {
public:
    // Free params
    double sigma, b, alphaS, muR, mC, rC;
    
    // Helping functions
    double smearedDelta (double r);

    // Potentials
    double Vv (double r); // Coulomb
    double dVv(double r); // Derivative of Vv
    double ddVv(double r); // 2-Derivative of Vv
    double Vs(double r); // Cornell linear confinement
    double dVs(double r); // Derivative of Vlin
    double Vss(double r, double xS, double xS1, double xS2); // Spin-spin interaction
    double Vsl(double r, double xJ, double xL, double xS); // Spin-orbit interaction
    double St(double xJ, double xL, double xS); // Tensor interaction
    double Vt(double r, double xJ, double xL, double xS); // Tensor interaction
};

void from_json(const json& j, EnvLin& p);

#endif
