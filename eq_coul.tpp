#include <cmath>
#include "types.hpp"

template <typename Env>
void EqCoul<Env>::operator() (const fldarr &u, fldarr &dudr, const double r) {
    double ruse = r < env.rC ? env.rC : r;
    dudr[0] = u[1];
    dudr[1] = -2*env.muR*(E - env.Vv(ruse) - (xL*xL - 1.)/4/2/env.muR/ruse/ruse)*(r < env.rC ? std::pow(env.rC, (xL+1.)/2) : u[0]);
}
