#include <cmath>
#include "types.hpp"

template <typename Env>
void EqQuark<Env>::operator() (const fldarr &u, fldarr &dudr, const double r) {
    double ruse = r < env.rC ? env.rC : r;
    dudr[0] = u[1];
    dudr[1] = -2*env.muR*(E - env.Vv(ruse) - env.Vs(ruse)
            - (xL*xL - 1.)/4/2/env.muR/ruse/ruse
            - env.Vss(ruse, xS, xS1, xS2)
            - env.Vsl(ruse, xJ, xL, xS)
            - env.Vt(ruse, xJ, xL, xS)
            )*(r < env.rC ? std::pow(env.rC, (xL+1.)/2) : u[0]);
}
