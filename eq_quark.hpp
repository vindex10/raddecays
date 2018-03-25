#ifndef EQ_QUARK_HPP
#define EQ_QUARK_HPP

#include "odeint_types.hpp"
#include "json_types.hpp"

template <class Env>
class EqQuark {
public:
    Env env;
    double xJ, xL, xS, xS1, xS2, E;

    void operator() (const fldarr &u, fldder &dudr, const double r);
    void initU (fldarr &u, double h); //define u at 0
    template<class cnt>
    void initTu (cnt &Tu, double h);
};

template <typename Env>
void from_json(const json &j, EqQuark<Env>& p);

#include "eq_quark.tpp"

#endif
