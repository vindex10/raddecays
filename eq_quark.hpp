#ifndef EQ_QUARK_HPP
#define EQ_QUARK_HPP

#include "types.hpp"

template <class Env>
class EqQuark {
public:
    Env env;
    double xJ, xL, xS, xS1, xS2, E;

    void operator() (const fldarr &u, fldarr &dudr, const double r);
};

#include "eq_quark.tpp"

#endif
