#ifndef EQ_COUL_HPP
#define EQ_COUL_HPP

#include "types.hpp"

template <class Env>
class EqCoul {
public:
    Env env;
    double xJ, xL, xS, xS1, xS2, E;

    void operator() (const fldarr &u, fldarr &dudr, const double r);
};

#include "eq_coul.tpp"

#endif
