#ifndef OBSERVE_U_HPP
#define OBSERVE_U_HPP

#include "types.hpp"

template <typename Eq, typename Observer>
class ObserveU {
public:
    Eq eq;
    double rMin, rMax, step;

    fldarr u;
    stepper stpr;
    Observer obs;

    ObserveU();
    void operator() ();
};

#include "observe_u.tpp"

#endif
