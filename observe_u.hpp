#ifndef OBSERVE_U_HPP
#define OBSERVE_U_HPP

#include "json_types.hpp"
#include "odeint_types.hpp"

template <typename Eq, typename Observer>
class ObserveU {
public:
    Eq eq;
    double rMax, step;

    fldarr u;
    stepper stpr;
    Observer obs;

    ObserveU();
    void operator() ();
};

template <typename Eq, typename Observer>
void from_json(const json &j, ObserveU<Eq, Observer>& p);

#include "observe_u.tpp"

#endif
