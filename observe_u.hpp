#ifndef OBSERVE_U_HPP
#define OBSERVE_U_HPP

#include <nlohmann/json.hpp>
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

template <typename Eq, typename Observer>
void from_json(const nlohmann::json &j, ObserveU<Eq, Observer>& p);

#include "observe_u.tpp"

#endif
