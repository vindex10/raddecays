#include <cmath>
#include <boost/numeric/odeint.hpp>


template <typename Eq, typename Obs>
void ObserveU<Eq, Obs>::operator() () {
    u[0] = 0;
    u[1] = eq.xL < 1.2 ? 1 : 0; // xL is double, so == is bad idea

    boost::numeric::odeint::integrate_const(stpr, eq, u, rMin, rMax, step, obs);
}

template <typename Eq, typename Obs>
ObserveU<Eq, Obs>::ObserveU () {}
