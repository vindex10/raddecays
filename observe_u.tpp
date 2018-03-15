#include <cmath>
#include <boost/numeric/odeint.hpp>


template <typename Eq, typename Obs>
void ObserveU<Eq, Obs>::operator() () {
    eq.initU(u, step);
    eq.initTu(stpr.prevdu, step);
    boost::numeric::odeint::integrate_adaptive(stpr, eq, u, rMin, rMax, step, obs);
}

template <typename Eq, typename Obs>
ObserveU<Eq, Obs>::ObserveU () {}
