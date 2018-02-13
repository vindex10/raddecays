#ifndef TYPES_HPP
#define TYPES_HPP

#include <complex>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <Eigen/Dense>

//typedef std::complex<double> fldvar;
typedef double fldvar;
typedef Eigen::Array<fldvar, 2, 1> fldarr;

namespace boost { namespace numeric { namespace odeint {
template<>
struct vector_space_norm_inf<fldarr>
{
    typedef double result_type;
    double operator()( const fldarr &p ) const
    {
        return std::real((p.conjugate()*p).sum());
    }
};
} } }

//typedef runge_kutta4<fldarr
//typedef runge_kutta_dopri5<fldarr
typedef boost::numeric::odeint::bulirsch_stoer_dense_out<fldarr
                                                         , double
                                                         , fldarr
                                                         , double
                                                         , boost::numeric::odeint::vector_space_algebra
                                                         > stepper; 

#endif
