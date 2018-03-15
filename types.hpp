#ifndef TYPES_HPP
#define TYPES_HPP

#include <complex>
#include <cmath>
#include <vector>
#include <algorithm>
#include <boost/numeric/odeint.hpp>
#include <Eigen/Dense>

//typedef std::complex<double> fldvar;
typedef double fldvar;
typedef Eigen::Array<fldvar, 2, 1> fldarr;
typedef fldvar fldder;

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

//see appendix of 1607.04696v3 for details of method
//system here provides value of T = y''/y (in defs of article above)
//stepper should be initialized with initial values of prevdu = T*y
template<class ST, class DT=double> class coulomb_modified_gowell
{
public:

    typedef ST state_type;
    typedef DT deriv_type;
    typedef double value_type;
    typedef double time_type;
    typedef unsigned short order_type;
    typedef boost::numeric::odeint::stepper_tag stepper_category;
    static order_type order( void ) { return 4; }

    std::vector<deriv_type> prevdu{0., 0.};

    template< class System >
    void do_step( System system , state_type &u , time_type t , time_type dt )
    {
        deriv_type TNow;
        system(u, TNow, t+dt);
        u[1] = ((2.*u[0]+5./6.*dt*dt*prevdu[0]) - (u[1]-1./12.*dt*dt*prevdu[1]))/(1. - 1./12.*dt*dt*TNow);
        std::swap(u[0], u[1]);
        prevdu[1] = prevdu[0];
        prevdu[0] = TNow*u[0];
    }
};

typedef coulomb_modified_gowell<fldarr, fldder> stepper;

//typedef boost::numeric::odeint::runge_kutta4<fldarr
//typedef boost::numeric::odeint::runge_kutta_fehlberg78<fldarr
//typedef boost::numeric::odeint::bulirsch_stoer_dense_out<fldarr
                                                         //, double
                                                         //, fldder
                                                         //, double
                                                         //, boost::numeric::odeint::vector_space_algebra
                                                         //> stepper; 
//typedef boost::numeric::odeint::result_of::make_controlled<boost::numeric::odeint::runge_kutta_dopri5<fldarr
                                                         //, double
                                                         //, fldarr
                                                         //, double
                                                         //, boost::numeric::odeint::vector_space_algebra
                                                         //> >::type stepper; 

#endif
