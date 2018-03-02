#include "types.hpp"

#include <fstream>
#include <random>
#include <vector>

#include "env_deng2016lin.hpp"
#include "eq_coul.hpp"
#include "theeigenval.hpp"
#include "scaneigens.hpp"

using namespace std;
using namespace boost::numeric::odeint;

int main() {

    // initialize random numbers for generation of points to scan
    seed_seq seed{time(0)};
    default_random_engine e(seed);
    normal_distribution<double> distr(0, 0.01);

    ScanEigens<TheEigenVal<EqCoul<EnvLin> >, vector<double>> scaneigs;
   
    scaneigs.theeval.eq.xJ = 3;
    scaneigs.theeval.eq.xL = 1;
    scaneigs.theeval.eq.xS = 3;
    scaneigs.theeval.eq.xS1 = 2;
    scaneigs.theeval.eq.xS2 = 2;
    scaneigs.theeval.eq.env.alphaS = 0.5461;
    scaneigs.theeval.eq.env.b = 0.1425;
    scaneigs.theeval.eq.env.mC = 1.4830;
    scaneigs.theeval.eq.env.muR = scaneigs.theeval.eq.env.mC/2;
    scaneigs.theeval.eq.env.sigma = 1.1384;
    scaneigs.theeval.eq.env.rC = 1E-3;
    scaneigs.theeval.cutscale = 300;
    scaneigs.theeval.intstep = 1;
    scaneigs.theeval.etol = 1E-8;
    scaneigs.theeval.estep = 1E-15;
    scaneigs.theeval.stpr = stepper(1E-5, 0.);
    //scaneigs.theeval.stpr = make_controlled(1, 0., runge_kutta_dopri5<fldarr
                                                         //, double
                                                         //, fldarr
                                                         //, double
                                                         //, vector_space_algebra
                                                         //>());
    
    ofstream outE("initE.dat");
    for (int i=0; i<1000; i++) {
        scaneigs.points.push_back(-abs(distr(e)));
        scaneigs.theeval.eq.E = scaneigs.points.back();
        outE << scaneigs.points.back();
        outE.flush();
        outE << " " << scaneigs.theeval.f() << endl;
    }
    outE.flush();
    outE.close();

    //scaneigs.lookup();

    //ofstream outSpec("spectrum.dat");
    //for (double E: scaneigs.found) {
        //outSpec << E << endl;
        //outSpec.flush();
    //}
    //outSpec.close();

    return 0;
}
