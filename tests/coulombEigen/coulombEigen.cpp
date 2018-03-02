#include <iostream>
#include <fstream>
#include <vector>
#include <boost/format.hpp>
#include "env_deng2016lin.hpp"
#include "eq_coul.hpp"
#include "theeigenval.hpp"

using namespace std;

int main() {

    TheEigenVal<EqCoul<EnvLin> > evals;
   
    evals.eq.xJ = 3;
    evals.eq.xL = 1;
    evals.eq.xS = 3;
    evals.eq.xS1 = 2;
    evals.eq.xS2 = 2;
    evals.eq.env.alphaS = 0.5461;
    evals.eq.env.b = 0.1425;
    evals.eq.env.mC = 1.4830;
    evals.eq.env.muR = evals.eq.env.mC/2;
    evals.eq.env.sigma = 1.1384;
    evals.eq.env.rC = 1E-6;
    evals.intstep = 1E-3;
    evals.stpr = stepper(1E-8, 0);

    vector<double> cutscales{50., 75., 150.};

    double minE = -0.026;
    double maxE = -0.006;
    int steps = 10000;
    double step =  (maxE - minE)/steps;

    for (double cutscale: cutscales) {
        ofstream fout((boost::format("asymp-%1%.dat") % (int)cutscale).str().c_str());
        
        evals.cutscale = cutscale;
        evals.eq.E = minE;
        for (int i = 0; i<steps; i++) {
           fout << evals.eq.E << "," << evals.f() << endl;
           evals.eq.E += step;
        }
        
        fout.close();
    }

    return 0;
}
