#include <iostream>
#include <fstream>
#include <vector>
#include <boost/format.hpp>
#include "env_deng2016lin.hpp"
#include "eq_quark.hpp"
#include "theeigenval.hpp"

using namespace std;

int main() {

    TheEigenVal<EqQuark<EnvLin> > evals;
   
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
    evals.eq.env.rC = 0.202;
    evals.intstep = 1E-2;
    evals.stpr = stepper(1E-7, 0);
    evals.etol = 1E-14;
    evals.estep = 1;

    vector<double> cutscales{14,16,18,20,22};

    double minE = 0.131-0.008;
    double maxE = 0.131+0.008;
    int steps = 1600;
    double step =  (maxE - minE)/steps;
    double fval;

    ofstream minEf("minE.dat");
    cout.precision(14);
    minEf.precision(14);
    for (double cutscale: cutscales) {
        ofstream fout((boost::format("asymp-%1%.dat") % (int)cutscale).str().c_str());
        evals.cutscale = cutscale;
        evals.eq.E = minE;
        for (int i = 0; i<steps; i++) {
           fout << evals.eq.E << "," << evals.f() << endl;
           evals.eq.E += step;
        }
        fout.close();
        
        evals.eq.E = (minE+maxE)/2. + 3E-3;
        fval = evals.findmin();
        minEf << cutscale << "," << evals.eq.E << endl;
        cout << cutscale << "," << evals.eq.E << "," << fval << endl;
    }
    minEf.close();

    return 0;
}
