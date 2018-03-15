#include <iostream>
#include <fstream>
#include <string>
#include <boost/format.hpp>
#include "env_deng2016lin.hpp"
#include "eq_coul.hpp"
#include "observe_u.hpp"
#include "observers.hpp"

using namespace std;

int main() {

    ObserveU<EqCoul<EnvLin>, uObserverToFile > observu;
   
    observu.eq.xJ = 3;
    observu.eq.xS = 3;
    observu.eq.xS1 = 2;
    observu.eq.xS2 = 2;
    observu.eq.env.alphaS = 0.5461;
    observu.eq.env.b = 0.1425;
    observu.eq.env.mC = 1.4830;
    observu.eq.env.muR = observu.eq.env.mC/2;
    observu.eq.env.sigma = 1.1384;
    observu.eq.env.rC = 1E-2;
    observu.rMin = 0;
    observu.rMax = 200;
    observu.step = 1E-4;
    observu.stpr = stepper();

    double R = 8./9.*observu.eq.env.muR*observu.eq.env.alphaS*observu.eq.env.alphaS;
    int maxN = 5;

    for (int n=1; n<=maxN; n++) {
        for (int l=0; l<n; l++) {
            ofstream fout((boost::format("points_%1%-%2%.dat") % n % l).str().c_str());
            observu.eq.xL = 2.*l+1.;
            observu.eq.E = -R/n/n;
            observu.obs.fout = &fout;

            fout << "r" << "|" << "u" << "|" << "du" << endl;
            observu();

            fout.close();
        }
    }
    

    return 0;
}
