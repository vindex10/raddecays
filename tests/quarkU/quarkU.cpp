#include <iostream>
#include <fstream>
#include <string>
#include <boost/format.hpp>
#include "env_deng2016lin.hpp"
#include "eq_quark.hpp"
#include "observe_u.hpp"
#include "observers.hpp"

using namespace std;

int main() {

    ObserveU<EqQuark<EnvLin>, uObserverToFile > observu;
   
    observu.eq.xJ = 3;
    observu.eq.xL = 1;
    observu.eq.xS = 3;
    observu.eq.xS1 = 2;
    observu.eq.xS2 = 2;
    observu.eq.env.alphaS = 0.5461;
    observu.eq.env.b = 0.1425;
    observu.eq.env.mC = 1.4830;
    observu.eq.env.muR = observu.eq.env.mC/2;
    observu.eq.env.sigma = 1.1384;
    observu.eq.env.rC = 0.202;
    observu.rMin = 0;
    observu.rMax = 15;
    observu.step = 1E-4;
    observu.stpr = stepper(1E-16, 0);
    observu.eq.E = 0.13002398333501;

    ofstream fout("points.dat");
    observu.obs.fout = &fout;

    fout << "r" << "|" << "u" << "|" << "du" << endl;
    observu();

    fout.close();

    return 0;
}
