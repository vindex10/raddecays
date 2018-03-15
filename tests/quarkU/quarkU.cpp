#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/format.hpp>
#include "env_deng2016scr.hpp"
#include "eq_quark.hpp"
#include "observe_u.hpp"
#include "observers.hpp"

using namespace std;

int main() {

    ObserveU<EqQuark<EnvScr>, uObserverToFile > observu;
   
    observu.eq.xJ = 3;
    observu.eq.xL = 1;
    observu.eq.xS = 3;
    observu.eq.xS1 = 2;
    observu.eq.xS2 = 2;
    observu.eq.env.alphaS = 0.3683;
    observu.eq.env.b = 0.2062;
    observu.eq.env.mC = 4.7572;
    observu.eq.env.mu = 0.05611;
    observu.eq.env.muR = observu.eq.env.mC/2;
    observu.eq.env.sigma = 3.1025;
    observu.eq.env.rC = 0.30467198680443935261;
    //observu.eq.env.alphaS = 0.5461;
    //observu.eq.env.b = 0.1425;
    //observu.eq.env.mC = 1.4830;
    //observu.eq.env.muR = observu.eq.env.mC/2;
    //observu.eq.env.sigma = 1.1384;
    //observu.eq.env.rC = 0.202;
    observu.rMin = 0;
    observu.rMax = 12;
    observu.step = 1E-2;
    observu.stpr = stepper();

    //observu.eq.E = 0.35148258; //chi_b0_1P
    //observu.eq.E = 1.09124138; // yps_4S
    observu.eq.E = 0.83844501; //yps_3S



    ofstream fout("points.dat");
    observu.obs.fout = &fout;

    fout.precision(14);
    fout << "r" << "|" << "u" << "|" << "du" << endl;

    observu();

    fout.close();

    return 0;
}
