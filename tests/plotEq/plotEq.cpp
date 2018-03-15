#include <iostream>
#include <fstream>
#include <string>
#include "types.hpp"
#include "env_deng2016scr.hpp"
#include "eq_quark.hpp"

using namespace std;

int main() {
    EqQuark<EnvScr> eq;
   
    eq.xJ = 1;
    eq.xL = 3;
    eq.xS = 3;
    eq.xS1 = 2;
    eq.xS2 = 2;
    eq.env.alphaS = 0.3683;
    eq.env.b = 0.2062;
    eq.env.mC = 4.7572;
    eq.env.muR = eq.env.mC/2;
    eq.env.mu = 0.05611;
    eq.env.sigma = 3.1025;
    eq.env.rC = 0.3052;
    eq.E = 0.35148258;

    ofstream fout("points.dat");

    int numpoints = 1000;
    double minr = 5;
    double maxr = 10;
    double r = minr;
    fldarr x;
    x[0] = 1;
    x[1] = 0;

    fldarr answ;

    fout << "r" << "," << "V"  << endl;
    double step = (maxr - minr)/numpoints;
    for (int i=0; i<numpoints; ++i) {
        eq(x, answ, r);
        fout << r << "," << answ[1] << endl;
        r += step;
    }

    fout.close();
    return 0;
}
