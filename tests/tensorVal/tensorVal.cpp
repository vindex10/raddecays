#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "env_deng2016lin.hpp"

using namespace std;

int main() {

    EnvLin env;
   
    ofstream fout("vals.dat");
    fout << "xJ,xL,xS,ST" << endl;
    
    int minxJ, maxxJ;
    for (int xs = 1; xs <= 3; xs+=2) {
        for (int xl=1; xl <= 11; xl+=2) {
            minxJ = min(abs(xs-xl)+1, abs(xs+xl-2)+1);
            maxxJ = max(abs(xs-xl)+1, abs(xs+xl-2)+1);
            for (int xj=minxJ; xj<=maxxJ; xj+=2) {
                fout << (xj-1.)/2. <<","<< (xl-1.)/2. << "," << (xs-1.)/2. << "," << env.St(xj, xl, xs) << endl;
            }
        }
    }

    fout.close();

    return 0;
}
