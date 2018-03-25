#include <iostream>
#include <fstream>
#include <vector>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include "json_types.hpp"
#include "env_deng2016lin.hpp"
#include "eq_coul.hpp"
#include "theeigenval.hpp"

using namespace std;
namespace fs = boost::filesystem;

int main() {

    TheEigenVal<EqCoul<EnvLin> > evals;

    ifstream cfgF("config.cfg");
    json cfg;
    cfgF >> cfg;
    cfgF.close();
   
    evals = cfg;
    evals.stpr = stepper();

    vector<double> cutscales = cfg["cutscales"].get<vector<double> >();

    double minE = cfg["minE"].get<double>();
    double maxE = cfg["maxE"].get<double>();
    int steps = cfg["steps"].get<int>();
    double step =  (maxE - minE)/steps;

    fs::create_directory("output");
    for (double cutscale: cutscales) {
        ofstream fout((boost::format("output/asymp-%g.dat") % (int)cutscale).str().c_str());
        
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
