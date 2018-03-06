#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include <boost/format.hpp>
#include <nlohmann/json.hpp>
#include "env_deng2016lin.hpp"
#include "eq_quark.hpp"
#include "theeigenval.hpp"

using namespace std;
using json = nlohmann::json;

int main() {
    
    json config;
    ifstream configFile("particles.cfg");
    configFile >> config;
    configFile.close();

    system("mkdir -p output");
    for (json::iterator particle=config.begin(); particle != config.end(); ++particle) {
        string title = particle.key();
        system(("mkdir -p output/" + title).c_str());

        json p = particle.value();
        TheEigenVal<EqQuark<EnvLin> > evals;
       
        evals.eq.xJ = p["eq"]["xJ"].get<double>();
        evals.eq.xL = p["eq"]["xL"].get<double>();
        evals.eq.xS = p["eq"]["xS"].get<double>();
        evals.eq.xS1 = p["eq"]["xS1"].get<double>();
        evals.eq.xS2 = p["eq"]["xS2"].get<double>();
        evals.eq.env.alphaS = p["eq"]["env"]["alphaS"].get<double>();
        evals.eq.env.b = p["eq"]["env"]["b"].get<double>();
        evals.eq.env.mC = p["eq"]["env"]["mC"].get<double>();
        evals.eq.env.muR = evals.eq.env.mC/2;
        evals.eq.env.sigma = p["eq"]["env"]["sigma"].get<double>();
        evals.eq.env.rC = p["eq"]["env"]["rC"].get<double>();
        evals.intstep = p["intstep"].get<double>();
        evals.stpr = stepper(p["intabsTol"].get<double>(), p["intrelTol"].get<double>());
        evals.etol = p["etol"].get<double>();
        evals.estep = p["estep"].get<double>();
        evals.ewindow = p["ewindow"].get<double>();

        vector<double> cutscales = p["cutscales"].get<vector<double> >();
        double minE = p["eq"]["E"].get<double>()-evals.ewindow;
        double maxE = p["eq"]["E"].get<double>()+evals.ewindow;
        int steps = p["steps"].get<int>();
        double step =  (maxE - minE)/steps;
        double fval;

        ofstream minEf(("output/" + title + "/minE.dat").c_str());
        cout.precision(p["outprec"].get<int>());
        minEf.precision(p["outprec"].get<int>());
        cout << p["eq"]["E"] << endl;
        for (double cutscale: cutscales) {
            ofstream fout(("output/"+title+"/asymp-"+to_string((int)cutscale)+".dat").c_str());
            evals.cutscale = cutscale;
            evals.eq.E = minE;
            for (int i = 0; i<steps; i++) {
               fout << evals.eq.E << "," << evals.f() << endl;
               evals.eq.E += step;
            }
            fout.close();
            
            evals.eq.E = p["eq"]["E"].get<double>();
            fval = evals.findmin();
            minEf << cutscale << "," << evals.eq.E << endl;
            cout << cutscale << "," << evals.eq.E << "," << fval << endl;
        }
        minEf.close();
    }

    return 0;
}