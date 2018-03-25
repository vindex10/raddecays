#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <cmath>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include "eq_quark.hpp"
#include "theeigenval.hpp"
#include "json_types.hpp"

using namespace std;
namespace fs = boost::filesystem;

int main(int argc, char* argv[]) {
    cout << "working in `" << prefix << "` mode" << endl;

    json config;
    ifstream configFile(argv[1]);
    configFile >> config;
    configFile.close();
    
    string cfgname = fs::path(argv[1]).stem().string();
    string title = prefix+"."+cfgname;
    cout << "cfg loaded: " << cfgname << endl;
    fs::create_directories(("output/"+title).c_str());
    fs::copy_file(argv[1], ("output/"+title+"/"+cfgname+".old.cfg").c_str(), fs::copy_option::overwrite_if_exists);

    fstream exclF(("output/"+title+"/exclude").c_str(), fstream::in|fstream::out|fstream::app);
    set<string> toExclude;
    exclF.seekg(0);
    string buf;
    cout <<endl;
    cout << "Exclusions:" << endl;
    while (exclF >> buf) {
        toExclude.insert(buf);
        cout << buf << endl;
    }
    cout << endl;
    exclF.clear();
    exclF.seekp(0, ios_base::end);

    for (json::iterator particle=config.begin(); particle != config.end(); ++particle) {
        if (toExclude.find(particle.key()) != toExclude.end()) {
            cout << "* excluding " << particle.key() << endl;
            continue;
        }
        string outdir = "output/"+title+"/data/"+particle.key();
        fs::create_directories(outdir.c_str());

        json p = particle.value();

#if defined(ENV_DENG2016LIN_HPP)
        TheEigenVal<EqQuark<EnvLin> > evals;
#elif defined(ENV_DENG2016SCR_HPP)
        TheEigenVal<EqQuark<EnvScr> > evals;
#endif

        evals = p;

        vector<double> cutscales = p["cutscales"].get<vector<double> >();
        double minE = p["eq"]["E"].get<double>()-evals.ewindow;
        double maxE = p["eq"]["E"].get<double>()+evals.ewindow;
        int steps = p["steps"].get<int>();
        double step =  (maxE - minE)/steps;
        double fval=0;

        ofstream minEf((outdir+"/minE.dat").c_str());
        cout.precision(p["outprec"].get<int>());
        minEf.precision(p["outprec"].get<int>());
        cout << particle.key() << "(" << p["eq"]["E"] << ")" << endl;
        double maxVal = DBL_MIN;
        for (double cutscale: cutscales) {
            ofstream fout((outdir+"/asymp-"+boost::str(boost::format("%g")%cutscale)+".dat").c_str());
            evals.cutscale = cutscale;
            evals.eq.E = minE;
            for (int i = 0; i<steps; i++) {
                double val = evals.f();
               fout << evals.eq.E << "," << val << endl;
               maxVal = val > maxVal ? val : maxVal;
               evals.eq.E += step;
            }
            fout.close();
            
            evals.eq.E = p["eq"]["E"].get<double>();
            if (argc == 2) {
                fval = evals.findmin();
            }
            minEf << cutscale << "," << evals.eq.E << endl;
            cout << cutscale << "," << evals.eq.E << "," << fval << endl;
        }
        config[particle.key()]["etol"] = pow(10, log10(maxVal) - p["damping"].get<double>());
        exclF << particle.key() << endl;
        minEf.close();
    }
    exclF.close();
    
    ofstream newcfg("output/"+title+"/"+cfgname+".cfg");
    newcfg << config;
    newcfg.close();

    return 0;
}
