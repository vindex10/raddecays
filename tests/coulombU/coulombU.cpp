#include <iostream>
#include <fstream>
#include <string>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <nlohmann/json.hpp>
#include "env_deng2016lin.hpp"
#include "eq_coul.hpp"
#include "observe_u.hpp"
#include "observers.hpp"

using namespace std;
using json = nlohmann::json;
namespace fs = boost::filesystem;

int main() {

    ObserveU<EqCoul<EnvLin>, uObserverToFile > observu;

    ifstream cfgF("config.cfg");
    json cfg;
    cfgF >> cfg;
    cfgF.close();
   
    observu = cfg;
    observu.stpr = stepper();

    double R = 8./9.*observu.eq.env.muR*observu.eq.env.alphaS*observu.eq.env.alphaS;
    int maxN = cfg["maxN"].get<int>();

    fs::create_directory("output");
    for (int n=1; n<=maxN; n++) {
        for (int l=0; l<n; l++) {
            ofstream fout((boost::format("output/points_%1%-%2%.dat") % n % l).str().c_str());
            observu.eq.xL = 2.*l+1.;
            observu.eq.E = -R/n/n;
            observu.obs.fout = &fout;

            fout << "r" << "," << "u" << "," << "du" << endl;
            observu();

            fout.close();
        }
    }
    

    return 0;
}
