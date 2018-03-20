#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <set>
#include <vector>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <nlohmann/json.hpp>
#include "eq_quark.hpp"
#include "observers.hpp"
#include "observe_u.hpp"

using namespace std;
using json = nlohmann::json;
namespace fs = boost::filesystem;

int main(int argc, char* argv[]) {

    cout.precision(14);
    cout << "working in `" << prefix << "` mode" << endl;

    string cfgname = fs::path(argv[1]).stem().string();
    cout << "cfg loaded: " << cfgname << endl;

    string title = prefix+"."+cfgname;
    fs::create_directories("output/"+title);
    fs::copy_file(argv[1], "output/"+title+"/config", fs::copy_option::overwrite_if_exists);
    
    fstream exclF(("output/"+title+"/"+"exclude").c_str(), fstream::in|fstream::out|fstream::app);
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

    ifstream pparamF(("../quarkEigen/output/"+title+"/config").c_str());
    json pparams;
    pparamF >> pparams;
    pparamF.close();


    for (json::iterator particle=pparams.begin(); particle != pparams.end(); ++particle) {
        if (toExclude.find(particle.key()) != toExclude.end()) {
            cout << "* excluding " << particle.key() << endl;
            continue;
        }
        cout << particle.key() << endl;

#if defined(ENV_DENG2016LIN_HPP)
        ObserveU<EqQuark<EnvLin>, uObserverToFile > observu;
#elif defined(ENV_DENG2016SCR_HPP)
        ObserveU<EqQuark<EnvScr>, uObserverToFile > observu;
#endif

        observu.eq = particle.value()["eq"];
        observu.rMax = particle.value()["cutscales"].get<vector<double>>().back();
        observu.step = particle.value()["intstep"].get<double>();

        //read fitted energy
        ifstream energF(("../quarkEigen/output/"+title+"/data/"+particle.key()+"/minE.dat").c_str());
        string item,tmp;
        while (getline(energF, tmp)) {
            item = tmp;
        }
        stringstream linestr(item);
        while (getline(linestr, item, ','));
        stringstream i;
        i.precision(14);
        i << item;
        i >> observu.eq.E;

        observu.stpr = stepper();

        fs::create_directory("output/"+title+"/data");
        string outF = "output/"+title+"/data/"+particle.key();
        ofstream fout(outF.c_str());
        observu.obs.fout = &fout;

        fout.precision(14);
        fout << "r" << "," << "u" << endl;
        observu();
        fout.close();

        exclF << particle.key() << endl;
    }
    exclF.close();
   

    return 0;
}
