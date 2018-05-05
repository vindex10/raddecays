#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <boost/filesystem.hpp>
#include "json_types.hpp"
#include "eq_quark.hpp"
#include "state.hpp"
#include "interaction.hpp"
#include "utils.hpp"

using namespace std;
namespace fs = boost::filesystem;

template <class Eq>
void init_state(State<Eq>& state, string statename, string datapath, vector<json*> &cfg) {
    for (json*& item: cfg) {
        state = item->at(statename.c_str());
    }

    ifstream dataf((datapath+statename).c_str());
    state.data = readPlot(dataf);
    dataf.close();

    state.init_spline();
    state.renorm();
    state.init_spline();
}

int main(int argc, char* argv[]) {

    cout.precision(14);

    string cfgname = fs::path(argv[1]).stem().string();
    cout << "cfg loaded: " << cfgname << endl;

    string title = PREFIX+"."+cfgname;
    string outdir = "output/"+title;
    fs::create_directories(outdir);

    json cfgP = gain_config(argv[1], (outdir+"/config").c_str());

    string system = cfgP["system"].get<string>();
    string prefix = system.substr(system.rfind("-")+1);
    
    json eigenP = gain_config(("../quarkEigen/output/"+prefix+"."+system+"/config").c_str(), (outdir+"/eigen_config").c_str());
    json uP = gain_config(("../quarkU/output/"+prefix+"."+system+"/config").c_str(), (outdir+"/u_config").c_str());

    string reportPath = outdir+"/data";
    touchCSV(reportPath.c_str()
            , "instate,outstate,width,R_width,widthE1,R_widthE1,widthELW,R_widthELW");
    ofstream report(reportPath.c_str(), ios_base::app);
    report << fixed << setprecision(16);

    for (json::iterator trans = cfgP["specs"].begin(); trans != cfgP["specs"].end(); ++trans) {
#if defined(ENV_DENG2016LIN_HPP)
        Interaction<EqQuark<EnvLin> > inter;
#elif defined(ENV_DENG2016SCR_HPP)
        Interaction<EqQuark<EnvScr> > inter;
#endif

        inter.alphaEM = cfgP["alphaEM"].get<double>();
        json params = trans.value();

        string datapath = "../quarkU/output/"+prefix+"."+system+"/data/";
        vector<json*> configs{&eigenP, &uP};
        string instateName = params["instate"].get<string>();
        string outstateName = params["outstate"].get<string>();

        init_state(inter.instate, instateName, datapath, configs);
        init_state(inter.outstate, outstateName, datapath, configs);

        cout << " -------------------------------- " << endl;
        cout << instateName << " -> " << outstateName << endl;
        cout << " -------------------------------- " << endl;

        double prev, cur;

        cout << "EJ+MJ width" << endl;
        prev = 0;
        for (int i=1; i <= 3; i++) {
            cur = inter.width(2.*i+1.);
            cout << i << "," << cur << " (" << inter.reduceWidth(cur) << ")," << cur - prev << endl;
            prev = cur;
        }

        cout << "EJ widths" << endl;
        prev = 0;
        for (int i=1; i <= 3; i++) {
            cur = inter.widthExJ(2.*i+1.);
            cout << i << "," << cur << " (" << inter.reduceWidth(cur) << ")," << cur-prev << endl;
            prev = cur;
        }

        cout << "MJ widths" << endl;
        prev = 0;
        for (int i=1; i <= 3; i++) {
        cur = inter.widthMxJ(2.*i+1.);
        cout << i << "," << cur << " (" << inter.reduceWidth(cur) << ")," << cur-prev << endl;
        prev = cur;
        }
        cout << endl;

        cout << "E1 h=0 " << inter.widthExJHel(3., 1.)  << "," << inter.widthExJHel(3., 1.)/inter.widthExJ(3.) << endl;
        cout << "tot1 h=0 " << inter.widthHel(3., 1.)  << "," << inter.widthHel(3., 1.)/inter.width(3.) << endl;
        cout << "E5 h=0 " << inter.widthExJHel(11., 1.)  << "," << inter.widthExJHel(11., 1.)/inter.widthExJ(11.) << endl;
        cout << "tot5 h=0 " << inter.widthHel(11., 1.)  << "," << inter.widthHel(11., 1.)/inter.width(11.) << endl;
        cout << "E1 h=2 " << inter.widthExJHel(3., 5.)  << "," << inter.widthExJHel(3., 5.)/inter.widthExJ(3.) << endl;
        cout << "tot1 h=2 " << inter.widthHel(3., 5.)  << "," << inter.widthHel(3., 5.)/inter.width(3.) << endl;
        cout << "E5 h=2 " << inter.widthExJHel(11., 5.)  << "," << inter.widthExJHel(11., 5.)/inter.widthExJ(11.) << endl;
        cout << "tot5 h=2 " << inter.widthHel(11., 5.)  << "," << inter.widthHel(11., 5.)/inter.width(11.) << endl;

        cout <<endl;

        report      << instateName
             << "," << outstateName
             << "," << inter.width(11.)
             << "," << inter.reduceWidth(inter.width(11.))
             << "," << inter.widthExJ(3.)
             << "," << inter.reduceWidth(inter.widthExJ(3.))
             << "," << inter.widthELW()
             << "," << inter.reduceWidth(inter.widthELW())
             << endl;

    }
    report.close();


    return 0;
}