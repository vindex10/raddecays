#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <set>
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
    
    // update eigenP with right energies
    for (json::iterator pcl = eigenP.begin(); pcl != eigenP.end(); ++pcl) {
        string pclname = pcl.key();
        ifstream Edat(("../quarkEigen/output/"+prefix+"."+system+"/data/"+pclname+"/minE.dat").c_str());
        eigenP[pclname]["eq"]["E"] = readPlot(Edat, false)[1].back();
    }

    json uP = gain_config(("../quarkU/output/"+prefix+"."+system+"/config").c_str(), (outdir+"/u_config").c_str());

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


    string reportPath;

    reportPath = outdir+"/widths";
    touchCSV(reportPath.c_str()
            , "instate,outstate,width,in0,in2,out0,out2");
    ofstream report(reportPath.c_str(), ios_base::app);
    report << fixed << setprecision(16);

    reportPath = outdir+"/widthsE1";
    touchCSV(reportPath.c_str()
            , "instate,outstate,width,in0,in2,out0,out2");
    ofstream reportE1(reportPath.c_str(), ios_base::app);
    
    reportPath = outdir+"/widthsELW";
    touchCSV(reportPath.c_str()
            , "instate,outstate,width,in0,in2,out0,out2");
    ofstream reportELW(reportPath.c_str(), ios_base::app);
    report << fixed << setprecision(16);
    
    reportPath = outdir+"/melsq";
    touchCSV(reportPath.c_str()
            , "instate,outstate,melsq,melsqE1,melsqELW");
    ofstream reportTech(reportPath.c_str(), ios_base::app);
    report << fixed << setprecision(16);

    for (json::iterator trans = cfgP["specs"].begin(); trans != cfgP["specs"].end(); ++trans) {
        json params = trans.value();
        string instateName = params["instate"].get<string>();
        string outstateName = params["outstate"].get<string>();
        string transLabel = instateName+"->"+outstateName;

        if (toExclude.find(transLabel) != toExclude.end()) {
            cout << "* excluding " << transLabel << endl;
            continue;
        }

        cout << transLabel << endl;

#if defined(ENV_DENG2016LIN_HPP)
        Interaction<EqQuark<EnvLin> > inter;
#elif defined(ENV_DENG2016SCR_HPP)
        Interaction<EqQuark<EnvScr> > inter;
#endif

        inter.alphaEM = cfgP["alphaEM"].get<double>();

        string datapath = "../quarkU/output/"+prefix+"."+system+"/data/";
        vector<json*> configs{&eigenP, &uP};

        init_state(inter.instate, instateName, datapath, configs);
        init_state(inter.outstate, outstateName, datapath, configs);

        double prev, cur;
        double xJmax = 25.;

        report      << instateName
             << "," << outstateName
             << "," << inter.width(xJmax)
             << "," << inter.widthHel(xJmax,1.)
             << "," << inter.widthHel(xJmax,5.)
             << "," << inter.widthHel(xJmax,1., true)
             << "," << inter.widthHel(xJmax,5., true)
             << endl;

        reportE1    << instateName
             << "," << outstateName
             << "," << inter.widthExJ(3.)
             << "," << inter.widthExJHel(3., 1.)
             << "," << inter.widthExJHel(3., 5.)
             << "," << inter.widthExJHel(3., 1., true)
             << "," << inter.widthExJHel(3., 5., true)
             << endl;

        reportELW   << instateName
             << "," << outstateName
             << "," << inter.widthELW()
             << "," << inter.widthELWHel(1.)
             << "," << inter.widthELWHel(5.)
             << "," << inter.widthELWHel(1., true)
             << "," << inter.widthELWHel(5., true)
             << endl;

        reportTech  << instateName
             << "," << outstateName
             << "," << inter.reduceWidth(inter.width(xJmax))
             << "," << inter.reduceWidth(inter.widthExJ(3.))
             << "," << inter.reduceWidth(inter.widthELW())
             << endl;


        exclF << transLabel << endl;
    }
    report.close();
    reportE1.close();
    reportELW.close();
    reportTech.close();


    return 0;
}
