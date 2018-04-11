#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <boost/filesystem.hpp>
#include <gsl/gsl_sf_coupling.h>
#include "json_types.hpp"

std::vector<std::vector<double> > readPlot(std::ifstream &f, bool hashead) {
    std::string line;
    if (hashead) {
        std::getline(f, line);
    }

    double val;
    std::string part;
    std::vector<double> x,y;
    while (std::getline(f, line)) {
        std::stringstream str(line);
        std::vector<double> tmp;
        while (std::getline(str, part, ',')) {
            std::stringstream strval(part);
            strval.precision(14);
            strval >> val;
            tmp.push_back(val);
        }
        x.push_back(tmp[0]);
        y.push_back(tmp[1]);
    }
    std::vector<std::vector<double> > out;
    out.push_back(x);
    out.push_back(y);
    return out;
}

#include <iostream>
double clebsch(double xJ1, double xm1, double xJ2, double xm2, double xJ, double xm) {
    return (std::lround((-xJ1+xJ2-xm+1.)/2.)%2 == 0 ? 1. : -1.)*sqrt(xJ)*gsl_sf_coupling_3j(std::lround(xJ1-1.), std::lround(xJ2-1.), std::lround(xJ-1.), std::lround(xm1-1.), std::lround(xm2-1.), std::lround(1.-xm));
}

json gain_config(const char* cfgpath, const char* storeto) {
    boost::filesystem::copy_file(cfgpath, storeto, boost::filesystem::copy_option::overwrite_if_exists);
    std::ifstream cfgF(storeto);
    json cfg;
    cfgF >> cfg;
    cfgF.close();
    return cfg;
}

int touchCSV(const char* path, const char* header) {
    if (!boost::filesystem::exists(path)) {
        std::ofstream f(path);
        f << header << std::endl;
        return 1;
    }
    return 0;
}
