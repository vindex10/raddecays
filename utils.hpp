#ifndef UTILS_HPP
#define UTILS_HPP

#include <fstream>
#include <vector>
#include "json_types.hpp"

std::vector<std::vector<double> > readPlot(std::ifstream &f, bool hashead=false);

double clebsch(double xJ1, double xm1, double xJ2, double xm2, double xJ, double xm);
json gain_config(const char* cfgpath, const char* storeto);
int touchCSV(const char* path, const char* header);

#endif
