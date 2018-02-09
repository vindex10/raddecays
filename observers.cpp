#include <fstream>
#include <iostream>
#include "types.hpp"
#include "observers.hpp"

uObserverToFile::uObserverToFile(std::ofstream* fout): fout(fout) {}

void uObserverToFile::operator() (const fldarr &u, double t) {
    *fout << t << "|" << u[0] << "|" <<  u[1] << std::endl;
}
