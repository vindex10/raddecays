#include <fstream>
#include <iostream>
#include "odeint_types.hpp"
#include "observers.hpp"

void uObserverToFile::operator() (const fldarr &u, double t) {
    *fout << t << "," << u[0] << std::endl;
}
