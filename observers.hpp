#ifndef OBSERVERS_HPP
#define OBSERVERS_HPP

#include <fstream>
#include "odeint_types.hpp"

class uObserverToFile {
    public:
        std::ofstream* fout;
        void operator() (const fldarr &u, double t);
};

#endif
