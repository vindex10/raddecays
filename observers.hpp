#include <fstream>
#include "types.hpp"

class uObserverToFile {
    private:
        std::ofstream* fout;
    public:
        uObserverToFile(std::ofstream* fout);
        void operator() (const fldarr &u, double t);
};
