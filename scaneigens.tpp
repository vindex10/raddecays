#include <iostream>
#include <fstream>

template <typename Eval, typename PointsIter>
void ScanEigens<Eval, PointsIter>::lookup() {
    std::ofstream fout("found.dat");
    for (double point: points) {
        theeval.eq.E = point;
        theeval.findmin();
        if (theeval.eq.E < 0) {
            found.push_back(theeval.eq.E);
            std::cout << theeval.eq.E << std::endl;
            fout << theeval.eq.E << std::endl;
            fout.flush();
            std::cin.get();
        }
        std::cout << "-" << std::endl;
    }
    fout.close();
}
