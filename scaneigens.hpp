#ifndef SCANEIGENS_HPP
#define SCANEIGENS_HPP

#include <vector>
#include "types.hpp"

template <typename Eval, typename PointsIter>
class ScanEigens {
    public:
        Eval theeval;
        PointsIter points;

        std::vector<double> found;

        void lookup();
};

#include "scaneigens.tpp"

#endif
