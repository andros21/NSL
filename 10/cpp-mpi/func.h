#include "../../RC/random.h"

void
randomOnCircle(unsigned int npoints, double radius = 5., string ofile = "circle");

void
randomInSquare(unsigned int npoints, double edge = 5., string ofile = "square");

#ifndef __Permutation__
# define __Permutation__
class Permutation {
private:
    vector<int> _workers;

public:
    Permutation() : _workers({ 0, 1, 2, 3 }){ }

    void
    doPermutation()
    {
        random_shuffle(_workers.begin(), _workers.end());
    };
    pair<int, int>
    get1Pair()
    {
        return pair<int, int>(_workers[0], _workers[1]);
    }

    pair<int, int>
    get2Pair()
    {
        return pair<int, int>(_workers[2], _workers[3]);
    }
};

#endif /* ifndef __Permutation__ */
