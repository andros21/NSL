#include "func.h"

// calculate chi test for a random uniform [0,1] sequence
double
chi_test(Random & rnd, unsigned int nthrow, unsigned int nblk)
{
    auto w(getRannyu(rnd, nthrow)); // random vector
    vector<double> chi(nblk);       // chi elements
    auto sub  = 1. / nblk;          // subintervals of [0,1]
    auto nsub = nthrow / nblk;      // number of elements for subintervals

    auto i(0);

    for (auto & chi_el : chi) {
        auto eq = [&](double & d){
              return ((sub * i <= d) and (d < sub * (i + 1)));
          };
        chi_el = (double) pow((count_if(w.begin(), w.end(), eq) - nsub), 2) / nsub;
        i     += 1;
    }
    return accumulate(chi.begin(), chi.end(), 0.0d);
}
