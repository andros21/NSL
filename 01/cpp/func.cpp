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
        auto eq = [&i, &sub](double & d){
              return ((sub * i <= d) and (d < sub * (i + 1)));
          };
        chi_el = (double) pow((count_if(w.begin(), w.end(), eq) - nsub), 2) / nsub;
        i     += 1;
    }
    return accumulate(chi.begin(), chi.end(), 0.0d);
}

// Simulate Buffon expirement, d = _width, L = _needle, N_tr = _nthrow
double
Buffon_exp(Random & rnd, double width, double needle, int nthrow)
{
    vector<double> center(getRannyu(rnd, nthrow, 0., width));
    vector<double> cosine(getCosine(rnd, nthrow));

    transform(cosine.begin(), cosine.end(), cosine.begin(), [&needle](double & c){
        return c * needle * 0.5;
    });

    auto right = center;
    transform(right.begin(), right.end(), cosine.begin(), right.begin(), std::plus<double>());
    auto left = center;
    transform(left.begin(), left.end(), cosine.begin(), left.begin(), std::minus<double>());

    auto test = [&width](double & d){
          return d >= width or d <= 0.;
      };

    auto N_dx = count_if(right.begin(), right.end(), test);
    auto N_sx = count_if(left.begin(), left.end(), test);

    return (2. * needle * nthrow) / ( width * (N_sx + N_dx));
} // Buffon_exp
