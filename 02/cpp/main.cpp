/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/

#include "../../RC/random.h"
#include "func.h"

using namespace std;

int
main(int argc, char * argv[])
{
    Random rnd;

    rnd.SetRandom(); // set random default seed

    /***************************
    *  Exercise 02.1
    ***************************/
    // cos function
    auto cosf = [](double x){
          return cos(M_PI * x * 0.5);
      };
    // cos Taylor 2 order function normalized
    double norm = 1. - (1. / 6.) * pow(M_PI * 0.5, 2);
    auto cosft = [&norm](double x){
          return (1. - 0.5 * pow(M_PI * x * 0.5, 2)) * pow(norm, -1);
      };

    // cos / cos Taylor function
    auto cosf_mod = [&cosf, &cosft](double x){
          return cosf(x) / cosft(x);
      };

    // Simulation #1 no important sampling
    auto I = [&rnd, &cosf](unsigned int L){
          vector<double> vct(getRannyuf1d(rnd, L, cosf));

          return M_PI * 0.5 * (accumulate(vct.begin(), vct.end(), 0.0d) / L);
      };

    // Simulation #2 yes important sampling
    auto I2 = [&rnd, &cosf_mod, &cosft, &norm](unsigned int L){
          vector<double> vct(getRannyuf1d(rnd, L, cosf_mod, cosft, pow(norm, -1)));

          return M_PI * 0.5 * (accumulate(vct.begin(), vct.end(), 0.0d) / L);
      };

    blockingMethod(I, 1000000, 100, "I1");
    blockingMethod(I2, 1000000, 100, "I2");

    /***************************
    *  Exercise 02.2
    ***************************/

    return 0;
} // main

/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/
