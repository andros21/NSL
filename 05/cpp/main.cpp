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
    *  Exercise 05.1.1
    ***************************/

    double a0(0.0529); // Bohr radius [nm]

    vector<double> pt_100 = { 0., 0., 0. }; // starting point for eign_100
    double var_100(a0 * 1.2);               // parameter for pt_100 and eign_100
    ;                                       // to reach 50% empirical rule

    // 100 eigenfunction take input a cartesian point cc=(x,y,z)
    auto eign_100 = [&a0](vector<double> cc){
          double r = sqrt(accumulate(cc.begin(), cc.end(), 0.0d, [](double a, double b){
            return a + pow(b, 2);
        }));

          return pow(pow(a0, -3. / 2.) / sqrt(M_PI) * exp(-r / a0), 2);
      };

    vector<double> pt_210 = { 0., 0., 0. }; // starting point for eign_210
    double var_210(a0 * 3);                 // parameter for pt_210 and eign_210
    ;                                       // to reach 50% empirical rule

    // 210 eigenfunction take input a cartesian point cc=(x,y,z)
    auto eign_210 = [&a0](vector<double> cc){
          double r = sqrt(accumulate(cc.begin(), cc.end(), 0.0d, [](double a, double b){
            return a + pow(b, 2);
        }));
          double cost = cc[2] / r;

          return pow(pow(a0, -5. / 2.) / 8. * sqrt(2. / M_PI) * r * exp(-r / (2 * a0)) * cost, 2);
      };

    /* 100 */

    // trivial transition matrix random uniform from (-var_100,var_100)
    // take input a cartesian point cc=(x,y,z)
    auto t1_unif_100 = [&rnd, &var_100](vector<double> cc){
          vector<double> ccn(getRannyu(rnd, cc.size(), -var_100, var_100));

          transform(ccn.begin(), ccn.end(), cc.begin(), ccn.begin(), plus<double>() );

          return ccn;
      };

    // <r>_{100} iterator for blocking method
    auto eign_100_blk = [&rnd, &eign_100, &t1_unif_100, &pt_100](unsigned int L){
          vector<vector<double> > vct(L);

          auto gen = [&rnd, &eign_100, &t1_unif_100, &pt_100](){
                // rnd.coutMetroRatio();
                pt_100 = rnd.Metropolis3d(eign_100, t1_unif_100, pt_100);

                return pt_100;
            };

          generate(vct.begin(), vct.end(), gen);

          auto acc = [](double & a, vector<double> & b){
                auto norm = [](double ac, double bc){
                      return ac + pow(bc, 2);
                  };

                return a + sqrt(accumulate(b.begin(), b.end(), 0.0d, norm));
            };

          return accumulate(vct.begin(), vct.end(), 0.0d, acc) / L;
      };


    blockingMethod(eign_100_blk, 5e4, 100, "r-100-eq"); // equilibration <r>_{100}
    blockingMethod(eign_100_blk, 1e5, 100, "r-100");    // simulation    <r>_{100}
    rnd.resetMetroRatio();

    /* 210 */

    // trivial transition matrix random uniform from (-var_210,var_210)
    // take input a cartesian point cc=(x,y,z)
    auto t1_unif_210 = [&rnd, &var_210](vector<double> cc){
          vector<double> ccn(getRannyu(rnd, cc.size(), -var_210, var_210));

          transform(ccn.begin(), ccn.end(), cc.begin(), ccn.begin(), plus<double>() );

          return ccn;
      };

    // <r>_{210} iterator for blocking method
    auto eign_210_blk = [&rnd, &eign_210, &t1_unif_210, &pt_210](unsigned int L){
          vector<vector<double> > vct(L);

          auto gen = [&rnd, &eign_210, &t1_unif_210, &pt_210](){
                // rnd.coutMetroRatio();
                pt_210 = rnd.Metropolis3d(eign_210, t1_unif_210, pt_210);

                return pt_210;
            };

          generate(vct.begin(), vct.end(), gen);

          auto acc = [](double & a, vector<double> & b){
                auto norm = [](double ac, double bc){
                      return ac + pow(bc, 2);
                  };

                return a + sqrt(accumulate(b.begin(), b.end(), 0.0d, norm));
            };

          return accumulate(vct.begin(), vct.end(), 0.0d, acc) / L;
      };


    blockingMethod(eign_210_blk, 5e4, 100, "r-210-eq"); // equilibration <r>_{210}
    blockingMethod(eign_210_blk, 1e5, 100, "r-210");    // simulation    <r>_{210}
    rnd.resetMetroRatio();

    /***************************
    *  Exercise 05.1.2
    ***************************/

    pt_100  = { 0., 0., 0. }; // starting point for eign_100
    var_100 = a0 * 0.77;      // parameter for pt_100 and eign_100
    ;                         // to reach 50% empirical rule

    pt_210  = { 0., 0., 0. }; // starting point for eign_210
    var_210 = a0 * 2;         // parameter for pt_210 and eign_210
    ;                         // to reach 50% empirical rule

    /* 100 */

    // trivial transition matrix gaus normal from with variance var_100
    // take input a cartesian point cc=(x,y,z)
    auto t1_gaus_100 = [&rnd, &var_100](vector<double> cc){
          vector<double> ccn(getGauss(rnd, cc.size(), 0., var_100));

          transform(ccn.begin(), ccn.end(), cc.begin(), ccn.begin(), plus<double>() );

          return ccn;
      };

    // <r>_{100} iterator for blocking method
    auto eign_gaus_100_blk = [&rnd, &eign_100, &t1_gaus_100, &pt_100](unsigned int L){
          vector<vector<double> > vct(L);

          auto gen = [&rnd, &eign_100, &t1_gaus_100, &pt_100](){
                // rnd.coutMetroRatio();
                pt_100 = rnd.Metropolis3d(eign_100, t1_gaus_100, pt_100);

                return pt_100;
            };

          generate(vct.begin(), vct.end(), gen);

          auto acc = [](double & a, vector<double> & b){
                auto norm = [](double ac, double bc){
                      return ac + pow(bc, 2);
                  };

                return a + sqrt(accumulate(b.begin(), b.end(), 0.0d, norm));
            };

          return accumulate(vct.begin(), vct.end(), 0.0d, acc) / L;
      };

    blockingMethod(eign_gaus_100_blk, 5e4, 100, "r-g-100-eq"); // equilibration <r>_{100}
    blockingMethod(eign_gaus_100_blk, 1e5, 100, "r-g-100");    // simulation    <r>_{100}
    rnd.resetMetroRatio();

    /* 210 */

    // trivial transition matrix gaus normal from with variance var_210
    // take input a cartesian point cc=(x,y,z)
    auto t1_gaus_210 = [&rnd, &var_210](vector<double> cc){
          vector<double> ccn(getGauss(rnd, cc.size(), 0., var_210));

          transform(ccn.begin(), ccn.end(), cc.begin(), ccn.begin(), plus<double>() );

          return ccn;
      };

    // <r>_{100} iterator for blocking method
    auto eign_gaus_210_blk = [&rnd, &eign_210, &t1_gaus_210, &pt_210](unsigned int L){
          vector<vector<double> > vct(L);

          auto gen = [&rnd, &eign_210, &t1_gaus_210, &pt_210](){
                // rnd.coutMetroRatio();
                pt_210 = rnd.Metropolis3d(eign_210, t1_gaus_210, pt_210);

                return pt_210;
            };

          generate(vct.begin(), vct.end(), gen);

          auto acc = [](double & a, vector<double> & b){
                auto norm = [](double ac, double bc){
                      return ac + pow(bc, 2);
                  };

                return a + sqrt(accumulate(b.begin(), b.end(), 0.0d, norm));
            };

          return accumulate(vct.begin(), vct.end(), 0.0d, acc) / L;
      };

    blockingMethod(eign_gaus_210_blk, 5e4, 100, "r-g-210-eq"); // equilibration <r>_{210}
    blockingMethod(eign_gaus_210_blk, 1e5, 100, "r-g-210");    // simulation    <r>_{210}
    rnd.resetMetroRatio();

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
