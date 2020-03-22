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
    *  Exercise 01.1.1, 01.1.2
    ***************************/
    // auto v(getRannyu(rnd, 100000)); // generate random vector
    // writeVector(v);                 // write to file random vector
    // blockingMethod(v, 0.5, 100);    // write to file blockingMethod arrays

    /****************************
    *  Exercise 01.1.3
    ****************************/
    // vector<double> chis(100);
    // auto gen = [&](){
    //       return chi_test(rnd, 10000, 100);
    //   };
    // generate(chis.begin(), chis.end(), gen);
    // writeVector(chis, "chi");

    /****************************
    *  Exercise 01.2
    ****************************/
    // vector<double> w;                              // appo vector for lambda function
    // int values[]   = { 1, 2, 10, 100 };            // N
    // string dists[] = { "Gauss", "Exp", "Cauchy" }; // type of distro
    // map<string, map<int, vector<double> > > sum;   // dict of sum vector
    //
    // for (const auto & dist : dists) {
    //     for (const auto & value : values) {
    //         auto gen2 = [&](){
    //               if (dist == "Gauss") w = getGauss(rnd, value);
    //               if (dist == "Exp") w = getExp(rnd, value);
    //               if (dist == "Cauchy") w = getCauchy(rnd, value);
    //
    //               return accumulate(w.begin(), w.end(), 0.0d);
    //           };
    //         sum[dist][value] = vector<double>(10000);
    //         generate(sum[dist][value].begin(), sum[dist][value].end(), gen2);
    //     }
    // }
    //
    // for (auto & p1 : sum) {
    //     for (auto & p2 : p1.second) {
    //         writeVector(p2.second, p1.first + "-" + to_string(p2.first));
    //     }
    // }

    /****************************
    *  Exercise 01.3
    ****************************/
    double d = 1.;
    double L = d - 0.1 * d;
    vector<double> pi(100);


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
