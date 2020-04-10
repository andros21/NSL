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
    *  Exercise 03.1.1
    ***************************/
    // map<string, double> par;
    // par["S0"]    = 100;  // Asset price at t=0
    // par["T"]     = 1;    // Delivery date
    // par["K"]     = 100;  // Strike price
    // par["r"]     = 0.1;  // Risk free interest rate
    // par["sigma"] = 0.25; // Volatility
    //
    // // S(t)
    // auto asset_price = [&rnd, &par] (double time){
    //       return par["S0"] * exp((par["r"] - 0.5 * pow(par["sigma"], 2)) + par["sigma"] * rnd.Gauss(0., time));
    //   };
    //
    // // Call numeric function
    // auto call = [&asset_price, &par] (){
    //       return exp(-par["r"] * par["T"]) * max(0., asset_price(par["T"]) - par["K"]);
    //   };
    //
    // // Simulation Call per blk
    // auto call_blk = [&call] (unsigned int L){
    //       vector<double> vct(L);
    //
    //       generate(vct.begin(), vct.end(), call);
    //
    //       return accumulate(vct.begin(), vct.end(), 0.0d) / L;
    //   };
    //
    // // Put numeric function
    // auto put = [&asset_price, &par] (){
    //       return exp(-par["r"] * par["T"]) * max(0., par["K"] - asset_price(par["T"]));
    //   };
    //
    // // Simulation Put per blk
    // auto put_blk = [&put] (unsigned int L){
    //       vector<double> vct(L);
    //
    //       generate(vct.begin(), vct.end(), put);
    //
    //       return accumulate(vct.begin(), vct.end(), 0.0d) / L;
    //   };
    //
    // blockingMethod(call_blk, 100000, 100, "call");
    // blockingMethod(put_blk, 100000, 100, "put");

    /***************************
    *  Exercise 03.1.2
    ***************************/

    // // S(t) discrete
    // auto asset_price_dsc = [&rnd, &par] (double time){
    //       double ST;
    //       double S_appo(par["S0"]);
    //       double dt = time / 100;
    //
    //       for (auto i = 0; i < 100; ++i) {
    //           ST = S_appo * exp((par["r"] - 0.5 * pow(par["sigma"],
    //               2)) * (dt * (i + 1) - dt * i) + par["sigma"] * rnd.Gauss() * sqrt(
    //                   (dt * (i + 1) - dt * i)));
    //           S_appo = ST;
    //       }
    //       return ST;
    //   };
    //
    // // Call numeric function
    // auto call_dsc = [&asset_price_dsc, &par] (){
    //       return exp(-par["r"] * par["T"]) * max(0., asset_price_dsc(par["T"]) - par["K"]);
    //   };
    //
    // // Simulation Call per blk
    // auto call_dsc_blk = [&call_dsc] (unsigned int L){
    //       vector<double> vct(L);
    //
    //       generate(vct.begin(), vct.end(), call_dsc);
    //
    //       return accumulate(vct.begin(), vct.end(), 0.0d) / L;
    //   };
    //
    // // Put numeric function
    // auto put_dsc = [&asset_price_dsc, &par] (){
    //       return exp(-par["r"] * par["T"]) * max(0., par["K"] - asset_price_dsc(par["T"]));
    //   };
    //
    // // Simulation Put per blk
    // auto put_dsc_blk = [&put_dsc] (unsigned int L){
    //       vector<double> vct(L);
    //
    //       generate(vct.begin(), vct.end(), put_dsc);
    //
    //       return accumulate(vct.begin(), vct.end(), 0.0d) / L;
    //   };
    //
    // blockingMethod(call_dsc_blk, 100000, 100, "call-dsc");
    // blockingMethod(put_dsc_blk, 100000, 100, "put-dsc");

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
