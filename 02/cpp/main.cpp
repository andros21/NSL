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
    *  Exercise 02.2.1
    ***************************/
    int NN(1000);                    // Number of RWs
    int steps(100);                  // Number of steps
    vector<vector<double> > DRW(NN); // Discrete RW
    vector<double> Dpt(3);           // Discrete RW starting point
    int trial, dir;

    vector<double> DRW_mean(steps);
    vector<double> DRW_dev(steps);

    // square sum
    auto acc = [](double a, double b){
          return a + pow(b, 2);
      };

    // step forward if 0 step reward if 1
    auto fwrw = [] (int t){
          if (t == 0) {
              return -1;
          } else {
              return 1;
          }
      };

    auto genD2 = [&rnd, &Dpt, &acc, &fwrw, &trial, &dir](){
          int trial = rnd.RannyuDiscrete(0, 1);
          int dir   = rnd.RannyuDiscrete(0, 2);

          Dpt[dir] += fwrw(trial);
          return accumulate(Dpt.begin(), Dpt.end(), 0, acc);
      };

    auto genD1 = [&genD2, &steps, &Dpt](){
          vector<double> vct(steps);

          Dpt    = { 0, 0, 0 };
          vct[0] = 0;
          generate(vct.begin() + 1, vct.end(), genD2);
          return vct;
      };

    // Simulate discrete RW
    generate(DRW.begin(), DRW.end(), genD1);

    // Calculate mean and dev of discrete RW for each step
    auto sum(0.), sum2(0.);
    for (int i = 0; i < steps + 1; ++i) {
        sum  = 0;
        sum2 = 0.;
        for (int j = 0; j < NN; ++j) {
            sum  += DRW[j][i];
            sum2 += pow(DRW[j][i], 2);
        }
        DRW_mean[i] = sqrt(sum / NN);
        DRW_dev[i]  = sqrt(sqrt((sum2 / NN - pow(DRW_mean[i], 2)) / NN));
    }

    writeVector(DRW_mean, "DRW_mean"); // write DRW mean
    writeVector(DRW_dev, "DRW_dev");   // write DRW dev

    /***************************
    *  Exercise 02.2.2
    ***************************/
    vector<vector<double> > CRW(NN); // Continuum RW
    vector<double> Cpt(3);           // Continuum RW starting point
    double theta, phi;               // theta phi angle

    vector<double> CRW_mean(steps);
    vector<double> CRW_dev(steps);

    auto genC2 = [&rnd, &Cpt, &acc, &theta, &phi](){
          theta   = rnd.Rannyu(0., M_PI);
          phi     = rnd.Rannyu(0., 2 * M_PI);
          Cpt[0] += sin(theta) * cos(phi);
          Cpt[1] += sin(theta) * sin(phi);
          Cpt[2] += cos(theta);

          return accumulate(Cpt.begin(), Cpt.end(), 0.0d, acc);
      };

    auto genC1 = [&genC2, &steps, &Cpt](){
          vector<double> vct(steps);

          Cpt    = { 0., 0., 0. };
          vct[0] = 0;
          generate(vct.begin() + 1, vct.end(), genC2);
          return vct;
      };

    // Simulate continuum RW
    generate(CRW.begin(), CRW.end(), genC1);

    // Calculate mean and dev of continuum RW for each step
    for (int i = 0; i < steps + 1; ++i) {
        sum  = 0;
        sum2 = 0.;
        for (int j = 0; j < NN; ++j) {
            sum  += CRW[j][i];
            sum2 += pow(CRW[j][i], 2);
        }
        CRW_mean[i] = sqrt(sum / NN);
        CRW_dev[i]  = sqrt(sqrt((sum2 / NN - pow(CRW_mean[i], 2)) / NN));
    }

    writeVector(CRW_mean, "CRW_mean"); // Write CRW mean
    writeVector(CRW_dev, "CRW_dev");   // Write CRW dev

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
