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
    *  Exercise 08.1
    ***************************/

    // initialization values
    double x0(0.);
    double mu(0.8);
    double sigma(0.5);
    double metro(2.7);
    bool hist(false);
    ofstream ofs("../out/x.out");

    auto t1_unif = [&rnd, &metro](double cc){
          double ccn(rnd.Rannyu(-metro, metro));

          return cc + ccn;
      };

    auto density = [&mu, &sigma](double x){
          return psi2(x, mu, sigma);
      };

    auto func_blk = [&rnd, &density, &t1_unif, &x0, &mu, &sigma, &ofs, &hist](unsigned int L){
          vector<double> vct(L);

          auto gen = [&rnd, &density, &t1_unif, &x0, &mu, &sigma, &ofs, &hist](){
                x0 = rnd.Metropolis1d(density, t1_unif, x0);
                if (hist) ofs << x0 << endl;

                return intg(x0, mu, sigma);
            };

          generate(vct.begin(), vct.end(), gen);

          return accumulate(vct.begin(), vct.end(), 0.0d) / L;
      };

    // Optimize mu and sigma, to reach best approx of GS energy
    double mu_eps(0.1);
    double sigma_eps(0.1);
    unsigned int nrange(1e2);
    auto mu_v(getRange(nrange, mu - mu_eps, mu + mu_eps));
    auto sigma_v(getRange(nrange, sigma - sigma_eps, sigma + sigma_eps));

    vector<double> sigma_opt;
    vector<double> mu_opt;

    for (auto el : sigma_v) {
        sigma = el;
        sigma_opt.push_back(func_blk(5e5));
    }
    auto sigma_min = distance(sigma_opt.begin(), min_element(sigma_opt.begin(), sigma_opt.end()));
    sigma = sigma_v[sigma_min];

    for (auto el : mu_v) {
        mu = el;
        mu_opt.push_back(func_blk(5e5));
    }
    auto mu_min = distance(mu_opt.begin(), min_element(mu_opt.begin(), mu_opt.end()));
    mu = mu_v[mu_min];

    cout << "mu: " << mu << ", sigma: " << sigma << endl;

    /***************************
    *  Exercise 08.2
    ***************************/

    rnd.resetMetroRatio();
    hist = true;
    blockingMethod(func_blk, 1e5, 100, "E-gs");
    cout << "Metropolis ratio: " << rnd.getMetroRatio() << endl;
    ofs.close();
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
