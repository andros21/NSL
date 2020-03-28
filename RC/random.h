/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <map>
#include <functional>

using namespace std;


#ifndef __Random__
# define __Random__

class Random {
private:
    int m1, m2, m3, m4, l1, l2, l3, l4, n1, n2, n3, n4;

protected:

public:
    Random();
    ~Random();
    void
    SetRandom(int *, int, int);
    void
    SetRandom(string iseed = "../../RC/seed.in");
    void
    SaveSeed();
    double
    Rannyu(void);
    double
    Rannyu(double min, double max);
    double
    Gauss(double mean, double sigma);
    double
    Exp(double lambda);
    double
    Cauchy(double mu, double lambda);
    double
    Sine();
    double
    Cosine();
};

vector<double>
getRannyu(Random & rnd, unsigned int size);
vector<double>
getRannyu(Random & rnd, unsigned int size, double min, double max);
vector<double>
getGauss(Random & rnd, unsigned int size, double mean = 0., double sigma = 1.);
vector<double>
getExp(Random & rnd, unsigned int size, double lambda = 1.);
vector<double>
getCauchy(Random & rnd, unsigned int size, double mu = 0., double lambda = 1.);
vector<double>
getCosine(Random & rnd, unsigned int size);
vector<double>
getSine(Random & rnd, unsigned int size);
void
writeVector(vector<double> & vct, string ofile = "rand");
void
blockingMethod_old(vector<double> & vct, double mean_teor, unsigned int nblk, string oname = "");
void
blockingMethod(function<double(unsigned int)> lbf, unsigned int nthrow, unsigned int nblk,
  string oname = "");

#endif // ifndef __Random__

/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/
