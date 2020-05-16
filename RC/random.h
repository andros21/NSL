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
#include <set>
#include <functional>
#include <iomanip>

using namespace std;


#ifndef __Random__
# define __Random__

class Random {
private:
    int m1, m2, m3, m4, l1, l2, l3, l4, n1, n2, n3, n4;
    unsigned int Macpt, Mttot;

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
    Gauss(double mean = 0., double sigma = 1.);
    double
    Exp(double lambda);
    double
    Cauchy(double mu, double lambda);
    double
    Sine();
    double
    Cosine();
    double
    HitMiss1d(function<double(double)> pdf, double pdfM, double a = 1., double b = 1.);
    double
    Rannyuf1d(function<double(double)> f1, double a = 0., double b = 1.);
    double
    Rannyuf1d(function<double(double)> f1, function<double(double)> pdf, double pdfM, double a = 0.,
      double b = 1.);
    int
    RannyuDiscrete(int a = 0, int b = 0);
    double
    Metropolis1d(function<double(double)> pdf, function<double(double)> t1, double pt = 0.);
    vector<double>
    Metropolis3d(function<double(vector<double>)> pdf, function<vector<double>(vector<double>)> t1,
      vector<double> pt = { 0., 0., 0. });
    void
    resetMetroRatio();
    double
    getMetroRatio();
    void
    coutMetroRatio();
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
vector<double>
getHitMiss1d(Random & rnd, unsigned int size, function<double(double)> pdf, double pdfM, double a = 0., double b = 1.);
vector<double>
getRannyuf1d(Random & rnd, unsigned int size, function<double(double)> f1, double a = 0., double b = 1.);
vector<double>
getRannyuf1d(Random & rnd, unsigned int size, function<double(double)> f1, function<double(double)> pdf, double pdfM,
  double a = 0., double b = 1.);
vector<int>
getRannyuDiscrete(Random &rnd, unsigned int size, int a = 0, int b = 1);
vector<vector<double> >
getMetropolis3d(Random &rnd, unsigned int size, function<double(vector<double>)> pdf,
  function<vector<double>(vector<double>)> t1, vector<double> pt = { 0., 0., 0. });
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
