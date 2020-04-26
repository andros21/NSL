/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/

#include "random.h"

using namespace std;

static string outPath = "../out/";

Random :: Random(){ }

Random :: ~Random(){ }

void
Random :: SaveSeed()
{
    ofstream WriteSeed;

    WriteSeed.open(outPath + "seed.out");
    if (WriteSeed.is_open()) {
        WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;
    } else { cerr << "ERROR: Unable to open random.out" << endl; }
    WriteSeed.close();
}

double
Random :: Gauss(double mean, double sigma)
{
    double s = this->Rannyu();
    double t = this->Rannyu();
    double x = sqrt(-2. * log(1. - s)) * cos(2. * M_PI * t);

    return mean + x * sigma;
}

double
Random :: Rannyu(double min, double max)
{
    return min + (max - min) * this->Rannyu();
}

double
Random :: Rannyu(void)
{
    const double twom12 = 0.000244140625;
    int i1, i2, i3, i4;
    double r;

    i1 = l1 * m4 + l2 * m3 + l3 * m2 + l4 * m1 + n1;
    i2 = l2 * m4 + l3 * m3 + l4 * m2 + n2;
    i3 = l3 * m4 + l4 * m3 + n3;
    i4 = l4 * m4 + n4;
    l4 = i4 % 4096;
    i3 = i3 + i4 / 4096;
    l3 = i3 % 4096;
    i2 = i2 + i3 / 4096;
    l2 = i2 % 4096;
    l1 = (i1 + i2 / 4096) % 4096;
    r  = twom12 * (l1 + twom12 * (l2 + twom12 * (l3 + twom12 * (l4))));

    return r;
}

void
Random :: SetRandom(int * s, int p1, int p2)
{
    m1 = 502;
    m2 = 1521;
    m3 = 4071;
    m4 = 2107;
    l1 = s[0] % 4096;
    l2 = s[1] % 4096;
    l3 = s[2] % 4096;
    l4 = s[3] % 4096;
    l4 = 2 * (l4 / 2) + 1;
    n1 = 0;
    n2 = 0;
    n3 = p1;
    n4 = p2;
}

// set random using default method specifing a seed file
void
Random :: SetRandom(string iseed)
{
    int seed[4];
    int p1, p2;
    ifstream Primes("../../RC/Primes");

    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    } else { cerr << "ERROR: Unable to open Primes" << endl; }
    Primes.close();

    ifstream input(iseed);
    string property;
    if (input.is_open()) {
        while (!input.eof() ) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                this->SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else { cerr << "ERROR: Unable to open " + iseed << endl; }
}

// Exponential random using inverse CDF method
double
Random :: Exp(double lambda)
{
    return -(1. / lambda) * (log(1 - this->Rannyu()));
}

// Cauchy random using inverse CDF method
double
Random :: Cauchy(double mu, double lambda)
{
    return mu + lambda * tan(M_PI * (this->Rannyu() - 0.5));
}

// Cosine of 2D random angle
double
Random :: Cosine()
{
    pair<double, double> pt;
    double rr;

    do {
        pt.first  = this->Rannyu(-1., 1.);
        pt.second = this->Rannyu(0., 1.);
        rr        = pow(pt.first, 2) + pow(pt.second, 2);
    } while (rr > 1.);

    return (pow(pt.first, 2) - pow(pt.second, 2)) / rr;
}

// Sine of 2D random angle
double
Random :: Sine()
{
    return sqrt(1. - pow(this->Cosine(), 2));
}

// Hit or Miss method 1-dim function
double
Random :: HitMiss1d(function<double(double)> pdf, double pdfM, double a, double b)
{
    pair<double, double> pt;

    do {
        pt.first  = this->Rannyu(a, b);
        pt.second = pdf(pt.first);
    } while (this->Rannyu() > pt.second / pdfM);

    return pt.first;
}

// Return f1(random_unif(a,b))
double
Random :: Rannyuf1d(function<double(double)> f1, double a, double b)
{
    return f1(this->Rannyu(a, b));
}

// Return f1(random(a,b)) distributed as pdf using HitMiss1d
double
Random :: Rannyuf1d(function<double(double)> f1, function<double(double)> pdf, double pdfM, double a, double b)
{
    return f1(this->HitMiss1d(pdf, pdfM, a, b));
}

// Return random discrete number from [a,b]
int
Random :: RannyuDiscrete(int a, int b)
{
    if (a < 0 or b < 0 or b < a) {
        cerr << "RannyuDiscrete error parameters" << endl;
    }

    set<double> intvs;
    int nintv = (b - a + 1);

    intvs.insert(0.);
    intvs.insert(1.);

    for (auto i = 1; i < nintv; ++i) {
        intvs.insert(0. + i * 1. / nintv);
    }

    auto rr = this->Rannyu();
    intvs.insert(rr);
    auto it = find(intvs.begin(), intvs.end(), rr);

    return (int) distance(intvs.begin(), it) - 1 + a;
}

vector<double>
Random :: Metropolis3d(function<double(vector<double>)> pdf, function<vector<double>(vector<double>)> t1,
  vector<double> pt)
{
    double alpha, rr;
    vector<double> ptn;

    while (true) {
        ptn = t1(pt);

        alpha = min(1., pdf(ptn) / pdf(pt));
        rr    = this->Rannyu();

        if (rr <= alpha) /*cout << 1 << endl;*/ break;
        // cout << 0 << endl;
    }

    return ptn;
}

// get random vector of _size = size, distributed unif [0,1]
vector<double>
getRannyu(Random & rnd, unsigned int size)
{
    auto gen = [&rnd](){
          return rnd.Rannyu();
      };

    vector<double> v(size);

    generate(v.begin(), v.end(), gen);

    return v;
}

// get random vector of _size = size, distributed unif [min, max]
vector<double>
getRannyu(Random & rnd, unsigned int size, double min, double max)
{
    auto gen = [&rnd, &min, &max](){
          return rnd.Rannyu(min, max);
      };

    vector<double> v(size);

    generate(v.begin(), v.end(), gen);

    return v;
}

// get random vector of _size = size, distributed gaus(mean, sigma)
vector<double>
getGauss(Random & rnd, unsigned int size, double mean, double sigma)
{
    auto gen = [&rnd, &mean, &sigma](){
          return rnd.Gauss(mean, sigma);
      };

    vector<double> v(size);

    generate(v.begin(), v.end(), gen);

    return v;
}

// get random vector of _size = size, distributed Exp(lambda)
vector<double>
getExp(Random & rnd, unsigned int size, double lambda)
{
    auto gen = [&rnd, &lambda](){
          return rnd.Exp(lambda);
      };

    vector<double> v(size);

    generate(v.begin(), v.end(), gen);

    return v;
}

// get random vector of _size = size, distributed Cauchy(mu, lambda)
vector<double>
getCauchy(Random & rnd, unsigned int size, double mu, double lambda)
{
    auto gen = [&rnd, &mu, &lambda](){
          return rnd.Cauchy(mu, lambda);
      };

    vector<double> v(size);

    generate(v.begin(), v.end(), gen);

    return v;
}

// get random vector of _size = size, cosine of distributed unfi 2D angle [0,2pi]
vector<double>
getCosine(Random & rnd, unsigned int size)
{
    auto gen = [&rnd](){
          return rnd.Cosine();
      };

    vector<double> v(size);

    generate(v.begin(), v.end(), gen);

    return v;
}

// get random vector of _size = size, sine of distributed unfi 2D angle [0,2pi]
vector<double>
getSine(Random & rnd, unsigned int size)
{
    auto gen = [&rnd](){
          return rnd.Sine();
      };

    vector<double> v(size);

    generate(v.begin(), v.end(), gen);

    return v;
}

// get random vector of _size = size, distributed as _pdf inside [a,b]
vector<double>
getHitMiss1d(Random & rnd, unsigned int size, function<double(double)> pdf, double pdfM, double a, double b)
{
    auto gen = [&rnd, &a, &b, &pdf, &pdfM](){
          return rnd.HitMiss1d(pdf, pdfM, a, b);
      };

    vector<double> v(size);

    generate(v.begin(), v.end(), gen);

    return v;
}

// get random vector of _size = size, _f1(random(a,b))
vector<double>
getRannyuf1d(Random & rnd, unsigned int size, function<double(double)> f1, double a, double b)
{
    auto gen = [&rnd, &a, &b, &f1](){
          return rnd.Rannyuf1d(f1, a, b);
      };

    vector<double> v(size);

    generate(v.begin(), v.end(), gen);

    return v;
}

vector<int>
getRannyuDiscrete(Random & rnd, unsigned int size, int a, int b)
{
    auto gen = [&rnd, &a, &b](){
          return rnd.RannyuDiscrete(a, b);
      };

    vector<int> v(size);

    generate(v.begin(), v.end(), gen);

    return v;
}

// get random vector of _size = size, _f1(random(a,b)) distributed as pdf using HitMiss1d
vector<double>
getRannyuf1d(Random & rnd, unsigned int size, function<double(double)> f1, function<double(double)> pdf, double pdfM,
  double a, double b)
{
    auto gen = [&rnd, &a, &b, &f1, &pdf, &pdfM](){
          return rnd.Rannyuf1d(f1, pdf, pdfM, a, b);
      };

    vector<double> v(size);

    generate(v.begin(), v.end(), gen);

    return v;
}

vector<vector<double> >
getMetropolis3d(Random &rnd, unsigned int size, function<double(vector<double>)> pdf,
  function<vector<double>(vector<double>)> t1,
  vector<double> pt)
{
    auto pta = pt;

    auto gen = [&rnd, &pdf, &t1, &pta](){
          pta = rnd.Metropolis3d(pdf, t1, pta);
          return pta;
      };

    vector<vector<double> > v(size);

    generate(v.begin(), v.end(), gen);

    return v;
}

// write the input vector to a file inside out dir
void
writeVector(vector<double> & vct, string ofile)
{
    ofstream ofs(outPath + ofile + ".out");

    ostream_iterator<double> output_iterator(ofs, "\n");

    if (ofs.is_open()) {
        copy(vct.begin(), vct.end(), output_iterator);
    } else { cerr << "ERROR: Unable to open " + ofile + ".out" << endl; }

    ofs.close();
}

/***************************************************************************
* BLOCKING METHOD ***
* TO CHECK NUMERIC SIMULATION CONVERGENCE ***
*
* Write to file _oname-mean-conv.csv inside out dir
*   mean + mean statistical uncertainty convergence of the random
*   vector vct divided in nblk
* Write to file _oname-var-conv.csv inside out dir
*   variance + variance statistical uncertainty convergence of the random
*   vector vct divided in nblk, using mean_teo as expetation value
***************************************************************************/
void
blockingMethod_old(vector<double> & vct, double mean_teor, unsigned int nblk, string oname)
{
    if (vct.size() % nblk != 0) {
        cerr << "ERROR: Throw size not divisible by block size exactly" << endl;
    }

    // number of elements per block
    unsigned int L = vct.size() / nblk;

    /****************************************
    *  the idea is too have dict for every block
    *  dict keys for i-th block
    *  - mean  : r_i
    *  - mean2 : (r_i)^2
    *  - var   : sigma^2_{r_i}
    *  - var2  : (sigma^2_{r_i})^2
    ****************************************/
    vector<map<string, double> > blks(nblk);

    /****************************************
    *  the idea is too have dict for every block
    *  dict keys for i-th block
    *  - mean     : cumulative average
    *  - mean2    : cumulative square average
    *  - err_mean : statistical uncertainty over mean
    *  - var      : cumulative average
    *  - var2     : cumulative square average
    *  - err_var  : statistical uncertainty over variance
    ****************************************/
    vector<map<string, double> > sum_blks(nblk);

    /*******************************************
    * calculate r_i, (r_i)^2
    * calculate sigma^2_{r_i}, (sigma^2_{r_i})^2
    *******************************************/
    auto i(0);
    // accumulate binary lambda operator for blk
    auto acc = [&mean_teor](double a, double b){
          return a + pow((b - mean_teor), 2);
      };
    for (auto & blk : blks) {
        blk["mean"]  = accumulate(vct.begin() + i * L, vct.begin() + (i + 1) * L, 0.0d) / L;
        blk["mean2"] = pow(blk["mean"], 2);
        blk["var"]   = accumulate(vct.begin() + i * L, vct.begin() + (i + 1) * L, 0.0d, acc) / L;
        blk["var2"]  = pow(blk["var"], 2);
        i += 1;
    }

    /****************************************
    *  calculate
    *    cumulative average,
    *    cumulative square average,
    *    statistical uncertainty on mean
    *  calculate
    *    cumulative average,
    *    cumulative square average,
    *    statistical uncertainty on variance
    ****************************************/
    i = 0;

    for (auto & sum_blk : sum_blks) {
        for (const auto pair : *blks.begin()) {
            // accumulate binary lambda operator for sum_blk
            auto acc2 = [&pair](double a, map<string, double> & b){
                  return a + b[pair.first];
              };
            sum_blk[pair.first] =
              accumulate(blks.begin(), blks.begin() + (i + 1), 0.0d, acc2) / (i + 1);
        }
        // better to use (i+1)-1 = i cause (i+1)->nblk
        if (i == 0) {
            sum_blk["err_mean"] = 0;
            sum_blk["err_var"]  = 0;
        } else {
            sum_blk["err_mean"] = sqrt((sum_blk["mean2"] - pow(sum_blk["mean"], 2)) / i);
            sum_blk["err_var"]  = sqrt((sum_blk["var2"] - pow(sum_blk["var"], 2)) / i);
        }
        i += 1;
    }

    string of1(oname + "-mean-conv.csv");
    string of2(oname + "-var-conv.csv");
    ofstream ofs(outPath + of1);

    if (ofs.is_open()) {
        for (auto & sum_blk : sum_blks) {
            ofs << sum_blk["mean"] << "," << sum_blk["err_mean"] << endl;
        }
    } else { cerr << "ERROR: Unable to open " + of1 << endl; }

    ofs.close();
    ofs = ofstream(outPath + of2);

    if (ofs.is_open()) {
        for (auto & sum_blk : sum_blks) {
            ofs << sum_blk["var"] << "," << sum_blk["err_var"] << endl;
        }
    } else { cerr << "ERROR: Unable to open " + of1 << endl; }
    ofs.close();
} // blockingMethod_old

/***************************************************************************
* BLOCKING METHOD ***
* TO CHECK NUMERIC SIMULATION CONVERGENCE ***
*
* Input:
* - _lbf, function of the expirement to simulate
* - _nthrow, number of random generation
* - _nblk, number of blocks or numer of simulations
* - _oname, name of the file csv in output
***************************************************************************/
void
blockingMethod(function<double(unsigned int)> lbf, unsigned int nthrow, unsigned int nblk, string oname)
{
    if (nthrow % nblk != 0) {
        cerr << "ERROR: Throw size not divisible by block size exactly" << endl;
    }

    string of1(oname + "-blk.csv");
    ofstream ofs(outPath + of1);

    vector<double> vals(nblk);      // r_i
    vector<double> means(nblk);     // cumulative average
    vector<double> errs(nblk);      // statistical uncertainty
    unsigned int L = nthrow / nblk; // generation per block

    // expirement simulator
    auto gen = [&lbf, &L](){
          return lbf(L);
      };

    generate(vals.begin(), vals.end(), gen); // calculate r_i
    vector<double> val2s(vals);              // (r_i)^2
    transform(val2s.begin(), val2s.end(), val2s.begin(), [](double & a){
        return pow(a, 2);
    });

    // calculate and write to file cumulative average, statistical uncertainty
    double sum_val, sum_val2;
    for (unsigned int i = 0; i < nblk; ++i) {
        sum_val  = accumulate(vals.begin(), vals.begin() + (i + 1), 0.0d) / (i + 1);
        sum_val2 = accumulate(val2s.begin(), val2s.begin() + (i + 1), 0.0d) / (i + 1);
        means[i] = sum_val;
        if (i == 0) {
            errs[i] = 0;
        } else {
            errs[i] = sqrt((sum_val2 - pow(sum_val, 2)) / i);
        }

        if (ofs.is_open()) {
            ofs << means[i] << "," << errs[i] << endl;
        } else { cerr << "ERROR: Unable to open " + of1 << endl; }
    }
} // blockingMethod

/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/
