/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/

#include "i1d.h"

static string outPath = "../out/";

I1D :: I1D(char method) : _rnd(), _method(method)
{
    // Set random default seed
    _rnd.SetRandom();

    // Read input file
    string input_file("input.dat");
    ifstream ReadInput(input_file);

    if (ReadInput.is_open()) {
        ReadInput >> _temp;
        ReadInput >> _nspin;
        ReadInput >> _JJ;
        ReadInput >> _hh;
        ReadInput >> _nblk;
        ReadInput >> _nstep;
    } else { cerr << "ERROR: Unable to open " + input_file << endl; }
    ReadInput.close();

    // Set random spin configuration
    this->setRandomSpin();
}

I1D :: ~I1D(){ };

vector<int>
I1D :: getSpin(){ return _mu; }

char
I1D :: getMethod(){ return _method; }

double
I1D :: getTemp(){ return _temp; }

unsigned int
I1D :: getPbc(int i)
{
    if (abs(i) >= _nspin) i = i - _nspin;
    else if (i < 0) i = i + _nspin;
    return i;
}

string
I1D :: getAvailableVar()
{
    string vars = "";

    for (auto & pp : _th)
        vars += pp.first;

    return vars;
}

double
I1D :: getMeasure(char var)
{
    return _th[var];
}

double
I1D :: getEnergyGap(unsigned int i)
{
    return 2 * _JJ * _mu[i] * ( _mu[this->getPbc(i - 1)] + _mu[this->getPbc(i + 1)] ) + 2 * _hh * _mu[i];
}

double
I1D :: getFieldh(){ return _hh; }

void
I1D :: setTemp(double temp){ _temp = temp; }

void
I1D :: setRandomSpin()
{
    _mu = getRannyuDiscrete(_rnd, _nspin);
    replace_if(_mu.begin(), _mu.end(), [](int & i){
        return i == 0;
    }, -1);
    this->Measure();
}

void
I1D :: setFieldh(double h){ _hh = h; }

void
I1D :: setMethod(char method){ if (method == 'G' or method == 'M') _method = method; }

void
I1D :: Equilibrate(unsigned int lowstep)
{
    for (unsigned int i = 0; i < lowstep; ++i)
        this->Move();
    this->Measure();
}

void
I1D :: Move()
{
    unsigned int k;

    if (_method == 'M') {
        double aa;
        k  = _rnd.RannyuDiscrete(0, _nspin - 1);
        aa = min(1., exp(-(1. / _temp) * this->getEnergyGap(k)));
        if (_rnd.Rannyu() < aa) _mu[k] = -_mu[k];
    }
    if (_method == 'G') {
        double pp;
        k  = _rnd.RannyuDiscrete(0, _nspin - 1);
        pp = 1. / (1. + exp(-(1. / _temp) * (this->getEnergyGap(k) / _mu[k])));
        if (_rnd.Rannyu() <= pp) _mu[k] = 1; else _mu[k] = -1;
    }
}

void
I1D :: Measure()
{
    double appo;
    double uu(0.);
    double cc(0.);
    double xx(0.);

    for (unsigned int i = 0; i < _nspin; ++i) {
        appo = -_JJ * _mu[i] * _mu[this->getPbc(i + 1)] - 0.5 * _hh * (_mu[i] + _mu[this->getPbc(i + 1)]);
        uu  += appo;
        cc  += pow(appo, 2);
        xx  += _mu[i];
    }

    _th['U'] = uu / (double) _nspin;
    _th['C'] = pow(1. / _temp, 2) * (cc / (double) _nspin - pow(_th['U'], 2));
    _th['X'] = (1. / _temp) * pow(xx, 2) / (double) _nspin;
    _th['M'] = xx / (double) _nspin;
}

void
I1D :: blockingMethod()
{
    if (_nstep % _nblk != 0) {
        cerr << "ERROR: Throw size not divisible by block size exactly" << endl;
    }

    map<char, vector<double> > th;
    string vars(this->getAvailableVar());

    for (unsigned int i = 1; i <= _nstep; ++i) {
        this->Move();
        this->Measure();
        for (auto ch : vars)
            th[ch].push_back(this->getMeasure(ch));
    }
    for (auto ch : vars) {
        string of1("-blk.csv");
        if (_hh != 0.) of1 = "-h" + of1;
        of1 = ch + of1;
        of1 = "-" + of1;
        of1 = _method + of1;
        ofstream ofs(outPath + of1, ios::app);

        unsigned int L = _nstep / _nblk; // generation per block
        vector<double> vals(_nblk);      // r_i

        for (unsigned int i = 0; i < _nblk; ++i)
            vals[i] = (double) accumulate(th[ch].begin() + i * L, th[ch].begin() + (i + 1) * L, 0.0d) / L;

        vector<double> means(_nblk); // cumulative average
        vector<double> errs(_nblk);  // statistical uncertainty

        vector<double> val2s(vals); // (r_i)^2
        transform(val2s.begin(), val2s.end(), val2s.begin(), [](double & a){
            return pow(a, 2);
        });

        // calculate and write to file cumulative average, statistical uncertainty
        double sum_val, sum_val2;
        for (unsigned int i = 0; i < _nblk; ++i) {
            sum_val  = accumulate(vals.begin(), vals.begin() + (i + 1), 0.0d) / (i + 1);
            sum_val2 = accumulate(val2s.begin(), val2s.begin() + (i + 1), 0.0d) / (i + 1);
            means[i] = sum_val;
            if (i == 0) {
                errs[i] = 0;
            } else {
                errs[i] = sqrt((sum_val2 - pow(sum_val, 2)) / i);
            }
        }
        if (ofs.is_open()) {
            ofs << this->getTemp() << "," << means[_nblk - 1] << "," << errs[_nblk - 1] << endl;
        } else { cerr << "ERROR: Unable to open " + of1 << endl; }
        ofs.close();
    }
} // blockingMethod

vector<double>
getTempRange(double ti, double tf, unsigned int size)
{
    vector<double> tv(size + 1);

    tv[0] = ti;
    for (unsigned int i = 1; i < tv.size(); ++i)
        tv[i] = tv[i - 1] + (tf - ti) / (double) size;

    return tv;
}

/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/
