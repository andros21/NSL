/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/

#include "moldyn.h"

static string outPath = "../out/";

MolDyn :: MolDyn(string state, bool restart, unsigned int ndim) : _state(state), _ndim(ndim), _seed(1), _r(_ndim), _ro(
        _ndim), _v(
        _ndim)
{
    // Read input file of particular state
    string input_file(_state + ".input");
    ifstream ReadInput(input_file);

    if (ReadInput.is_open()) {
        ReadInput >> _temp;
        ReadInput >> _npart;
        ReadInput >> _rho;
        _vol = (double) _npart / _rho;
        _box = pow(_vol, 1.0 / _ndim);
        ReadInput >> _rcut;
        ReadInput >> _delta;
        ReadInput >> _nbin;
        ReadInput >> _nstep;
        ReadInput >> _pprint;
        ReadInput >> _psave;
    } else { cerr << "ERROR: Unable to open " + input_file << endl; }
    ReadInput.close();

    // Get radius range
    _rrange = this->getRadiusRange();

    // Decide to restart from previous simulation or randomize
    if (restart) {
        this->Restart(); // Restart from previous simulation
    } else {
        cout << "INFO: Use random vel for start conf " + _state << endl;
        string conf_file = "../../MD/fcc.tsv"; // Take fcc as starting conf
        this->setPosFromFile(conf_file);       // Take fcc as starting conf
        this->setVelRandom();                  // Set Random Velocities
        this->setPosOldFromVel();              // Setold positions from velocities
        this->Measure();
    }
}

// NULL decostractor
MolDyn :: ~MolDyn(){ }

// Write last configuration _r and lastlast configuration _ro
void
MolDyn :: End()
{
    string ofile("r-" + _state + ".tsv");
    ofstream WriteConf(ofile);

    // Write _r  at the end
    if (WriteConf.is_open()) {
        for (unsigned int j = 0; j < _npart; ++j) {
            for (unsigned int i = 0; i < _ndim; ++i) {
                if (i == _ndim - 1) WriteConf << _r[i][j] / _box;
                else WriteConf << _r[i][j] / _box << "\t";
            }
            WriteConf << endl;
        }
    } else { cerr << "ERROR: Unable to open " + ofile << endl; }
    WriteConf.close();


    ofile = "ro-" + _state + ".tsv";
    WriteConf.open(ofile, ios::out);
    // Write _ro at the end
    if (WriteConf.is_open()) {
        for (unsigned int j = 0; j < _npart; ++j) {
            for (unsigned int i = 0; i < _ndim; ++i) {
                if (i == _ndim - 1) WriteConf << _ro[i][j] / _box;
                else WriteConf << _ro[i][j] / _box << "\t";
            }
            WriteConf << endl;
        }
    } else { cerr << "ERROR: Unable to open " + ofile << endl; }
    WriteConf.close();
}

// Restart from previous simulation and adjust it to match _temp
void
MolDyn :: Restart()
{
    cout << "INFO: Use last conf as start conf " + _state << endl;

    // Read from previous simulation
    string conf_file("r-" + _state + ".tsv");

    this->setPosFromFile(conf_file);

    // Read from previous simulation
    conf_file = "ro-" + _state + ".tsv";
    this->setPosOldFromFile(conf_file);

    // Iterate 7 times the restart -> simulation(1000)
    for (auto i = 0; i < 7; ++i) {
        this->Simulate(1000);
        this->setVelReFactor();
        this->setPosOldFromVel();
    }

    // set current meausure
    this->Measure();
}

// Return a copy of positions
vector<vector<double> >
MolDyn :: getPos()
{
    return _r;
}

// Return a copy of old positions
vector<vector<double> >
MolDyn :: getPosOld()
{
    return _ro;
}

// Return a copy of velocities
vector<vector<double> >
MolDyn :: getVel()
{
    return _v;
}

// Load conf vector<(_ndim)vector(_npart)> from file
void
MolDyn :: setVctFromFile(vector<vector<double> > & vct, string ifile)
{
    for (unsigned int i = 0; i < _ndim; ++i) {
        if (vct[i].size() != 0)
            vct[i].clear();
    }
    ifstream ReadConf(ifile);

    if (ReadConf.is_open()) {
        double appo;
        ReadConf >> appo;
        while (!ReadConf.eof()) {
            for (unsigned int i = 0; i < _ndim; ++i) {
                vct[i].push_back(appo * _box);
                ReadConf >> appo;
            }
        }
    } else { cerr << "ERROR: Unable to open " + ifile << endl; }
    ReadConf.close();
}

// Load pos vector<(_ndim)vector(_npart)> from file
void
MolDyn :: setPosFromFile(string ifile)
{
    this->setVctFromFile(this->_r, ifile);
}

// Load pos old vector<(_ndim)vector(_npart)> from file
void
MolDyn :: setPosOldFromFile(string ifile)
{
    this->setVctFromFile(this->_ro, ifile);
}

// Initialize random velocities and refactor to reach the correct _temp
void
MolDyn :: setVelRandom()
{
    srand(_seed);
    vector<double> empty(_npart);

    // Randomize velocities
    for (unsigned int i = 0; i < _ndim; ++i) {
        _v[i] = empty;
        auto gen = [](){
              return rand() / double(RAND_MAX) - 0.5;
          };
        generate(_v[i].begin(), _v[i].end(), gen);
    }

    // Adjust velocities drift
    for (unsigned int i = 0; i < _ndim; ++i) {
        double drift = accumulate(_v[i].begin(), _v[i].end(), 0.0d) / _npart;
        auto tra = [&drift](double & el){
              return el - drift;
          };
        transform(_v[i].begin(), _v[i].end(), _v[i].begin(), tra);
    }

    // Refactor velocity as function of _temp
    this->setVelReFactor();
} // MolDyn::setRandomVel

// Refactor _v as function of _temp
void
MolDyn :: setVelReFactor()
{
    double fs = sqrt(_ndim * _temp / (this->getNormVel() / (double) _npart));

    for (unsigned int i = 0; i < _ndim; ++i) {
        auto tra = [&fs](double & el){
              return el * fs;
          };
        transform(_v[i].begin(), _v[i].end(), _v[i].begin(), tra);
    }
}

// Derive ignote _ro from _r and _v
void
MolDyn :: setPosOldFromVel()
{
    for (unsigned int i = 0; i < _ndim; ++i) {
        _ro[i] = _r[i];
        auto tra = [this](double & rr, double & vv){
              return this->getPbc(rr - vv * _delta);
          };
        transform(_ro[i].begin(), _ro[i].end(), _v[i].begin(), _ro[i].begin(), tra);
    }
}

// Get sum of vel norm
double
MolDyn :: getNormVel()
{
    double sum = 0.;
    auto acc = [](double & a, double & b){
          return a + pow(b, 2);
      };

    for (unsigned int i = 0; i < _ndim; ++i) {
        sum += accumulate(_v[i].begin(), _v[i].end(), 0.0d, acc);
    }
    return sum;
}

vector<double>
MolDyn :: getRadiusRange()
{
    vector<double> vct(_nbin + 1);

    for (double i = 0.; i <= _nbin; ++i)
        vct[i] = (i * _box * 0.5 / (double) _nbin);

    return vct;
}

double
MolDyn :: getPbc(double r)
{
    return r - _box * rint(r / _box);
}

// Make an step using Verlet algorithm
void
MolDyn :: Move()
{
    // new positions
    vector<vector<double> > rf(_ndim);

    // lambda function for transformation
    auto minus2 = [](double & a, double & b){
          return 2. * a - b;
      };

    auto plus2 = [this](double & a, double & b){
          return this->getPbc(a + b * pow(_delta, 2));
      };

    auto tra = [this](double & a, double & b){
          return this->getPbc(a - b) / (2. * _delta);
      };

    // positions new from positions and positions old
    for (unsigned int i = 0; i < _ndim; ++i) {
        rf[i] = _r[i];
        transform(rf[i].begin(), rf[i].end(), _ro[i].begin(), rf[i].begin(), minus2);
    }

    auto ff = this->Force();

    for (unsigned int i = 0; i < _ndim; ++i) {
        transform(rf[i].begin(), rf[i].end(), ff[i].begin(), rf[i].begin(), plus2);
    }

    // velocity new from positions new and positions old
    for (unsigned int i = 0; i < _ndim; ++i) {
        _v[i] = rf[i];
        transform(_v[i].begin(), _v[i].end(), _ro[i].begin(), _v[i].begin(), tra);
    }


    // positions old are now positions and positions are now positions new
    for (unsigned int i = 0; i < _ndim; ++i) {
        _ro[i] = _r[i];
        _r[i]  = rf[i];
    }
} // MolDyn::Move

// Make a meausure
void
MolDyn :: Measure(bool corl)
{
    _th['K'] = 0.5 * this->getNormVel() / (double) _npart; // Kinetic energy per particle
    _th['T'] = (2.0 / (double) _ndim) * _th['K'];          // Temperature
    _th['V'] = this->getPotential(corl) / (double) _npart; // Potential energy per particle
    _th['E'] = _th['K'] + _th['V'];                        // Total energy per paticle
}

void
MolDyn :: Simulate(int nstep, bool corl, bool verbose, bool result, bool xyz)
{
    int nconf(1);

    if (nstep == 0) nstep = this->_nstep;

    for (int istep = 1; istep <= nstep; ++istep) {
        this->Move();
        if ((istep % (int) (_pprint * nstep) == 0) && verbose) cout << "Number of time-steps: " << istep << endl;
        if ((istep % (int) (_psave * nstep) == 0) && result) {
            this->Measure(corl);
            this->getMeasureToFile();
            if (xyz) this->writeVctToFile(_r, nconf);
            ++nconf;
        }
    }
    if (result != true) this->Measure();
}

// Return avaiable var meausurament
string
MolDyn :: getAvailableVar()
{
    string vars = "";

    for (auto & pp : _th)
        vars += pp.first;

    return vars;
}

// Return a particular meausure
double
MolDyn :: getMeasure(char var)
{
    return _th[var];
}

void
MolDyn :: getMeasureToFile(string var)
{
    if (var == "") var = this->getAvailableVar();
    string ofile(outPath + "th-" + _state + ".tsv");
    ofstream WriteResult(ofile, ios::app);

    if (WriteResult.is_open()) {
        for (auto car : var) {
            if (car == var[var.size() - 1]) WriteResult << _th[car];
            else WriteResult << _th[car] << "\t";
        }
        WriteResult << endl;
    } else { cerr << "ERROR: Unable to open " + ofile << endl; }
    WriteResult.close();
}

// Compute forces as vector<(_ndim)vector(_npart)>
vector<vector<double> >
MolDyn::Force()
{
    vector<vector<double> > ff(_ndim);
    vector<double> dr(_ndim, 0.);
    vector<double> fr(_ndim, 0.);
    double sdr, sfr;

    for (unsigned int j = 0; j < _npart; ++j) {
        for (unsigned int i = 0; i < _ndim; ++i)
            fr[i] = 0.;
        for (unsigned int k = 0; k < _npart; ++k) {
            if (j != k) {
                for (unsigned int i = 0; i < _ndim; ++i)
                    dr[i] = this->getPbc(_r[i][j] - _r[i][k]);
                sdr = sqrt(accumulate(dr.begin(), dr.end(), 0.0d, [](double &a, double &b){
                    return a + pow(b, 2);
                }));
                if (sdr < _rcut) {
                    sfr = (48.0 / pow(sdr, 14) - 24.0 / pow(sdr, 8));
                    for (unsigned int i = 0; i < _ndim; ++i)
                        fr[i] += dr[i] * sfr;
                }
            }
        }
        for (unsigned int i = 0; i < _ndim; ++i)
            ff[i].push_back(fr[i]);
    }

    return ff;
}

// Compute potential energy
double
MolDyn::getPotential(bool corl)
{
    double vv(0.);
    vector<double> dr(_ndim, 0.);
    double sdr;
    unsigned int cf;
    vector<double> gg(_nbin, 0.0d);

    for (unsigned int j = 0; j < _npart - 1; ++j) {
        for (unsigned int k = j + 1; k < _npart; ++k) {
            for (unsigned int i = 0; i < _ndim; ++i)
                dr[i] = this->getPbc(_ro[i][j] - _ro[i][k]);
            sdr = sqrt(accumulate(dr.begin(), dr.end(), 0.0d, [](double &a, double &b){
                return a + pow(b, 2);
            }));
            if (corl) {
                // update histo of g(r) [not normalized]
                if (sdr < _box * 0.5) {
                    cf      = (unsigned int) ((sdr * _nbin) / (_box * 0.5));
                    gg[cf] += 2;
                }
            }
            if (sdr < _rcut) vv += 4.0 / pow(sdr, 12) - 4.0 / pow(sdr, 6);
        }
    }

    if (corl) _gg = gg;

    return vv;
}

void
MolDyn :: writeVctToFile(vector<vector<double> > vct, int blk, bool pbc)
{
    auto nu = to_string(blk);

    while (nu.length() < to_string(_nstep / (int) (_psave * _nstep)).size() - 1) {
        nu = "0" + nu;
    }

    auto xyz_file = outPath + "xyz/" + nu + ".xyz";
    ofstream WriteXYZ(xyz_file);

    if (WriteXYZ.is_open()) {
        if (pbc) {
            for (unsigned int j = 0; j < _npart; ++j) {
                for (unsigned int i = 0; i < _ndim; ++i) {
                    if (i == _ndim - 1) WriteXYZ << this->getPbc(vct[i][j]);
                    else WriteXYZ << this->getPbc(vct[i][j]) << "\t";
                }
                WriteXYZ << endl;
            }
        } else {
            for (unsigned int j = 0; j < _npart; ++j) {
                for (unsigned int i = 0; i < _ndim; ++i)
                    WriteXYZ << vct[i][j] << "\t";
                WriteXYZ << endl;
            }
        }
    } else { cerr << "ERROR: Unable to open " + xyz_file << endl; }
    WriteXYZ.close();
}

void
MolDyn :: blockingMethod(unsigned int nblk, bool corl, bool bth)
{
    if (_nstep % nblk != 0) {
        cerr << "ERROR: Throw size not divisible by block size exactly" << endl;
    }

    map<char, vector<double> > th;
    vector<vector<double> > gg(_nbin);
    string vars(this->getAvailableVar());
    double ggnorm;

    for (auto i = 1; i <= _nstep; ++i) {
        this->Move();
        this->Measure(corl);
        if (bth) {
            for (auto ch : vars)
                th[ch].push_back(this->getMeasure(ch));
        }
        if (corl) {
            // load and normalized
            for (int n = 0; n < _nbin; ++n) {
                ggnorm = ((4. * M_PI) / 3) * _npart * _rho * (pow(_rrange[n + 1], 3) - pow(_rrange[n], 3));
                gg[n].push_back(_gg[n] / ggnorm);
            }
        }
    }
    if (bth) {
        for (auto ch : vars) {
            string of1("-blk.csv");
            of1 = _state + of1;
            of1 = "-" + of1;
            of1 = ch + of1;
            ofstream ofs(outPath + of1);

            unsigned int L = _nstep / nblk; // generation per block
            vector<double> vals(nblk);      // r_i

            for (unsigned int i = 0; i < nblk; ++i)
                vals[i] = (double) accumulate(th[ch].begin() + i * L, th[ch].begin() + (i + 1) * L, 0.0d) / L;

            vector<double> means(nblk); // cumulative average
            vector<double> errs(nblk);  // statistical uncertainty

            vector<double> val2s(vals); // (r_i)^2
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
            ofs.close();
        }
    }
    if (corl) {
        string of1("-blk.csv");
        of1 = _state + of1;
        of1 = "-" + of1;
        of1 = "gr" + of1;
        ofstream ofs(outPath + of1);
        for (int n = 0; n < _nbin; ++n) {
            // ofstream ofs(outPath + of1, ios::app);
            unsigned int L = _nstep / nblk; // generation per block
            vector<double> vals(nblk);      // r_i

            for (unsigned int i = 0; i < nblk; ++i)
                vals[i] = (double) accumulate(gg[n].begin() + i * L, gg[n].begin() + (i + 1) * L, 0.0d) / L;

            vector<double> means(nblk); // cumulative average
            vector<double> errs(nblk);  // statistical uncertainty

            vector<double> val2s(vals); // (r_i)^2
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
            }
            if (ofs.is_open()) {
                ofs << _rrange[n] << "," << means[nblk - 1] << "," << errs[nblk - 1] << endl;
            } else { cerr << "ERROR: Unable to open " + of1 << endl; }
            // ofs.close();
        }
        ofs.close();
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
