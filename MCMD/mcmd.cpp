/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/

#include "mcmd.h"

static string outPath = "../out/";

McMd :: McMd(string state, bool restart, unsigned int ndim) : _state(state), _ndim(ndim),
    _r(_ndim), Macpt(0), Mttot(0)
{
    // set random
    _rnd.SetRandom();

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

    // Tail corrections for potential energy and pressure
    _vtail = (8.0 * M_PI * _rho) / (9.0 * pow(_rcut, 9)) - (8.0 * M_PI * _rho) / (3.0 * pow(_rcut, 3));
    _ptail = (32.0 * M_PI * _rho) / (9.0 * pow(_rcut, 9)) - (16.0 * M_PI * _rho) / (3.0 * pow(_rcut, 3));

    // Get radius range
    _rrange = this->getRadiusRange();

    // Decide to restart from previous simulation or randomize
    if (restart) {
        this->Restart(); // Restart from previous simulation
    } else {
        cout << "INFO: Use fcc for start conf " + _state << endl;
        string conf_file = "../../MCMD/fcc.tsv"; // Take fcc as starting conf
        this->setPosFromFile(conf_file);         // Take fcc as starting conf
        this->Measure();
    }
}

// NULL decostractor
McMd :: ~McMd(){ }

// Write last configuration _r
void
McMd :: End()
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
}

// Restart from previous simulation and adjust it to match _temp
void
McMd :: Restart()
{
    cout << "INFO: Use last conf as start conf " + _state << endl;

    // Read from previous simulation
    string conf_file("r-" + _state + ".tsv");

    this->setPosFromFile(conf_file);

    // set current meausure
    this->Measure();
}

// Return a copy of positions
vector<vector<double> >
McMd :: getPos()
{
    return _r;
}

// Load conf vector<(_ndim)vector(_npart)> from file
void
McMd :: setVctFromFile(vector<vector<double> > & vct, string ifile)
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
McMd :: setPosFromFile(string ifile)
{
    this->setVctFromFile(this->_r, ifile);
}

double
McMd :: getPbc(double r)
{
    return r - _box * rint(r / _box);
}

vector<double>
McMd :: getRadiusRange()
{
    vector<double> vct(_nbin + 1);

    for (double i = 0.; i <= _nbin; ++i)
        vct[i] = (i * _box * 0.5 / (double) _nbin);

    return vct;
}

double
McMd :: getMetroRatio()
{
    if (Mttot == 0) {
        return 0.;
    }
    return (double) Macpt / Mttot;
}

void
McMd :: coutMetroRatio()
{
    cout << "\rMetropolis Ratio: " << round(this->getMetroRatio() * 100) << " " << flush;
}

void
McMd :: resetMetroRatio()
{
    Macpt = 0;
    Mttot = 0;
}

double
McMd :: getEnergyGap(vector<double> & vct, unsigned int ix)
{
    double dr;
    double ene(0.);
    vector<double> dd(_ndim);

    for (unsigned int j = 0; j < _npart; ++j) {
        if (j != ix) {
            // distance ix-i in pbc
            for (unsigned int i = 0; i < _ndim; ++i)
                dd[i] = this->getPbc(vct[i] - _r[i][j]);

            dr = sqrt(accumulate(dd.begin(), dd.end(), 0.0d, [](double & a, double & b){
                return a + pow(b, 2);
            }));

            if (dr < _rcut) {
                ene += 1.0 / pow(dr, 12) - 1.0 / pow(dr, 6);
            }
        }
    }

    return 4.0 * ene;
}

// Make an step using Metropolis Algorithm
void
McMd :: Move()
{
    double p;
    unsigned int o;
    vector<double> rold(_ndim), rnew(_ndim);

    for (unsigned int j = 0; j < _npart; ++j) {
        // Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
        o = (unsigned int) _rnd.Rannyu(0, _npart - 1);

        for (unsigned int i = 0; i < _ndim; ++i) {
            rold[i] = _r[i][o];                                                 // Old
            rnew[i] = this->getPbc(_r[i][o] + _delta * (_rnd.Rannyu() - 0.5) ); // New
        }

        // Metropolis test
        p = exp((1. / _temp) * (this->getEnergyGap(rold, o) - this->getEnergyGap(rnew, o)));
        Mttot++;
        if (p >= _rnd.Rannyu()) {
            Macpt++;
            // Update
            for (unsigned int i = 0; i < _ndim; ++i)
                _r[i][o] = rnew[i];
        }
    }
} // McMd::Move

// Make a meausure
void
McMd :: Measure(bool corl)
{
    double v(0.);
    double w(0.);
    double dr;
    unsigned int cf;
    vector<double> dd(_ndim);
    vector<double> gg(_nbin, 0.0d);

    // cycle over pairs of particles
    for (unsigned int j = 0; j < _npart - 1; ++j) {
        for (unsigned int k = j + 1; k < _npart; ++k) {
            // distance i-j in pbc
            for (unsigned int i = 0; i < _ndim; ++i) {
                dd[i] = this->getPbc(_r[i][j] - _r[i][k]);
            }

            dr = sqrt(accumulate(dd.begin(), dd.end(), 0.0d, [](double & a, double & b){
                return a + pow(b, 2);
            }));
            if (corl) {
                // update histo of g(r) [not normalized]
                if (dr < _box * 0.5) {
                    cf      = (unsigned int) ((dr * _nbin) / (_box * 0.5));
                    gg[cf] += 2;
                }
            }

            if (dr < _rcut) {
                v += 1.0 / pow(dr, 12) - 1.0 / pow(dr, 6);
                w += 1.0 / pow(dr, 12) - 0.5 / pow(dr, 6);
            }
        }
    }

    if (corl) _gg = gg;
    _th['V'] = (4.0 * v) / (double) _npart + _vtail;
    _th['P'] = _rho * _temp + ((48.0 * w / 3.0) + (double) _npart * _ptail) / _vol;
} // McMd::Measure

void
McMd :: Simulate(int nstep, bool corl, bool verbose, bool result, bool xyz)
{
    int nconf(1);

    if (nstep == 0) nstep = this->_nstep;

    for (int istep = 1; istep <= nstep; ++istep) {
        this->Move();
        if ((istep % (int) (_pprint * nstep) == 0) && verbose) {
            cout << "\r" << _state << ": " << "Number of time-steps: " << istep << ", " << "Metropolis Ratio: "
                 << round(
                this->getMetroRatio() * 100) << " " << flush;
        }
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
McMd :: getAvailableVar()
{
    string vars = "";

    for (auto & pp : _th)
        vars += pp.first;

    return vars;
}

// Return a particular meausure
double
McMd :: getMeasure(char var)
{
    return _th[var];
}

void
McMd :: getMeasureToFile(string var)
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

void
McMd :: writeVctToFile(vector<vector<double> > vct, int blk, bool pbc)
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
McMd :: blockingMethod(unsigned int nblk, bool file, bool corl, bool bth)
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
            if (file) {
                string of2 = "-" + _state;
                of2 = ch + of2;
                writeVector(th[ch], of2);
            }
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
