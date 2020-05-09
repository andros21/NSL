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
#include <map>

using namespace std;


#ifndef __MolDyn__
# define __MolDyn__


class MolDyn {
private:
    string _state;                         // state of matter
    unsigned int _npart, _ndim;            // particle number, dimension number
    double _temp, _vol, _rho, _box, _rcut; // thermodynamical state
    int _nstep, _seed;                     // simulation variable
    double _delta, _pprint, _psave;        // simulation variable
    vector<vector<double> > _r;            // positional configuration
    vector<vector<double> > _ro;           // positional configuration step before
    vector<vector<double> > _v;            // velocity configuration
    map<char, double> _th;                 // thermodynamical variable results

public:
    MolDyn(string state = "solid", bool restart = false, unsigned int ndim = 3);
    ~MolDyn();
    void
    End();
    void
    Restart();
    vector<vector<double> >
    getPos();
    vector<vector<double> >
    getPosOld();
    vector<vector<double> >
    getVel();
    void
    setVctFromFile(vector<vector<double> > & vct, string ifile);
    void
    setPosFromFile(string ifile);
    void
    setPosOldFromFile(string ifile);
    void
    setVelRandom();
    void
    setVelReFactor();
    void
    setPosOldFromVel();
    double
    getNormVel();
    double
    getPotential();
    double
    getPbc(double r);
    void
    Move();
    void
    Measure();
    void
    Simulate(int nstep = 0, bool verbose = false, bool result = false, bool xyz = false);
    string
    getAvailableVar();
    double
    getMeasure(char var);
    void
    getMeasureToFile(string var = "");
    vector<vector<double> >
    Force();
    void
    writeVctToFile(vector<vector<double> > vct, int blk, bool pbc = true);
    void
    blockingMethod(unsigned int nblk);
};


#endif // ifndef __MolDyn__

/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/
