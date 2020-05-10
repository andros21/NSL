/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/

#include "../RC/random.h"

using namespace std;


#ifndef __McMd__
# define __McMd__


class McMd {
private:
    Random _rnd;                                           // random obj
    string _state;                                         // state of matter
    unsigned int _npart, _ndim;                            // particle number, dimension number
    double _temp, _vol, _rho, _box, _rcut, _vtail, _ptail; // thermodynamical state
    int _nstep, _nbin;                                     // simulation variable
    double _delta, _pprint, _psave;                        // simulation variable
    vector<double> _rrange;                                // simulation variable
    vector<vector<double> > _r;                            // positional configuration
    map<char, double> _th;                                 // thermodynamical variable results
    vector<double> _gg;                                    // correlation variable
    unsigned int Macpt, Mttot;                             // Metropolis check vars

public:
    McMd(string state = "solid", bool restart = false, unsigned int ndim = 3);
    ~McMd();
    void
    End();
    void
    Restart();
    vector<vector<double> >
    getPos();
    void
    setVctFromFile(vector<vector<double> > & vct, string ifile);
    void
    setPosFromFile(string ifile);
    double
    getPbc(double r);
    double
    getEnergyGap(vector<double> & vct, unsigned int ix);
    vector<double>
    getRadiusRange();
    double
    getMetroRatio();
    void
    coutMetroRatio();
    void
    resetMetroRatio();
    void
    Move();
    void
    Measure(bool corl = false);
    void
    Simulate(int nstep = 0, bool corl = false, bool verbose = false, bool result = false, bool xyz = false);
    string
    getAvailableVar();
    double
    getMeasure(char var);
    void
    getMeasureToFile(string var = "");
    void
    writeVctToFile(vector<vector<double> > vct, int blk, bool pbc = true);
    void
    blockingMethod(unsigned int nblk, bool file =  false, bool corl = false);
};


#endif // ifndef __McMd__

/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/
