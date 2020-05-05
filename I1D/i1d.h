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


#ifndef __I1D__
# define __I1D__


class I1D {
private:
    Random _rnd;                // Random object
    char _method;               // simulation method: Metropolis or Gibbs
    unsigned int _nspin;        // particle number
    double _JJ, _temp, _hh;     // simulation variable
    unsigned int _nstep, _nblk; // simulation variable
    vector<int> _mu;            // spin configuration
    map<char, double> _th;      // thermodynamical variable results

public:
    I1D(char method = 'M');
    ~I1D();
    vector<int>
    getSpin();
    char
    getMethod();
    double
    getTemp();
    unsigned int
    getPbc(int i);
    double
    getEnergyGap(unsigned int i);
    double
    getFieldh();
    void
    setRandomSpin();
    void
    setFieldh(double h);
    void
    setTemp(double temp);
    void
    setMethod(char method);
    void
    Move();
    void
    Measure();
    void
    Equilibrate(unsigned int lowstep = 1000);
    void
    Simulate();
    string
    getAvailableVar();
    double
    getMeasure(char var);
    void
    blockingMethod();
};

vector<double>
getTempRange(double ti, double tf, unsigned int size);


#endif // ifndef __I1D__

/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/
