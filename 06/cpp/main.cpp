/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/

#include "../../I1D/i1d.h"
#include "func.h"

using namespace std;

int
main(int argc, char * argv[])
{
    /***************************
    *  Exercise 06.1.1
    ***************************/

    I1D isi;
    double tf(0.5);
    unsigned int t_size(10);

    /****
    *  after first run, uncomment method below
    *  and re-run to calculate the same thermodynamic variables with h!=0
    ****/
    // isi.setFieldh(0.02);

    for (auto t : getTempRange(isi.getTemp(), tf, t_size)) {
        isi.setTemp(t);
        isi.Equilibrate(1000);
        isi.blockingMethod();
    }

    /***************************
    *  Exercise 06.1.2
    ***************************/

    isi.setTemp(2.);
    isi.setRandomSpin();
    isi.setMethod('G');

    for (auto t : getTempRange(isi.getTemp(), tf, t_size)) {
        isi.setTemp(t);
        isi.Equilibrate(1000);
        isi.blockingMethod();
    }

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
