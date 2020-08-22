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

    I1D isi;                 // Load class, using default algo: Metropolis
    double tf(0.5);          // Final temperature
    unsigned int t_size(10); // Temperature array size

    /****
    *  after first run, uncomment method below
    *  and re-run to calculate the same thermodynamic variables with h!=0
    ****/
    // isi.setFieldh(0.02);

    for (auto t : getTempRange(isi.getTemp(), tf, t_size)) {
        isi.setTemp(t);         // Set temperature every step
        isi.Equilibrate(50000); // Equilibrate before real simulation
        isi.blockingMethod();   // Blocking method at current temp
    }

    /***************************
    *  Exercise 06.1.2
    ***************************/

    isi.setTemp(2.);     // Reset initial temperature
    isi.setRandomSpin(); // Reset spins
    isi.setMethod('G');  // Set Gibbs algo for movements

    for (auto t : getTempRange(isi.getTemp(), tf, t_size)) {
        isi.setTemp(t);         // Set temperature every step
        isi.Equilibrate(50000); // Equilibrate before real simulation
        isi.blockingMethod();   // Blocking method at current temperature
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
