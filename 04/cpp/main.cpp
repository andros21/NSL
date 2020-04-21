/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/

#include "../../MD/moldyn.h"
#include "func.h"

using namespace std;

int
main(int argc, char * argv[])
{
    /***************************
    *  Exercise 04.1
    ***************************/

    // MOLECOLAR DYNAMIC {SOLID, LIQUID, GAS}
    MolDyn mds("solid");
    MolDyn mdl("liquid");
    MolDyn mdg("gas");

    vector<MolDyn> vmd = { mds, mdl, mdg };

    for (auto & md : vmd) {
        md.Simulate(1000); // Initial simulation 1000 steps, random velocities
        md.End();          // Write final configuration to files
        md.Restart();      // Restart using these configuration files,
        ;                  // reaching equilibrium with temperature
        ;                  // to avoid initial oscillations
    }

    /***************************
    *  Exercise 04.2, 04.3
    ***************************/

    for (auto & md : vmd) {
        md.blockingMethod(100); // BlockingMethod with nblk 100
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
