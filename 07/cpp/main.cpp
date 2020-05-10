/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/

#include "../../MCMD/mcmd.h"
#include "func.h"

using namespace std;

int
main(int argc, char * argv[])
{
    /***************************
    *  Exercise 07.1
    ***************************/

    // MONTECARLO MOLECOLAR DYNAMIC {SOLID, LIQUID, GAS}
    McMd mds("solid");
    McMd mdl("liquid");
    McMd mdg("gas");

    vector<McMd> vmd = { mds, mdl, mdg };

    for (auto & md : vmd) {
        md.Simulate(5e4, false, true, true); // Initial simulation 1000 steps, random velocities
        md.resetMetroRatio();                // Reset Metropolis ratio
        md.End();                            // Write final configuration to files
        cout << endl;
        // md.Restart();      // Restart using these configuration files,
        ; // reaching equilibrium with temperature
        ; // to avoid initial oscillations
    }

    /***************************
    *  Exercise 07.4
    ***************************/

    // for (auto & md : vmd) {
    //     md.blockingMethod(100); // BlockingMethod with nblk 100
    // }

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
