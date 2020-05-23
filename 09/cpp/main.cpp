/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/

#include "../../TSP/tsp.h"
#include "func.h"

int
main(int argc, char * argv[])
{
    unsigned int nIter(10);

    // randomOnCircle(32);
    // randomInSquare(32);

    Population popOC("circle", 50);
    Population popIS("square", 50);

    /***************************
    *  Exercise 09.1.1
    ***************************/

    GeneticAlgo GAOC(0.1, 0.5);

    vector<double> bestL1OC(nIter);
    vector<double> bestHalfAvgL1OC(nIter, 0.0d);

    for (unsigned int i = 0; i < nIter; ++i) {
        popOC       = GAOC.evolvePopulation(popOC);
        bestL1OC[i] = popOC.getTour(0).getDistance();
        // for (unsigned int j = 0; j < (unsigned int) (0.5 * popOC.getPopulationSize()); ++j)
        //     bestHalfAvgL1OC[i] += popOC.getTour(j).getDistance();
        // bestHalfAvgL1OC[i] /= (0.5 * popOC.getPopulationSize());
    }

    popOC.getTour(0).writeTourToFile("OC");
    writeVector(bestL1OC, "bL1OC");
    // writeVector(bestHalfAvgL1OC, "bHAL1OC");

    /***************************
    *  Exercise 09.1.2
    ***************************/

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
