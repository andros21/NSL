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
    /***************************
    *  Exercise 09.1.1
    ***************************/

    unsigned int nIter = 45;

    randomOnCircle(32);

    Population popOC("circle", 1000);
    vector<double> bestL1OC(nIter);
    vector<double> bestHalfAvgL1OC(nIter, 0.0d);

    GeneticAlgo GAOC(0.2, 0.6);

    for (unsigned int i = 0; i < nIter; ++i) {
        popOC       = GAOC.evolvePopulation(popOC);
        bestL1OC[i] = popOC.getTour(0).getDistance();
        for (unsigned int j = 0; j < (unsigned int) (0.5 * popOC.getPopulationSize()); ++j)
            bestHalfAvgL1OC[i] += popOC.getTour(j).getDistance();
        bestHalfAvgL1OC[i] /= (0.5 * popOC.getPopulationSize());
    }

    popOC.getTour(0).writeTourToFile("OC");
    writeVector(bestL1OC, "bL1OC");
    writeVector(bestHalfAvgL1OC, "bHAL1OC");

    /***************************
    *  Exercise 09.1.2
    ***************************/

    randomInSquare(32);

    Population popIS("square", 1000);
    vector<double> bestL1IS(nIter);
    vector<double> bestHalfAvgL1IS(nIter, 0.0d);

    GeneticAlgo GAIS(0.2, 0.6);

    for (unsigned int i = 0; i < nIter; ++i) {
        popIS       = GAIS.evolvePopulation(popIS);
        bestL1IS[i] = popIS.getTour(0).getDistance();
        for (unsigned int j = 0; j < (unsigned int) (0.5 * popIS.getPopulationSize()); ++j)
            bestHalfAvgL1IS[i] += popIS.getTour(j).getDistance();
        bestHalfAvgL1IS[i] /= (0.5 * popIS.getPopulationSize());
    }

    popIS.getTour(0).writeTourToFile("IS");
    writeVector(bestL1IS, "bL1IS");
    writeVector(bestHalfAvgL1IS, "bHAL1IS");


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
