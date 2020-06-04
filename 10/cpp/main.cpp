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
    *  Exercise 10.1.1
    ***************************/

    double temp        = 1.1;
    unsigned int nIter = 1000;
    unsigned int nChan = 100;

    randomOnCircle(32);

    Tour tourOC("circle");
    tourOC.measureFit();

    vector<double> bestL1OC(nIter);

    SimulatedAnnealingAlgo SAAOC(temp);

    for (unsigned int i = 1; i <= nIter; ++i) {
        tourOC = SAAOC.evolveTour(tourOC);
        // SAAOC.coutMetroRatio();
        bestL1OC[i - 1] = tourOC.getDistance();
        if (i % nChan == 0) {
            temp = temp - 0.1;
            SAAOC.setTemp(temp);
        }
    }

    tourOC.writeTourToFile("OC");
    writeVector(bestL1OC, "bL1OC");

    /***************************
    *  Exercise 10.1.2
    ***************************/

    temp  = 1.1;
    nIter = 1000;
    nChan = 100;

    randomInSquare(32);

    Tour tourIS("square");
    tourIS.measureFit();

    vector<double> bestL1IS(nIter);

    SimulatedAnnealingAlgo SAAIS(temp);

    for (unsigned int i = 1; i <= nIter; ++i) {
        tourIS = SAAIS.evolveTour(tourIS);
        // SAAOC.coutMetroRatio();
        bestL1IS[i - 1] = tourIS.getDistance();
        if (i % nChan == 0) {
            temp = temp - 0.1;
            SAAIS.setTemp(temp);
        }
    }

    tourIS.writeTourToFile("IS");
    writeVector(bestL1IS, "bL1IS");


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
