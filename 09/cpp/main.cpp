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

    // randomOnCircle(32);
    //
    // Population popOC("circle", 50);
    //
    // GeneticAlgo GAOC(0.5, 0.4);
    //
    // GAOC.evolvePopulation(popOC, "OC", 300, true, true, true);

    /***************************
    *  Exercise 09.1.2
    ***************************/

    // randomInSquare(32);

    Population popIS("square", 50);

    GeneticAlgo GAIS(0.8, 0.5);

    GAIS.evolvePopulation(popIS, 700, "IS", true, true, true);

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
