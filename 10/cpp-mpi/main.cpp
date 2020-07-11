/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/

#include "../../TSP/tsp-mpi.h"
#include "func.h"

int
main(int argc, char * argv[])
{
    /***************************
    *  Exercise 10.2
    ***************************/

    unsigned int nIter = 50;
    unsigned int nMigr = 5;

    randomInSquare(32);

    Permutation permu;

    Population ipopIS("square", 1000);
    vector<Population> popIS = { ipopIS, ipopIS, ipopIS, ipopIS };

    vector<double> ibestL1IS(nIter);
    vector<vector<double> > bestL1IS = { ibestL1IS, ibestL1IS, ibestL1IS, ibestL1IS };

    vector<GeneticAlgo> GAIS = { GeneticAlgo(0.2, 0.6, "seed0.in"),
                                 GeneticAlgo(0.2, 0.6, "seed1.in"),
                                 GeneticAlgo(0.2, 0.6, "seed2.in"),
                                 GeneticAlgo(0.2, 0.6, "seed3.in") };

    pair<int, int> w1_pair;
    pair<int, int> w2_pair;

    mpi::environment env(argc, argv);
    mpi::communicator world;

    for (unsigned int i = 0; i < nIter; ++i) {
        if ((i + 1) % nMigr == 0) {
            if (world.rank() == 0) {
                permu.doPermutation();
                w1_pair = permu.get1Pair();
                w2_pair = permu.get2Pair();
                for (int j = 1; j < world.size(); ++j) {
                    world.send(j, 0, w1_pair);
                    world.send(j, 0, w2_pair);
                }
            } else {
                world.recv(0, 0, w1_pair);
                world.recv(0, 0, w2_pair);
            }
            if (world.rank() == w1_pair.first) {
                world.send(w1_pair.second, 0, popIS[w1_pair.first]);
                world.recv(w1_pair.second, 0, popIS[w1_pair.first]);
            } else if (world.rank() == w1_pair.second) {
                Population store1(popIS[w1_pair.second]);
                world.recv(w1_pair.first, 0, popIS[w1_pair.second]);
                world.send(w1_pair.first, 0, store1);
            } else if (world.rank() == w2_pair.first) {
                world.send(w2_pair.second, 0, popIS[w2_pair.first]);
                world.recv(w2_pair.second, 0, popIS[w2_pair.first]);
            } else if (world.rank() == w2_pair.second) {
                Population store2(popIS[w2_pair.second]);
                world.recv(w2_pair.first, 0, popIS[w2_pair.second]);
                world.send(w2_pair.first, 0, store2);
            } else { cout << "Error" << endl; }
        }
        for (int j = 0; j < world.size(); ++j) {
            if (world.rank() == j) {
                popIS[j]       = GAIS[j].evolvePopulation(popIS[j]);
                bestL1IS[j][i] = popIS[j].getTour(0).getDistance();
            }
        }
    }

    for (int j = 0; j < world.size(); ++j) {
        if (world.rank() == j) {
            popIS[j].getTour(0).writeTourToFile("IS-mpi-" + to_string(j));
            writeVector(bestL1IS[j], "bL1IS-mpi-" + to_string(j));
        }
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
