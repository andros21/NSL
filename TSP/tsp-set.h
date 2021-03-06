/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/

using namespace std;

#include "../RC/random.h"

#ifndef __City__
# define __City__

class City {
private:
    pair<double, double> _pt;
public:
    City() : _pt(){ };
    City(double x, double y) : _pt(x, y){ };
    City(const pair<double, double> & pt) : _pt(pt){ };

    pair<double, double>
    getXY() const { return _pt; }

    double
    distanceL2(City city){ return pow(_pt.first - city.getXY().first, 2) + pow(_pt.second - city.getXY().second, 2); }

    double
    distanceL1(City city){ return sqrt(distanceL2(city)); }
};

ostream&
operator << (ostream& os, const City& city);

bool
operator == (const City& city1, const City& city2);

#endif // ifndef __City__

#ifndef __Tour__
# define __Tour__

class Tour {
private:
    vector<City> _tour;
    double _fitness, _distance;

public:

    Tour(string ifile);

    Tour(unsigned int size) : _tour(size), _fitness(0.), _distance(0.){ }

    Tour(const Tour& tour) : _fitness(0.), _distance(0.){ *this = tour; }

    void
    addCity(City city){ _tour.push_back(city); this->resetFit(); }

    void
    setCity(unsigned int index, City city){ _tour[index] = city; this->resetFit(); }

    City
    getCity(unsigned int index) const { return _tour[index]; }

    bool
    containCity(City city){ return find(_tour.begin(), _tour.end(), city) != _tour.end(); }

    void
    resetFit(){ _fitness = 0.; _distance = 0.; }

    unsigned int
    numberOfCities() const { return _tour.size(); }

    void
    writeTourToFile(string ofile) const;

    void
    generateIndividual();

    double
    getFitness();

    double
    getDistance();

    void
    shuffleTour(){ random_shuffle(_tour.begin() + 1, _tour.end()); this->resetFit(); }

    void
    mutate(Random & rnd);

    Tour
    crossover(const Tour& tour, Random & rnd, double crossoverRate);
};

Tour
shuffleTour(const Tour& tour);

struct TourCompare {
    bool
    operator () (Tour lhs, Tour rhs)
    {
        return lhs.getFitness() > rhs.getFitness();
    }
};

#endif // ifndef __Tour__

#ifndef __Population__
# define __Population__

class Population {
private:
    set<Tour, TourCompare> _tours;

public:

    Population() : _tours(){ };

    Population(string inialTourFile, unsigned int populationSize);

    Population(const Tour& initialTour, unsigned int populationSize);

    void
    addTour(const Tour& tour){ _tours.insert(tour); }

    Tour
    getTour(unsigned int index) const { auto it = next(_tours.begin(), index); return *it; }

    // Tour
    // getTour(set<Tour>::iterator it) const { return *it; }

    unsigned int
    getPopulationSize() const { return _tours.size(); }

    void
    mutate(Random & rnd, double mutationRate);
};

#endif // ifndef __Population__

#ifndef __GeneticAlgo__
# define __GeneticAlgo_

class GeneticAlgo {
private:
    Random _rnd;
    double _mutationRate, _crossoverRate;

public:

    GeneticAlgo(double mutationRate, double crossoverRate) : _mutationRate(mutationRate), _crossoverRate(crossoverRate)
    { _rnd.SetRandom(); }

    Tour
    tourSelection(const Population& pop);

    void
    evolvePopulation(const Population& pop, unsigned int nIter, string name, bool saveLastConf = true,
      bool writeBestToFile        = false,
      bool writeHalfBestAvgToFile = false);
};

#endif // ifndef __GeneticAlgo__

/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/
