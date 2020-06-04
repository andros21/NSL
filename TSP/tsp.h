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

    void
    coutTour() const;

    bool
    containCity(City city){ return find(_tour.begin(), _tour.end(), city) != _tour.end(); }

    void
    resetFit(){ _fitness = 0.; _distance = 0.; }

    void
    measureFit();

    unsigned int
    numberOfCities() const { return _tour.size(); }

    void
    writeTourToFile(string ofile) const;

    void
    generateIndividual();

    double
    getFitness() const { return _fitness; }

    double
    getDistance() const { return _distance; }

    void
    shuffleTour(Random & rnd);

    void
    mutate(char type, unsigned int n, unsigned int m);

    Tour
    crossover(const Tour& tour, unsigned int n, unsigned int m);
};

struct TourCompare {
    bool
    operator () (const Tour & lhs, const Tour & rhs)
    {
        return lhs.getFitness() > rhs.getFitness();
    }
};

#endif // ifndef __Tour__

#ifndef __Population__
# define __Population__

class Population {
private:
    vector<Tour> _tours;

public:

    Population() : _tours(){ };

    Population(string inialTourFile, unsigned int populationSize);

    Population(const Population& pop){ *this = pop; }

    void
    addTour(const Tour& tour){ _tours.push_back(tour); }

    Tour
    getTour(unsigned int index) const { return _tours[index]; }

    vector<Tour>
    getTours() const { return _tours; }

    void
    setTour(const Tour& tour, unsigned index){ _tours[index] = tour; }

    unsigned int
    getPopulationSize() const { return _tours.size(); }

    void
    sortPopulation();
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

    Population
    evolvePopulation(const Population& pop);
};

#endif // ifndef __GeneticAlgo__


#ifndef __GeneticAlgo__
# define __GeneticAlgo_

class SimulatedAnnealingAlgo {
private:
    Random _rnd;
    double _temp;
    unsigned int Macpt, Mttot;
public:
    SimulatedAnnealingAlgo(double temp) : _temp(temp), Macpt(0), Mttot(0){ _rnd.SetRandom(); }

    double
    getTemp() const { return _temp; }

    void
    setTemp(double temp){ _temp = temp; }

    double
    getMetroRatio();

    void
    coutMetroRatio();

    void
    resetMetroRatio();

    char
    tourSelection(double lp, double lc);

    Tour
    evolveTour(const Tour& tour);
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
