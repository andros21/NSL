/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/

#include "tsp.h"

static string outPath = "../out/";

ostream
&
operator << (ostream& os, const City& city)
{
    os << city.getXY().first << "," << city.getXY().second;
    return os;
}

bool
operator == (const City& city1, const City& city2)
{
    return city1.getXY() == city2.getXY();
}

Tour :: Tour(string ifile) : _tour(), _fitness(0.), _distance(0.)
{
    ifstream Cities(outPath + ifile + ".tsv");
    pair<double, double> parser;

    if (Cities.is_open()) {
        Cities >> parser.first >> parser.second;
        while (!Cities.eof()) {
            _tour.push_back(City(parser));
            Cities >> parser.first >> parser.second;
        }
    } else { cerr << "ERROR: Unable to open " << ifile << ".tsv" << endl; }
    Cities.close();
}

void
Tour :: writeTourToFile(string ofile) const
{
    ofstream ofs(outPath + ofile + ".csv");

    if (ofs.is_open()) {
        for (auto el : _tour)
            ofs << el << endl;
    } else { cerr << "ERROR: Unable to open " + ofile + ".out" << endl; }

    ofs.close();
}

double
Tour :: getFitness()
{
    if (_fitness == 0.)
        _fitness = 1. / getDistance();
    return _fitness;
}

double
Tour :: getDistance()
{
    if (_distance == 0.) {
        double distance(0.);
        for (auto it = _tour.begin(); it != _tour.end() - 1; it++)
            distance += it->distanceL1(*(it + 1));
        distance += _tour.end()->distanceL1(*(_tour.begin()));
        _distance = distance;
    }
    return _distance;
}

void
Tour :: mutate(Random & rnd)
{
    for (auto it = ( _tour.begin() + 1); it != _tour.end(); ++it)
        iter_swap(it, _tour.begin() + rnd.RannyuDiscrete(1, _tour.size() - 1));

    this->resetFit();
}

Tour
Tour :: crossover(const Tour& tour, Random & rnd, double crossoverRate)
{
    if (_tour.size() != tour.numberOfCities() )
        cerr << "ERR: For crossing over tour size must be equal" << endl;

    if (rnd.Rannyu() < crossoverRate) {
        Tour child(_tour.size());

        unsigned int startPos = rnd.RannyuDiscrete(1, _tour.size() - 1);
        unsigned int endPos   = rnd.RannyuDiscrete(1, _tour.size() - 1);

        while (startPos == endPos) {
            endPos = rnd.RannyuDiscrete(1, _tour.size() - 1);
        }

        if (startPos > endPos) {
            unsigned int appoPos = startPos;
            startPos = endPos;
            endPos   = appoPos;
        }

        Tour sub_tour(endPos + 1 - startPos);

        for (unsigned int i = startPos; i <= endPos; ++i) {
            child.setCity(i, _tour[i]);
            sub_tour.setCity(i - startPos, _tour[i]);
        }

        unsigned int j = 0;

        for (unsigned int i = 0; i < tour.numberOfCities(); ++i) {
            if (i < startPos || i > endPos) {
                while (sub_tour.containCity(tour.getCity(j))) ++j;
                child.setCity(i, tour.getCity(j));
                ++j;
            }
        }
        // cout << startPos << "," << endPos << "|" << sub_tour.numberOfCities() << "|" << endl;

        return child;
    } else { return tour; }
} // Tour::crossover

Tour
shuffleTour(const Tour& tour)
{
    Tour Ctour(tour);

    Ctour.shuffleTour();

    return Ctour;
}

Population :: Population(string initialTourFile, unsigned int populationSize)
{
    Tour tour(initialTourFile);

    for (unsigned int i = 0; i < populationSize; ++i)
        _tours.insert(shuffleTour(tour));
}

Population :: Population (const Tour& initialTour, unsigned int populationSize)
{
    for (unsigned int i = 0; i < populationSize; ++i)
        _tours.insert(shuffleTour(initialTour));
}

void
Population :: mutate(Random & rnd, double mutationRate)
{
    for (auto it = (next(_tours.begin(), 1)); it != _tours.end(); ++it) {
        if (rnd.Rannyu() < mutationRate) {
            _tours.erase(it);
            Tour tour(this->getTour(it));
            do
                tour.mutate(rnd);
            while (!_tours.insert(tour).second);
        }
    }
}

Tour
GeneticAlgo :: tourSelection(const Population& pop)
{
    return pop.getTour(_rnd.RannyuDiscrete(1, pop.getPopulationSize() - 1));
}

Population
GeneticAlgo :: evolvePopulation(const Population& pop)
{
    Population newpop;

    newpop.addTour(pop.getTour(0));

    while (newpop.getPopulationSize() < pop.getPopulationSize()) {
        Tour parent1(this->tourSelection(pop));
        Tour parent2(this->tourSelection(pop));
        Tour child(parent1.crossover(parent2, _rnd, _crossoverRate));
        newpop.addTour(child);
    }


    // cout << newpop.getTour(0).getCity(0) << endl;
    // newpop.mutate(_rnd, _mutationRate);

    return newpop;
}

/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/
