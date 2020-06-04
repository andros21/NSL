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
Tour :: coutTour() const
{
    for (auto el : _tour)
        cout << el << "|";
    cout << endl;
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

void
Tour :: measureFit()
{
    if (_distance == 0.) {
        double distance(0.);
        for (auto it = _tour.begin(); it != _tour.end() - 1; it++)
            distance += it->distanceL1(*(it + 1));
        _distance = distance;
    }
    _fitness = 1. / _distance;
}

void
Tour :: mutate(char type, unsigned int n, unsigned int m)
{
    if (type == '1')
        iter_swap(_tour.begin() + n, _tour.begin() + m);
    if (type == '2')
        reverse(_tour.begin() + n, _tour.begin() + m);

    this->resetFit();
}

Tour
Tour :: crossover(const Tour& tour, unsigned int n, unsigned int m)
{
    if (_tour.size() != tour.numberOfCities() )
        cerr << "ERR: For crossing over tour size must be equal" << endl;

    if (n > m) {
        unsigned int appoPos = n;
        n = m;
        m = appoPos;
    }

    Tour child(_tour.size());
    Tour sub_tour(m + 1 - n);

    for (unsigned int i = n; i <= m; ++i) {
        child.setCity(i, _tour[i]);
        sub_tour.setCity(i - n, _tour[i]);
    }

    unsigned int j = 0;

    for (unsigned int i = 0; i < tour.numberOfCities(); ++i) {
        if (i < n || i > m) {
            while (sub_tour.containCity(tour.getCity(j))) ++j;
            child.setCity(i, tour.getCity(j));
            ++j;
        }
    }

    return child;
} // Tour::crossover

Population :: Population(string initialTourFile, unsigned int populationSize)
{
    Tour tour(initialTourFile);
    Random rnd;

    rnd.SetRandom();

    auto nc = tour.numberOfCities();
    for (unsigned int i = 0; i < populationSize; ++i) {
        tour.mutate('1', rnd.RannyuDiscrete(1, nc - 2), rnd.RannyuDiscrete(1, nc - 2));
        _tours.push_back(tour);
    }

    this->sortPopulation();

    // for (auto & tour : _tours) {
    //     auto FF = tour.getFitness();
    //     cout << FF << "| " << count_if(_tours.begin(), _tours.end(), [&FF](Tour & TT){
    //         return TT.getFitness() == FF;
    //     }) << endl;
    // }
}

void
Population :: sortPopulation()
{
    for (auto & tour : _tours)
        tour.measureFit();
    sort(_tours.begin(), _tours.end(), TourCompare() );
}

Tour
GeneticAlgo :: tourSelection(const Population& pop)
{
    auto rr = _rnd.Rannyu();

    for (unsigned int i = 1; i < pop.getPopulationSize(); ++i)
        if (pop.getTour(i).getFitness() < rr)
            return pop.getTour(i);

    return pop.getTour(1);
}

Population
GeneticAlgo :: evolvePopulation(const Population& pop)
{
    Population newpop;

    newpop.addTour(pop.getTour(0));
    auto nc = newpop.getTour(0).numberOfCities();

    while (newpop.getPopulationSize() < pop.getPopulationSize()) {
        Tour parent1(this->tourSelection(pop));
        Tour parent2(this->tourSelection(pop));
        if (_rnd.Rannyu() < _crossoverRate) {
            auto nn = _rnd.RannyuDiscrete(1, nc - 2);
            auto mm = _rnd.RannyuDiscrete(1, nc - 2);
            Tour child1(parent1.crossover(parent2, nn, mm));
            Tour child2(parent2.crossover(parent1, nn, mm));
            newpop.addTour(child1);
            newpop.addTour(child2);
        } else {
            newpop.addTour(parent1);
            newpop.addTour(parent2);
        }
    }

    for (unsigned int i = 1; i < newpop.getPopulationSize(); ++i) {
        Tour newtour(newpop.getTour(i));
        if (_rnd.Rannyu() < _mutationRate) {
            newtour.mutate('1', _rnd.RannyuDiscrete(1, nc - 2), _rnd.RannyuDiscrete(1, nc - 2));
            newpop.setTour(newtour, i);
        }
        if (_rnd.Rannyu() < _mutationRate) {
            newtour.mutate('2', _rnd.RannyuDiscrete(1, nc - 2), _rnd.RannyuDiscrete(1, nc - 2));
            newpop.setTour(newtour, i);
        }
    }

    newpop.sortPopulation();

    // cout << newpop.getTour(0).getFitness() << "|" << newpop.getTour(1).getFitness() << "|"
    //      << newpop.getTour(2).getFitness() << endl;
    // newpop.getTour(0).coutTour();


    // for (unsigned int i = 0; i < newpop.getPopulationSize(); ++i)
    //     cout << newpop.getTour(i).getDistance() << "|";
    // cout << endl;

    return newpop;
} // GeneticAlgo::evolvePopulation

double
SimulatedAnnealingAlgo :: getMetroRatio()
{
    if (Mttot == 0) {
        return 0.;
    }
    return (double) Macpt / Mttot;
}

void
SimulatedAnnealingAlgo :: coutMetroRatio()
{
    cout << "\rMetropolis Ratio: " << round(this->getMetroRatio() * 100) << " " << flush;
}

void
SimulatedAnnealingAlgo :: resetMetroRatio()
{
    Macpt = 0;
    Mttot = 0;
}

char
SimulatedAnnealingAlgo :: tourSelection(double lp, double lc)
{
    double alpha, rr;

    alpha = min(1., exp(-(1. / _temp) * (lc - lp)));
    rr    = _rnd.Rannyu();

    Mttot++;
    if (rr <= alpha) {
        Macpt++;
        return 'c';
    } else {
        return 'p';
    }
}

Tour
SimulatedAnnealingAlgo :: evolveTour(const Tour& tour)
{
    auto nc = tour.numberOfCities();
    Tour parent(tour);
    Tour child(parent);

    if (_rnd.Rannyu() <= 0.5)
        child.mutate('1', _rnd.RannyuDiscrete(1, nc - 2), _rnd.RannyuDiscrete(1, nc - 2));
    else
        child.mutate('2', _rnd.RannyuDiscrete(1, nc - 2), _rnd.RannyuDiscrete(1, nc - 2));
    parent.measureFit();
    child.measureFit();
    if (tourSelection(parent.getDistance(), child.getDistance()) == 'c')
        parent = child;
    // cout << Macpt << "|" << Mttot << endl;

    return parent;
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
