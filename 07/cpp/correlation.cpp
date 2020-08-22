/****************************************************************
 *****************************************************************
 *  _/    _/  _/_/_/  _/      Numerical Simulation Laboratory
 * _/_/  _/ _/       _/       Physics Department
 * _/  _/_/    _/    _/       Universita' degli Studi di Milano
 * _/    _/       _/ _/       Prof. D.E. Galli
 * _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/

#include <iostream>
#include <iterator>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

// read vector from file
vector<double>
readVector(string iname)
{
    ifstream is("../out/" + iname);
    istream_iterator<double> start(is), end;
    vector<double> vct(start, end);

    return vct;
}

// write the input vector to a file inside out dir
void
writeVector(vector<double> & vct, string ofile)
{
    ofstream ofs("../out/" + ofile + ".out");

    ostream_iterator<double> output_iterator(ofs, "\n");

    if (ofs.is_open()) {
        copy(vct.begin(), vct.end(), output_iterator);
    } else { cerr << "ERROR: Unable to open " + ofile + ".out" << endl; }

    ofs.close();
}

// calculate auto-corr for the given
// thermodynamic var and state,
// then write vct to file
void
acorrelation(string cth, string cstate)
{
    vector<double> TT(readVector(cth + "-" + cstate + ".out"));
    vector<double> CC(TT.size());

    unsigned int tmax = TT.size() - 1;

    for (unsigned int i = 0; i <= (tmax - 1); ++i) {
        cout << "\rStep: " << i << flush;
        vector<double> sums(5, 0.0d);
        for (unsigned int j = 0; j < tmax; ++j) {
            if (j < (tmax - i)) {
                sums[0] += TT[j] * TT[j + i];
                sums[1] += TT[j];
                sums[2] += TT[j + i];
            }
            sums[3] += pow(TT[j], 2);
            sums[4] += TT[j];
        }
        CC[i] = ( (1. / (tmax - i)) * sums[0]
          - ( (1. / (tmax - i)) * sums[1] * (1. / (tmax - i)) * sums[2] ) )
          / ( (1. / tmax) * sums[3] - pow((1. / tmax) * sums[4], 2) );
    }

    writeVector(CC, "ac-" + cth + "-" + cstate + "-head");
}

int
main()
{
    vector<string> vars   = { "P", "V" };
    vector<string> states = { "solid", "liquid", "gas" };

    for (auto state : states) {
        cout << state << endl;
        for (auto var : vars) {
            cout << var << endl;
            acorrelation(var, state);
            cout << endl;
        }
    }
    return 0;
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
