#include "func.h"

static string outPath = "../out/";

void
randomOnCircle(unsigned int npoints, double radius, string ofile)
{
    Random rnd;

    rnd.SetRandom();

    double theta;
    pair<double, double> spt;
    ofstream Cities(outPath + ofile + ".tsv");

    if (Cities.is_open()) {
        for (unsigned int i = 0; i < npoints; ++i) {
            theta = rnd.Rannyu(0., 2 * M_PI);
            if (i == 0) spt = { radius * cos(theta), radius * sin(theta) };
            Cities << radius * cos(theta) << "\t" << radius * sin(theta) << endl;
        }
        Cities << spt.first << "\t" << spt.second << endl;
    } else { cerr << "ERROR: Unable to open " << ofile << ".tsv" << endl; }
    Cities.close();
}

void
randomInSquare(unsigned int npoints, double edge, string ofile)
{
    Random rnd;

    rnd.SetRandom();

    pair<double, double> xy;
    pair<double, double> spt;
    ofstream Cities(outPath + ofile + ".tsv");

    if (Cities.is_open()) {
        for (unsigned int i = 0; i < npoints; ++i) {
            xy = { rnd.Rannyu(0., edge), rnd.Rannyu(0., edge) };
            if (i == 0) spt = { xy.first, xy.second };
            Cities << xy.first << "\t" << xy.second << endl;
        }
        Cities << spt.first << "\t" << spt.second << endl;
    } else { cerr << "ERROR: Unable to open " << ofile << ".tsv" << endl; }
    Cities.close();
}
