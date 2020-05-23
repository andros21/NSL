#include "func.h"

static string outPath = "../out/";

void
randomOnCircle(unsigned int npoints, double radius, string ofile)
{
    Random rnd;

    rnd.SetRandom();

    double theta;
    ofstream Cities(outPath + ofile + ".tsv");

    if (Cities.is_open()) {
        for (unsigned int i = 0; i < npoints; ++i) {
            theta = rnd.Rannyu(0., 2 * M_PI);
            Cities << radius * cos(theta) << "\t" << radius * sin(theta) << endl;
        }
    } else { cerr << "ERROR: Unable to open " << ofile << ".tsv" << endl; }
    Cities.close();
}

void
randomInSquare(unsigned int npoints, double edge, string ofile)
{
    Random rnd;

    rnd.SetRandom();

    ofstream Cities(outPath + ofile + ".tsv");

    if (Cities.is_open()) {
        for (unsigned int i = 0; i < npoints; ++i) {
            Cities << rnd.Rannyu(0., edge) << "\t" << rnd.Rannyu(0., edge) << endl;
        }
    } else { cerr << "ERROR: Unable to open " << ofile << ".tsv" << endl; }
    Cities.close();
}
