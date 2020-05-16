#include "func.h"

double
psi2(double x, double mu, double sigma)
{
    return (1.0 / 4.0) * M_SQRT2
           * (exp(2 * mu * x
           / pow(sigma,
           2)) + 2
           * exp(mu * x
           / pow(sigma,
           2)) + 1)
           * exp(-x * (mu + (1.0 / 2.0) * x)
             / pow(sigma, 2)) / (sqrt(M_PI) * sigma * (exp((1.0 / 2.0) * pow(mu, 2) / pow(sigma, 2)) + 1));
}

double
intg(double x, double mu, double sigma)
{
    return 1.0
           * (-0.125
           * pow(mu,
           2)
           * exp(mu * x
           / pow(sigma,
           2)) - 0.125
           * pow(mu,
           2) + 0.25 * mu * x
           * exp(mu * x
           / pow(sigma,
           2)) - 0.25 * mu * x + 1.0
           * pow(sigma,
           4)
           * pow(x,
           4)
           * exp(mu * x
           / pow(sigma,
           2)) + 1.0
           * pow(sigma,
           4)
           * pow(x,
           4) - 2.5
           * pow(sigma,
           4)
           * pow(x,
           2)
           * exp(mu * x
           / pow(sigma,
           2)) - 2.5
           * pow(sigma,
           4)
           * pow(x,
           2) + 0.25
           * pow(sigma,
           2)
           * exp(mu * x
           / pow(sigma,
           2)) + 0.25
           * pow(sigma,
           2) - 0.125
           * pow(x,
           2) * exp(mu * x / pow(sigma, 2)) - 0.125 * pow(x, 2)) / (pow(sigma, 4) * (exp(mu * x / pow(sigma, 2)) + 1));
} // intg

vector<double>
getRange(unsigned int NN, double min, double max)
{
    vector<double> vct(NN + 1);

    for (unsigned int i = 0; i <= NN; ++i)
        vct[i] = min + i * abs(max - min) / NN;

    return vct;
}
