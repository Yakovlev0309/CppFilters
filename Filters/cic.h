
#ifndef CIC_H
#define CIC_H

#include "signalgen.h"

struct Integrator
{
    Integrator() : yn(0.0), ynm(0.0) {}
    std::complex<double> update(std::complex<double> inp)
    {
        ynm = yn;
        yn = ynm + inp;
        return yn;
    }

    std::complex<double> yn;
    std::complex<double> ynm;
};

struct Comb
{
    Comb() : xn(0.0), xnm(0.0) {}
    std::complex<double> update(std::complex<double> inp)
    {
        xnm = xn;
        xn = inp;
        return xn - xnm;
    }

    std::complex<double> xn;
    std::complex<double> xnm;
};

class CIC
{
public:
    CIC();

    static std::vector<std::complex<double>> filter(const std::vector<std::complex<double>> &signal, int decimation, int stages);
};

#endif // CIC_H
