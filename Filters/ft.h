
#ifndef FT_H
#define FT_H

#include "signalgen.h"

class FT
{
public:
    FT();

    static std::vector<std::complex<double>> dft(const std::vector<std::complex<double>> &signal);
    static std::vector<std::complex<double>> forwardDft(const std::vector<std::complex<double>> &signal);

    static std::vector<std::complex<double>> fft(std::vector<std::complex<double> > a, bool invert);

    template<class T>
    static void fft_impl(std::vector<std::complex<T>> &out, size_t len, bool forward);
    static int fftWithoutWindow(int len, std::complex<double> *in, std::complex<double> *out, bool forward);

};

#endif // FT_H
