#ifndef FFT_H
#define FFT_H

#include <math.h>
#include <complex>
#include <functional>
#include <vector>

namespace Fft
{
    template <class T>
    std::vector<std::complex<T> > fft(const std::vector<T>& data, bool forward = true);

    template <class T>
    std::vector<std::complex<T> > fft(const std::vector<std::complex<T> >& data, bool forward = true);

    template <class T>
    void fft(const std::vector<T>& data, std::vector<std::complex<T> >&out, size_t len, bool forward = true);

    template <class T>
    void fft(const std::vector<std::complex<T> >& data, std::vector<std::complex<T> >&out, size_t len, bool forward = true);
}


namespace Dft
{
    void dftNear(const std::vector<std::complex<double> > &signal, std::vector<double>& out, double peakCarrier, int halfWin, double carrierStep, double sampleRate, std::function<double (int)> &func);
}


#endif // FFT_H
