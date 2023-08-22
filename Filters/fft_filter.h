
#ifndef FFT_FILTER_H
#define FFT_FILTER_H

#include "signalgen.h"


enum FILTRATION_TYPE {
    low   = 0,
    high  = 1,
    bpass = 2,
    none  = 3
};

class FFT_Filter
{
public:
    FFT_Filter();

    static int filtration(FILTRATION_TYPE ftype, int len, std::complex<double> *in, std::complex<double> *out, double low, double high, double samplingRate, int windowSize);
    static int fftWithoutWindow(int len, std::complex<double> *in, std::complex<double> *out, bool forward);

    static std::vector<std::complex<double>> filter(const std::vector<std::complex<double>> &signal, double cutoffFreq, double sampleRate);
};

#endif // FFT_FILTER_H
