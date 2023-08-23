
#ifndef FIR_H
#define FIR_H

#include "signalgen.h"


class FIR
{
public:
    FIR();
    
    static std::vector<std::complex<double>> compensatePhaseDelay(const std::vector<std::complex<double>> &signal, int filterSize);

    static std::vector<double> getLowPassFilterCoeffs(int filterSize, double cutoffFreq, double sampleRate);
    static std::vector<double> getHighPassFilterCoeffs(int filterSize, double cutoffFreq, double sampleRate);

    static std::vector<double> calculateKernels(int filterOrder, double cutoffFreq, double samplingFreq);

    static std::vector<double> calculate_fir_filter_coefficients(int N, double cutoff_freq);
    static std::vector<std::complex<double>> apply_fir_filter(const std::vector<std::complex<double>>& x, const std::vector<double>& h);

    static std::vector<double> getFilterCoeffs(int filterSize, double passbandFreq, double stopbandFreq, double sampleRate);
    static std::vector<double> getFilterCoeffs(int filterSize, double cutoffFreq, double sampleRate);
    static std::vector<std::complex<double>> filter(const std::vector<std::complex<double>> &signal, const std::vector<double> &coeffs);
    template <class T>
    static std::vector<T> slidingWindowAverage(const std::vector<T> &array, int win);
};

#endif // FIR_H
