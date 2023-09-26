
#ifndef FIR_H
#define FIR_H

#include "signalgen.h"


class FIR
{
public:
    FIR();
    
    static int pow2(int len);
    static bool isntPow2(int len);

    static std::vector<std::complex<double>> compensatePhaseDelay(const std::vector<std::complex<double>> &signal, int filterSize);
    static std::vector<std::complex<double>> compensateDelay(const std::vector<std::complex<double>> &signal, int filterSize, int oldLen);

    static std::vector<double> getLowPassFilterCoeffs(int filterSize, double cutoffFreq, double sampleRate);
    static std::vector<double> getHighPassFilterCoeffs(int filterSize, double cutoffFreq, double sampleRate);
    static std::vector<double> getBandPassFilterCoeffs(int filterSize, double lowCutoffFreq, double highCutoffFreq, double sampleRate);

    static std::vector<double> calculateKernels(int filterOrder, double cutoffFreq, double samplingFreq);

    static std::vector<double> calculate_fir_filter_coefficients(int N, double cutoff_freq);
    static std::vector<std::complex<double>> applyFilter(const std::vector<std::complex<double>>& signal, const std::vector<double>& coeffs);
    static std::vector<std::complex<double>> applyFilterInSpectrum(const std::vector<std::complex<double>>& signal, const std::vector<double>& coeffs);

    static std::vector<double> getFilterCoeffs(int filterSize, double passbandFreq, double stopbandFreq, double sampleRate);
    static std::vector<double> getFilterCoeffs(int filterSize, double cutoffFreq, double sampleRate);
    static std::vector<std::complex<double>> filter(const std::vector<std::complex<double>> &signal, const std::vector<double> &coeffs);
    template <class T>
    static std::vector<T> slidingWindowAverage(const std::vector<T> &array, int win);
};

#endif // FIR_H
