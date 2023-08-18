
#ifndef SIGNALGEN_H
#define SIGNALGEN_H

#include <vector>
#include <complex>

class SignalGen
{
public:
    SignalGen();

    static std::vector<std::complex<double>> getSin(int samplesCount, double ampl, double f, double sampleRate);
    static std::vector<std::complex<double>> addSomeNoise(const std::vector<std::complex<double>> &signal, double freq, double ampl, int noiseCount);
};

#endif // SIGNALGEN_H
