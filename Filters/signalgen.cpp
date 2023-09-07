
#include "signalgen.h"
#include <QDebug>
#include <cmath>
#include <ctime>

SignalGen::SignalGen()
{

}

std::vector<std::complex<double>> SignalGen::getSin(int samplesCount, double ampl, double freq, double sampleRate)
{
    std::vector<std::complex<double>> signal;
    double real, imag, time;
    for(int t = 0; t < samplesCount; t++)
    {
        time = t / sampleRate;
        real = ampl * sin(2 * M_PI * freq * time);
        imag = ampl * cos(2 * M_PI * freq * time);
        signal.push_back({real, imag});
    }
    return signal;
}

std::vector<std::complex<double>> SignalGen::addSomeNoise(const std::vector<std::complex<double>> &signal, double ampl, double freq, int noiseCount, double freqFactor)
{
    if (noiseCount == 0)
        return signal;
    srand(time(NULL));
    std::vector<std::complex<double>> res = signal;
    double real, imag;
    for (int i = 0; i < noiseCount; i++)
    {
        for (int j = 0; j < signal.size(); j++)
        {
            real = (double)rand() / RAND_MAX * ampl * sin(2 * M_PI * (double)rand() / RAND_MAX * freq * freqFactor);
            imag = (double)rand() / RAND_MAX * ampl * cos(2 * M_PI * (double)rand() / RAND_MAX * freq * freqFactor);
            res[j] += std::complex<double>(real, imag);
        }
    }
    return res;
}
