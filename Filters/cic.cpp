
#include "cic.h"

CIC::CIC()
{

}

std::vector<std::complex<double>> CIC::filter(const std::vector<std::complex<double>> &signal, int decimation, int stages)
{
    double gain = pow(decimation * 1, stages);

    int c_stages = stages;
    int i_stages = stages;


    std::vector<Integrator> intes(i_stages);
    std::vector<Comb> combs(c_stages);

    std::vector<std::complex<double>> filtered;

    for (int s = 0, v = 0; s < signal.size() ; s++)
    {
        std::complex<double> z = signal[v];
        for (int i = 0; i < i_stages; i++)
        {
            z = intes[i].update(z);
        }
        if (s % decimation == 0)
        {
            std::complex<double> j;
            for (int c = 0; c < c_stages; c++)
            {
                z = combs[c].update(z);
                j = z;
            }
            filtered.push_back(j / gain);
        }
    }
    return filtered;
}
