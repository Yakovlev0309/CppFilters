
#ifndef CIC_H
#define CIC_H

#include <vector>
#include <complex>

class CIC
{
public:
    CIC();

    static std::vector<std::complex<double>> moving_average_filter(const std::vector<std::complex<double>>& input, int window_size);
    static std::vector<std::complex<double>> cic_decimation_filter(const std::vector<std::complex<double>>& input, int decimation_factor, int num_stages, int differential_delay);
    static std::vector<std::complex<double>> cic_interpolation_filter(const std::vector<std::complex<double>>& input, int interpolation_factor, int num_stages, int differential_delay);
};

#endif // CIC_H
