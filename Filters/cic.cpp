
#include "cic.h"

CIC::CIC()
{

}

std::vector<std::complex<double>> CIC::moving_average_filter(const std::vector<std::complex<double>>& input, int window_size)
{
    std::vector<std::complex<double>> output(input.size() - window_size + 1);
    std::complex<double> sum;
    for (int i = 0; i < output.size(); ++i) {
        sum = 0.0;
        for (int j = 0; j < window_size; ++j) {
            sum+= input[i + j];
        }
        output[i] = sum / double(window_size);
    }

    return output;
}

std::vector<std::complex<double>> CIC::cic_decimation_filter(const std::vector<std::complex<double>>& input, int decimation_factor, int num_stages, int differential_delay)
{
    std::vector<std::complex<double>> output(input.size() / decimation_factor);
    std::complex<double> sum;
    for (int i = 0; i < output.size(); ++i) {
        sum = input[i * decimation_factor];
        for (int j = 1; j < num_stages; ++j) {
            sum += input[i * decimation_factor - j * differential_delay];
        }
        output[i] = sum;
    }

    return output;
}

std::vector<std::complex<double>> CIC::cic_interpolation_filter(const std::vector<std::complex<double>>& input, int interpolation_factor, int num_stages, int differential_delay) {
    std::vector<std::complex<double>> output(input.size() * interpolation_factor);
    std::complex<double> sum;
    for (int i = 0; i < output.size(); ++i) {
        sum = input[i / interpolation_factor];
        for (int j = 1; j < num_stages; ++j) {
            sum += input[i / interpolation_factor - j * differential_delay];
        }
        output[i] = sum;
    }

    return output;
}
