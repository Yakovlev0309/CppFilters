
#include "fir.h"
#include "fft.hpp"
#include <cmath>

double blackmanWindow(int n, int windowSize)
{
    // При alpha = 0.16
    double a0 = 0.42, a1 = 0.5, a2 = 0.08;
    double angle = M_PI * n / (windowSize - 1);
    return a0 - a1 * cosl(2 * angle) + a2 * cosl(4 * angle);
}

double blackmanHarrisWindow(int n, int windowSize)
{
    double a0 = 0.35875, a1 = 0.48829, a2 = 0.14128, a3 = 0.01168;
    double angle = M_PI * n / (windowSize - 1);
    return a0 - a1 * cos(2 * angle) + a2 * cos(4 * angle) - a3 * cos(6 * angle);
}

double hammingWindow(int n, int windowSize)
{
    double a0 = 0.54, a1 = 0.46;
    double angle = M_PI * n / (windowSize - 1);
    return a0 - a1 * cos(2 * angle);
}

double nuttallWindow(int n, int windowSize)
{
    double a0 = 0.355768, a1 = 0.487396, a2 = 0.144232, a3 = 0.012604;
    double angle = M_PI * n / (windowSize - 1);
    return a0 - a1 * cos(2 * angle) + a2 * cos(4 * angle) - a3 * cos(6 * angle);
}

double blackmanNuttallWindow(int n, int windowSize)
{
    double a0 = 0.3635819, a1 = 0.4891775, a2 = 0.1365995, a3 = 0.0106411;
    double angle = M_PI * n / (windowSize - 1);
    return a0 - a1 * cos(2 * angle) + a2 * cos(4 * angle) - a3 * cos(6 * angle);
}

double flatTopWindow(int n, int windowSize)
{
    double a0 = 0.21557895, a1 = 0.41663158, a2 = 0.277263158, a3 = 0.083578947, a4 = 0.006947368;
    double angle = M_PI * n / (windowSize - 1);
    return a0 - a1 * cos(2 * angle) + a2 * cos(4 * angle) - a3 * cos(6 * angle) + a4 * cos(8 * angle);
}

FIR::FIR()
{

}

std::vector<std::complex<double>> FIR::compensateDelay(const std::vector<std::complex<double>> &signal, int filterSize, int oldLen)
{
    std::vector<std::complex<double>> withoutDelays = signal;
    auto begin = withoutDelays.begin();
    auto end = withoutDelays.begin() + filterSize / 2;
    withoutDelays.erase(begin, end);

    begin = withoutDelays.begin() + oldLen;
    end = withoutDelays.end();
    withoutDelays.erase(begin, end);

    return withoutDelays;
}

std::vector<double> FIR::getLowPassFilterCoeffs(int filterSize, double cutoffFreq, double sampleRate)
{
    std::vector<double> coefficients(filterSize);

    // Вычисляем центральный индекс
    double center = filterSize / 2.0;

    // Вычисляем значение сдвига частоты среза
    double normalizedCutoffFreq = cutoffFreq / sampleRate;

    double offset, coeff;
    // Вычисляем коэффициенты фильтра
    for (int n = 0; n < filterSize; n++)
    {
        // Вычисляем смещение относительно центра
        offset = n - center;

        if (offset == 0)
            coeff = 2 * normalizedCutoffFreq;
        else
            coeff = sin(2 * M_PI * normalizedCutoffFreq * offset) / (M_PI * offset);

        coeff *= blackmanWindow(n, filterSize);

        // Добавляем коэффициент в массив
        coefficients[n] = coeff;
    }

    return coefficients;
}

std::vector<double> FIR::getHighPassFilterCoeffs(int filterSize, double cutoffFreq, double sampleRate)
{
//    cutoffFreq += sampleRate / 2.0;

    std::vector<double> coefficients(filterSize);

    // Вычисляем центральный индекс
    double center = filterSize / 2.0;

    // Вычисляем значение сдвига частоты среза
    double normalizedCutoffFreq = cutoffFreq / sampleRate;

    double offset, coeff;
    // Вычисляем коэффициенты фильтра
    for (int n = 0; n < filterSize; n++)
    {
        // Вычисляем смещение относительно центра
        offset = n - center;

        if (offset == 0)
            coeff = 2 * normalizedCutoffFreq;
        else
            coeff = sin(2 * M_PI * normalizedCutoffFreq * offset) / (M_PI * offset);

        coeff *= pow(-1, n);

        coeff *= blackmanWindow(n, filterSize);

        // Добавляем коэффициент в массив
        coefficients[n] = coeff;
    }

    return coefficients;
}

std::vector<double> FIR::getBandPassFilterCoeffs(int filterSize, double lowCutoffFreq, double highCutoffFreq, double sampleRate)
{
    std::vector<double> coefficients(filterSize);

    // Вычисляем центральный индекс
    int center = filterSize / 2;

    // Вычисляем значение сдвига частот среза
    double lowNormalizedCutoffFreq = lowCutoffFreq / sampleRate;
    double highNormalizedCutoffFreq = highCutoffFreq / sampleRate;

    double offset, coeff;
    // Вычисляем коэффициенты фильтра
    for (int n = 0; n < filterSize; n++)
    {
        // Вычисляем смещение относительно центра
        offset = n - center;

        if (offset == 0)
        {
            coeff = 2 * (highNormalizedCutoffFreq - lowNormalizedCutoffFreq);
        }
        else
        {
            coeff = (sin(2 * M_PI * highNormalizedCutoffFreq * offset) - sin(2 * M_PI * lowNormalizedCutoffFreq * offset)) / (M_PI * offset);
        }

        coeff *= hammingWindow(n, filterSize);

        // Добавляем коэффициент в массив
        coefficients[n] = coeff;
    }

    return coefficients;
}

std::vector<double> FIR::calculateKernels(int filterOrder, double cutoffFreq, double samplingFreq)
{
    std::vector<double> kernels;

    // Calculate the normalized cutoff frequency
    double normalizedCutoff = cutoffFreq / samplingFreq;

    // Calculate the filter mid-point index
    int midPoint = filterOrder / 2;

    // Calculate the filter coefficients using the sinc function
    for (int n = 0; n < filterOrder; ++n)
    {
        if (n == midPoint)
        {
            kernels.push_back(2 * normalizedCutoff);
        }
        else
        {
            double kernel = std::sin(2 * M_PI * normalizedCutoff * (n - midPoint)) / (M_PI * (n - midPoint));
            kernels.push_back(kernel);
        }
    }

    return kernels;
}

std::vector<double> FIR::calculate_fir_filter_coefficients(int N, double cutoff_freq)
{
    std::vector<double> h(N);
    for (int n = 0; n < N; ++n) {
        if (n == (N-1)/2) {
            h[n] = 2 * M_PI * cutoff_freq;
        } else {
            h[n] = sin(2 * M_PI * cutoff_freq * (n - (N-1)/2)) / (M_PI * (n - (N-1)/2));
        }
    }
    return h;

//    std::vector<double> taps(N);
//    for (int n = 0; n < N; ++n) {
//        taps[n] = sin(2 * M_PI * cutoff_freq * n) / (M_PI * n) * 2 * cutoff_freq;
//    }
//    taps[0] = 2 * cutoff_freq;
//    return taps;
}

std::vector<std::complex<double>> FIR::applyFilter(const std::vector<std::complex<double>>& signal, const std::vector<double>& coeffs)
{
    int len = signal.size();
    int N = coeffs.size();
    int newLen = len + N - 1; // Новое количество отсчётов с учётом задержки
    std::vector<std::complex<double>> y(newLen);
    for (int n = 0; n < newLen; ++n) {
        y[n] = 0.0;
        for (int k = 0; k < N; ++k) {
            if (n - k >= 0 && n - k < len) {
                y[n] += coeffs[k] * signal[n - k];
            }
        }
    }
    return y;
}

int FIR::pow2(int len){

    char pow = 0;
    int num = 1;
    while(num < len && num*2 - num <= len-num)
    {
        num *= 2;
        pow++;
    }

    return pow;
}

bool FIR::isntPow2(int len){

    return (pow(2, pow2(len)) != len);
}

std::vector<std::complex<double>> FIR::applyFilterInSpectrum(const std::vector<std::complex<double>>& signal, const std::vector<double>& coeffs)
{
    int N = coeffs.size();
    int len = signal.size();
    int newLen = pow(2, ceil(log2(len + N - 1)));

    // Добавление нулей входному сигналу и фильтру
    std::vector<std::complex<double>> inWithZeroes(newLen);
    std::vector<std::complex<double>> coeffsWithZeroes(newLen);
    for (int i = 0; i < len; ++i)
    {
        inWithZeroes[i] = signal[i];
        if (i < N)
        {
            coeffsWithZeroes[i] = coeffs[i];
        }
    }
    auto inFft = Fft::fft(inWithZeroes, true);
    auto coeffsFft = Fft::fft(coeffsWithZeroes, true);

    // Умножение преобразованных сигналов
    std::vector<std::complex<double>> outputFft(newLen);
    for (int i = 0; i < newLen; ++i)
    {
        outputFft[i] = inFft[i] * coeffsFft[i];
    }
    return Fft::fft(outputFft, false);
}

std::vector<double> FIR::getFilterCoeffs(int filterSize, double passbandFreq, double stopbandFreq, double sampleRate)
{
    //    int newSizeDegree = (int)ceil(log((double)sizeIn)/log(2));
    //    sizeIn = pow(2,newSizeDegree);

    // sampleRate - частота дискретизации входных данных
    // passbandFreq - частота полосы пропускания
    // attenuationFreq - частота полосы затухания

    std::vector<double> H(filterSize);    // Импульсная характеристика фильтра
    std::vector<double> H_id(filterSize); // Идеальная импульсная характеристика
    std::vector<double> W(filterSize);    // Весовая функция

    // Расчет импульсной характеристики фильтра
    double Fc = (passbandFreq + stopbandFreq) / (2 * sampleRate);

    for (int i = 0; i < filterSize; i++)
    {
        if (i == 0)
        {
            H_id[i] = 2 * M_PI * Fc;
        }
        else
        {
            H_id[i] = sin(2 * M_PI * Fc * i) / (M_PI * i);
        }
        // Весовая функция Блекмена
        W[i] = 0.42 + 0.5 * cosl((2 * M_PI * i) /(filterSize - 1)) + 0.08 * cosl((4 * M_PI * i) / (filterSize - 1));
        H[i] = H_id[i] * W[i];
    }

    // Нормировка импульсной характеристики
    double SUM = 0;
    for (int i = 0; i < filterSize; i++)
    {
        SUM += H[i];
    }
    for (int i = 0; i < filterSize; i++)
    {
        H[i] /= SUM; // Сумма коэффициентов равна 1
    }
    return H;
}

std::vector<double> FIR::getFilterCoeffs(int filterSize, double cutoffFreq, double sampleRate)
{
    std::vector<double> H(filterSize);    // Импульсная характеристика фильтра
    std::vector<double> H_id(filterSize); // Идеальная импульсная характеристика
    std::vector<double> W(filterSize);    // Весовая функция

    double normalizedCutoffFreq = cutoffFreq / sampleRate;

    for (int i = 0; i < filterSize; i++)
    {
        if (i == 0)
        {
            H_id[i] = 2 * M_PI * normalizedCutoffFreq;
        }
        else
        {
            H_id[i] = sin(2 * M_PI * normalizedCutoffFreq * i) / (M_PI * i);
        }
        // Весовая функция Блекмена
        W[i] = 0.42 + 0.5 * cosl((2 * M_PI * i) /(filterSize - 1)) + 0.08 * cosl((4 * M_PI * i) / (filterSize - 1));
        H[i] = H_id[i] * W[i];
    }

    // Нормировка импульсной характеристики
    double SUM = 0;
    for (int i = 0; i < filterSize; i++)
    {
        SUM += H[i];
    }
    for (int i = 0; i < filterSize; i++)
    {
        H[i] /= SUM; // Сумма коэффициентов равна 1
    }
    return H;
}

std::vector<std::complex<double>> FIR::filter(const std::vector<std::complex<double> > &signal, const std::vector<double> &coeffs)
{
    int sizeIn = signal.size();
    int filterSize = coeffs.size();
    // Фильтрация входных данных
    std::vector<std::complex<double>> out(sizeIn);
    for (int i = 0; i < sizeIn; i++)
    {
        out[i] = 0.0;
        for (int j = 0; j < filterSize - 1; j++) // Формула фильтра
            if(i - j >= 0)
                out[i] += coeffs[j] * signal[i - j];
    }
    return out;
}

template<class T>
std::vector<T> FIR::slidingWindowAverage(const std::vector<T> &array, int win)
{
    size_t sz = array.size();
    size_t newSz = sz;
    std::vector<T> res(newSz,T());

    double curWin;

    for (int i = 0; i < newSz; i++)
    {
        T sum = T();
        curWin = 0;
        for (int j = 0; j < win; j++)
        {
            int curInd = i+j-(win/2);
//            if (bs == BS_CYCLE && curInd < 0)
//                curInd += newSz;
//            if (bs == BS_CYCLE && curInd >= newSz)
//                curInd -= newSz;

            if (curInd >= 0 && curInd < newSz) {
                curWin++;
                sum += array[curInd];
            }
        }
        res[i] = sum/curWin;

    }
    return res;
}

template std::vector<std::complex<double>> FIR::slidingWindowAverage(const std::vector<std::complex<double>> &array, int win);
template std::vector<double> FIR::slidingWindowAverage(const std::vector<double> &array, int win);
