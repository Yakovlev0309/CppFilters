
#include "fir.h"
#include <cmath>

FIR::FIR()
{

}

std::vector<std::complex<double>> FIR::compensatePhaseDelay(const std::vector<std::complex<double>> &signal, int filterSize)
{
//    int size = signal.size();
//    int halfFilterSize = filterSize / 2;
//    std::vector<std::complex<double>> withoutDelays(size);
//    for (int i = halfFilterSize; i < size - halfFilterSize; i++)
//    {
//        withoutDelays[i - halfFilterSize] = signal[i];
//    }
//    return withoutDelays;

    std::vector<std::complex<double>> withoutDelays = signal;
    auto begin = withoutDelays.begin();
    auto end = withoutDelays.begin() + filterSize / 2;
    withoutDelays.erase(begin, end);

    begin = withoutDelays.end() - filterSize / 2;
    end = withoutDelays.end() - 1;
    withoutDelays.erase(begin, end);

    return withoutDelays;
}

std::vector<double> FIR::getLowPassFilterCoeffs(int filterSize, double cutoffFreq, double sampleRate)
{
    std::vector<double> coefficients(filterSize);

    // Вычисляем центральный индекс
    int center = filterSize / 2;

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

        // Добавляем коэффициент в массив
        coefficients[n] = coeff;
    }

    return coefficients;
}

std::vector<double> FIR::getHighPassFilterCoeffs(int filterSize, double cutoffFreq, double sampleRate)
{
    std::vector<double> coefficients(filterSize);

    // Вычисляем центральный индекс
    int center = filterSize / 2;

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

std::vector<std::complex<double>> FIR::apply_fir_filter(const std::vector<std::complex<double>>& x, const std::vector<double>& h)
{
    int N = x.size();
    int M = h.size();
    std::vector<std::complex<double>> y(N + M - 1, 0);
    for (int n = 0; n < N + M - 1; ++n) {
        for (int k = 0; k < M; ++k) {
            if (n - k >= 0 && n - k < N) {
                y[n] += h[k] * x[n - k];
            }
        }
    }
    return y;
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
