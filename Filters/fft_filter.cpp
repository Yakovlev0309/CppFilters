
#include "fft_filter.h"
#include "fft.hpp"

FFT_Filter::FFT_Filter()
{

}

int pow2(int len)
{

    char pow = 0;
    int num = 1;
    while(num < len && num*2 - num <= len-num)
    {
        num *= 2;
        pow++;
    }

    return pow;
}

int nextpow2(int len){

    char pow = 0;
    int num = 1;
    while(num < len)
    {
        num *= 2;
        pow++;
    }

    return pow;
}

bool isntPow2(int len)
{

    return (pow(2, pow2(len)) != len);
}

std::complex<double> plusd(std::complex<double> a, std::complex<double> b){

//    std::complex<double> out;
//    out = {a.real() + b.real(), a.imag() + b.imag()};
//    return out;
    return a + b;
}

std::complex<double> minusd(std::complex<double> a, std::complex<double> b){

//    std::complex<double> out;
//    out = {a.real() - b.real(), a.imag() - b.imag()};
//    return out;
    return a - b;
}

std::complex<double> multd(std::complex<double> a, std::complex<double> b){

//    std::complex<double> out;
//    out = {a.real() * b.real() - a.imag() * b.imag(), a.real() * b.imag() + a.imag() * b.real()};
//    return out;
    return a * b;
}

int FFT_Filter::fftWithoutWindow(int len, std::complex<double> *in, std::complex<double> *out, bool forward)
{
    if(NULL == in || NULL == out)
    {
        return -1;
    }

    if(isntPow2(len))
    {
        return -2;
    }
    if(out != in)
    {
        memcpy((void*)out,(void*)in,len*sizeof(std::complex<double>));
    }

    int pow = nextpow2(len);
    long i1, j, k, i2, l1, l2;

    std::complex<double> t,u,c;
    std::complex<double> * outd = (std::complex<double> *)calloc(len,sizeof(std::complex<double>));

    for(int i=0;i<len;i++)
    {
        outd[i] = out[i];
    }

    i2 = len >> 1;
    j = 0;

    for(long i = 0; i < len-1; i++)
    {
        if (i < j)
        {
            std::swap(outd[i],outd[j]);
        }

        k = i2;
        while(k <= j)
        {
            j -= k;
            k >>= 1;
        }
        j += k;
    }

    c = {-1.0, 0.0};
    l2 = 1;
    for(long l = 0; l < pow; l++)
    {
        l1 = l2;
        l2 <<= 1;
        u = {1.0, 0.0};
        for(long j = 0; j < l1; j++)
        {
            for(long i = j; i < len; i += l2)
            {
                i1 = i + l1;
                t=multd(u,outd[i1]);
                outd[i1]=minusd(outd[i],t);
                outd[i]=plusd(outd[i],t);
            }

            u=multd(u,c);
        }

        double real, imag;
        (forward) ? (imag = -sqrt((1.0 - c.real()) / 2.0)) : (imag = sqrt((1.0 - c.real()) / 2.0));
        real = sqrt((1.0 + c.real()) / 2.0);
        c = {real, imag};
    }

    if(!forward)
    {
#pragma omp parallel for
        for(int i = 0; i < len; i++)
        {
            outd[i] /= len;
        }


        for(int i=0;i<len;i++)
        {
            out[i] = outd[i];
        }

        free(outd);
        return 1;
    }
}

int FFT_Filter::filtration(FILTRATION_TYPE ftype, int len, std::complex<double> *in, std::complex<double> *out, double low, double high, double samplingRate, int windowSize)
{
    // Фильтрация сигнала в одноканальном режиме
    int window_size = windowSize;


    if(nullptr == in || nullptr == out || low < -(samplingRate/2) || high <= -(samplingRate/2) || high < low || samplingRate <=0 || high > (samplingRate/2))
    {
        return -1;
    }

    if(isntPow2(window_size))
    {
        return -2;
    }

    double HzToSample = samplingRate/window_size;

    low+=(samplingRate/2.);//так как на вход могут поступать отрицательные частоты, а отрицательных индексов у нас нет
    high+=(samplingRate/2.);//так как на вход могут поступать отрицательные частоты, а отрицательных индексов у нас нет

    int sampleLowPos  = int(ceil(low/HzToSample));//нижняя частота среза в положительной части спектра
    int sampleHighPos = int(ceil(high/HzToSample));//верхняя частота среза в положительной части спектра
    int sampleLowNeg  = int(window_size - ceil((low/HzToSample)));//нижняя частота среза в отрицательной части спектра
    int sampleHighNeg = int(window_size - ceil(high/HzToSample) + 1);//верхняя частота среза в отрицательной части спектра

    if(FILTRATION_TYPE::bpass == ftype)//для полосового фильтра нет необходимости зеркалить полосу пропускания относительно нуля
    {
        sampleLowNeg = -1;
        sampleHighNeg = -1;
    }

    int threads = 1;
    memcpy((void*)out,(void*)in,len*sizeof(std::complex<double>));//копируем данные для обработки

//    threads = omp_get_num_procs();
//    omp_set_dynamic(0);
//    omp_set_num_threads(threads);

    int slices = len / window_size;
    int cycles = slices/threads;
    if(slices - (cycles*threads) != 0)
    {
        cycles++;
    }

    for(int cycle = 0; cycle < cycles; cycle++)
    {
#pragma omp parallel for
        for(int j = 0; j < threads; j++)
        {
            int offset = threads*window_size*cycle + j*window_size;
            if((offset+window_size)>len)
            {
                continue;
            }

            fftWithoutWindow(window_size, &out[offset], &out[offset], true);
            for(int i = 0; i < window_size; i++)//применяем коэффициенты фильтра
            {
                int index = 0;//для определения границ полосы пропускания необходимо пересчитать индексы, так как после БПФ спектр имеет нулевые частоты по краям а не в центре
                if(i<window_size/2)
                {
                    index = i + window_size/2;
                }
                else
                {
                    index = i - window_size/2;
                }

                if(!((i>=(sampleLowPos) && i<(sampleHighPos-1))||(i>=(sampleHighNeg) && i<(sampleLowNeg-1))))
                {
                    out[index+offset] = 0.0;
                }
            }
            fftWithoutWindow(window_size, &out[offset], &out[offset], false);
        }
    }
    return 1;
}

std::vector<std::complex<double> > FFT_Filter::filter(const std::vector<std::complex<double> > &signal, double cutoffFreq, double sampleRate)
{
    auto spectrum = Fft::fft(signal, true);
    for (int i = 0; i < spectrum.size(); i++)
    {
        if (abs(spectrum[i]) <= cutoffFreq)
        {
            spectrum[i] = 0;
        }
    }
    auto filtered = Fft::fft(spectrum, false);
    return filtered;
}
