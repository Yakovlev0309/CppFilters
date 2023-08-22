
#include "ft.h"
#include <cmath>

FT::FT()
{

}

std::vector<std::complex<double>> FT::dft(const std::vector<std::complex<double>> &signal)
{
    int N = signal.size();
    std::vector<std::complex<double>> res(N);
    for (int k = 0; k < N; ++k)
    {
        for (int n = 0; n < N; ++n)
        {
            double angle = 2 * M_PI * k * n / N;
            res[k] += signal[n] * cos(angle) - signal[n] * sin(angle);
            res[k] += signal[n] * cos(angle) + signal[n] * sin(angle);
        }
    }
    return res;
}

std::vector<std::complex<double>> FT::forwardDft(const std::vector<std::complex<double>> &signal)
{
    std::vector<std::complex<double>> res(signal.size());

    for (int k = 0; k < signal.size(); ++k)
    {
        res[k] = {0.0, 0.0};

        for (int n = 0; n < signal.size(); ++n) {
            double angle = 2 * M_PI * k * n / signal.size();
            double cosine = cos(angle);
            double sine = sin(angle);

            res[k] = {res[k].real() + (signal[n].real() * cosine) - (signal[n].imag() * sine),
                      res[k].imag() + (signal[n].real() * sine) + (signal[n].imag() * cosine)};
        }

        // Нормализация
        res[k] /= signal.size();
    }
    return res;
}

std::vector<std::complex<double>> FT::fft(std::vector<std::complex<double>> a, bool invert)
{
    int n = (int) a.size();
    for (int i=1, j=0; i<n; ++i) {
        int bit = n >> 1;
        for (; j>=bit; bit>>=1)
        {
            j -= bit;
        }
        j += bit;
        if (i < j)
            swap (a[i], a[j]);
    }

    for (int len=2; len<=n; len<<=1) {
        double ang = 2*M_PI/len * (invert ? -1 : 1);
        std::complex<double> wlen (cos(ang), sin(ang));
        for (int i=0; i<n; i+=len) {
            std::complex<double> w (1);
            for (int j=0; j<len/2; ++j) {
                std::complex<double> u = a[i+j],  v = a[i+j+len/2] * w;
                a[i+j] = u + v;
                a[i+j+len/2] = u - v;
                w *= wlen;
            }
        }
    }
    if (invert)
        for (int i=0; i<n; ++i)
            a[i] /= n;
    return a;
}

template<class T>
void FT::fft_impl(std::vector<std::complex<T>> &out, size_t len, bool forward)
{
    unsigned long i1, j, k, i2, l1, l2;
    int spow = log2(len);

    std::complex<T> t,u,c;

    i2 = len >> 1;
    j = 0;

    for(unsigned long i = 0; i < len-1; i++){
        if (i < j)
            std::swap(out[i], out[j]);
        k = i2;
        while(k <= j){
            j -= k;
            k >>= 1;
        }
        j += k;
    }
    c.real(-1.L);
    c.imag(0.L);
    l2 = 1;
    for(long l = 0; l < spow; l++){
        l1 = l2;
        l2 <<= 1;
        u.real(1.L);
        u.imag(0.L);
        for(unsigned long j = 0; j < l1; j++){
            for(unsigned long i = j; i < len; i += l2){
                i1 = i + l1;
                t=u * out[i1];
                out[i1] = out[i] - t;
                out[i] = out[i] + t;
            }

            u=u *c;
        }

        c.imag((forward?-1:1)*sqrt((1.0 - c.real()) / 2.0));
        c.real(sqrt((1.0 + c.real()) / 2.0));
    }
    if(!forward)
        for(unsigned int i = 0; i < len; i++)
            out[i] /= len;

}

template void FT::fft_impl(std::vector<std::complex<double>> &out, size_t len, bool forward);
