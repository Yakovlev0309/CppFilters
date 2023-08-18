#include "fft.hpp"
#include <cmath>

template<typename T>
std::vector<std::complex<T> > Fft::fft(const std::vector<T>& data, bool forward)
{
    int sz = data.size();
    if (sz <= 0)
        throw std::logic_error("trying to get fft from empty array");
    int newSizeDegree = (int)ceil(log((double)sz)/log(2));
    int newSize = pow(2,newSizeDegree);

    std::vector<std::complex<T> > res(newSize);
    fft(data,res,newSize,forward);
    return res;
}

template<typename T>
std::vector<std::complex<T> > Fft::fft(const std::vector<std::complex<T> > &data, bool forward)
{
    int sz = data.size();
    if (sz <= 0)
        throw std::logic_error("trying to get fft from empty array");
    int newSizeDegree = (int)ceil(log((double)sz)/log(2));
    int newSize = pow(2,newSizeDegree);

    std::vector<std::complex<T> > res(newSize);
    fft(data,res,newSize,forward);
    return res;
}

template<class T>
void fft_impl(std::vector<std::complex<T> > &out, size_t len, bool forward){
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

template<class T>
void Fft::fft(const std::vector<T> &data, std::vector<std::complex<T> > &out, size_t len, bool forward)
{
    int dataSz = data.size();
    for(unsigned long i = 0; i < dataSz; i++){
        out[i] = data[i];
    }
    fft_impl(out,len,forward);
}

template<class T>
void Fft::fft(const std::vector<std::complex<T> > &data, std::vector<std::complex<T> > &out, size_t len, bool forward)
{
    int dataSz = data.size();
    for(unsigned long i = 0; i < dataSz; i++){
        out[i] = data[i];
    }
    fft_impl(out,len,forward);
}

template std::vector<std::complex<double> > Fft::fft<double >(const std::vector<double>&,bool);
template std::vector<std::complex<double> > Fft::fft<double >(const std::vector<std::complex<double> >&,bool);
template void Fft::fft<double >(const std::vector<double>&, std::vector<std::complex<double> >&,size_t,bool);
template void Fft::fft<double >(const std::vector<std::complex<double> >&, std::vector<std::complex<double> >&,size_t,bool);

void Dft::dftNear(const std::vector<std::complex<double> > &signal, std::vector<double>& out, double peakCarrier, int halfWin, double carrierStep, double sampleRate, std::function<double (int)> &func)
{

    bool forward = true;
    int len = signal.size();

    int resSize = halfWin*2+1;
//    std::vector<double> res(halfWin*2+1,0);

    double tmp = 1.;
    (forward) ? (tmp = -1.) : (tmp = tmp);

    double leftSide = ( double(-halfWin)*carrierStep+peakCarrier ) / sampleRate * M_PI;
    double rightSide =( double(halfWin)*carrierStep+peakCarrier ) / sampleRate * M_PI;

    func = [rightSide,resSize,sampleRate,leftSide](int i){
        double res = ( double(i)/(resSize-1)*rightSide + double(resSize-i-1)/(resSize-1) * leftSide ) / M_PI * sampleRate;
        if (res > sampleRate)
            res -= 2.*sampleRate;
        if (res < -sampleRate)
            res += 2.*sampleRate;
        return res;
    };

    for(int i=0;i<resSize;i++)
    {
        double outRe = 0;
        double outIm = 0;
        double mul = tmp * ( double(i)/(resSize-1)*rightSide + double(resSize-i-1)/(resSize-1) * leftSide );
        for(int j=0;j<len;j++)
        {
            double inRe=signal[j].real();
            double inIm=signal[j].imag();
            outRe+=inRe * cos(mul * j) - inIm * sin(mul * j);
            outIm+=inRe * sin(mul * j) + inIm * cos(mul * j);
        }

        if(!forward)
        {
            outRe /= len;
            outIm /= len;
        }
        out[i] = std::abs(std::complex<double>(outRe,outIm));

    }
}
