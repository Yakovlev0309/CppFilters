/* (C) 2016 EsonJohn(zhangyibupt@163.com) */

#ifndef CIC_OLD_H
#define CIC_OLD_H

#include "signalgen.h"

class CIC_OLD
{
public:
	// R:Decimation factor;
	// N:Number of sections;
    // M:Differential delay;
    CIC_OLD(int decimationFactor, int numberOfSections, int diferrencialDelay);

    // destructor
    ~CIC_OLD();

	// the actual filter function
	// the parameter input shuld be R-dimensional vector(R continuous samples) and the parameter length should be R
	// the output is double, corresponding to the downsampled output
    double filter(std::vector<double> input, int length);

	// reset the buffer
	void reset();

private:
    int R, N, M;
    std::vector<double> buffer_integrator;	// buffer of the integrator part
    std::vector<std::vector<double>> buffer_comb;		// buffer of the comb part
    std::vector<int> offset_comb;

};

#endif
