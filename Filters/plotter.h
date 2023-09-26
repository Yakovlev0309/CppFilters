
#ifndef PLOTTER_H
#define PLOTTER_H

#include "qcustomplot.h"
#include <complex>

class Plotter
{
public:
    Plotter();

    void initPlot(QCustomPlot *plot, std::string &xLabel, std::string &yLabel);
    void plot(std::vector<std::complex<double>> &signal, QCustomPlot *plot, const QColor &realColor = Qt::green, const QColor &imagColor = Qt::red, int width=1);
    void plot(std::vector<double> &signal, QCustomPlot *plot, const QColor &color = Qt::blue, int width = 1);

};

#endif // PLOTTER_H
