
#include "plotter.h"

Plotter::Plotter()
{

}

void Plotter::initPlot(QCustomPlot *plot, std::string &xLabel, std::string &yLabel)
{
    plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables | QCP::iSelectAxes);
    plot->setSelectionRectMode(QCP::SelectionRectMode::srmNone);

    plot->setAutoAddPlottableToLegend(false);
    //    plot->legend->setVisible(true);

    plot->xAxis->setLabel(xLabel.c_str());
    plot->yAxis->setLabel(yLabel.c_str());

    plot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignLeft|Qt::AlignTop);

    //    plot->addGraph(ui->signalPlot->xAxis, ui->signalPlot->yAxis);
}

void Plotter::plot(std::vector<std::complex<double>> &signal, QCustomPlot *plot, const QColor &realColor, const QColor &imagColor, int width)
{
    QCPGraph *real = new QCPGraph(plot->xAxis, plot->yAxis);
    QCPGraph *imag = new QCPGraph(plot->xAxis, plot->yAxis);
    real->addToLegend();
    real->setPen(QPen(realColor, width));
    imag->addToLegend();
    imag->setPen(QPen(imagColor, width));

    int dataIndex = 0;
    auto end = signal.end();
    for (auto begin = signal.begin(); begin != end; ++begin)
    {
        real->addData(dataIndex++, begin->real());
        imag->addData(dataIndex++, begin->imag());
    }

    if (signal.size() > 0)
    {
        plot->xAxis->setRange(0, dataIndex);
        double min = 0, max = 0;
        for (auto s : signal)
        {
            if (s.real() < min)
                min = s.real();
            if (s.imag() < min)
                min = s.imag();
            if (s.real() > max)
                max = s.real();
            if (s.imag() > max)
                max = s.imag();
        }
        plot->yAxis->setRange(min, max);
    }
}

void Plotter::plot(std::vector<double> &signal, QCustomPlot *plot, const QColor &color, int width)
{
    QCPGraph *graph = new QCPGraph(plot->xAxis, plot->yAxis);
    graph->addToLegend();
    graph->setPen(QPen(color, width));

    int dataIndex = 0;
    auto end = signal.end();
    int len = signal.size();
    for (auto begin = signal.begin(); begin != end; ++begin)
    {
        graph->addData(dataIndex++ - len / 2, *begin);
    }

    if (signal.size() > 0)
    {
        plot->xAxis->setRange(-len / 2, len / 2);
        plot->yAxis->setRange(*std::min_element(signal.begin(), signal.end()), *std::max_element(signal.begin(), signal.end()));
    }
}
