
#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDebug>

#include "signalgen.h"
#include "ft.h"
#include "fir.h"
#include "cic.h"
#include "cic_old.h"
#include "qcustomplot.h"

#include "fft.hpp"
#include "fft_filter.h"

enum class FilterType
{
    LOW_PASS,
    HIGH_PASS
};

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow

{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    virtual void keyPressEvent(QKeyEvent *event) override;
    virtual void keyReleaseEvent(QKeyEvent *event) override;

    void on_clearGraphs_clicked();

    void on_slidingWindowAverage_clicked();

    void on_lpFirFilter_clicked();

    void on_hpFirFilter_clicked();

    void on_lpFftFilter_clicked();

    void on_lpBatterworthFilter_clicked();

    void on_cicFilter_clicked();

    void on_cicDecimator_clicked();

    void on_cicInterpolator_clicked();

private:
    void updateValues();

    void initPlot(QCustomPlot *plot, std::string &xLabel, std::string &yLabel);
    void clearGraphs();
    void replot();

    std::vector<double> absComplex(const std::vector<std::complex<double>> &samples);

    void firFilter(const FilterType &filterType);

private:
    Ui::MainWindow *ui;

    void plot(std::vector<std::complex<double>> &signal, QCustomPlot *plot, const QColor &realColor = Qt::green, const QColor &imagColor = Qt::red, int width=1);
    void plot(std::vector<double> &signal, QCustomPlot *plot, const QColor &color = Qt::blue, int width = 1);

    QCustomPlot *comparePlot;
    QCustomPlot *signalPlot;
    QCustomPlot *filterPlot;
    QCustomPlot *filteredPlot;
    QCustomPlot *fSignalPlot;
    QCustomPlot *fFilterPlot;
    QCustomPlot *fFilteredPlot;
    QCustomPlot *mSignalPlot;
    QCustomPlot *mFilteredPlot;
    std::vector<QCustomPlot*> plots;

    QButtonGroup *signalType;

    int samplesCount;
    double ampl;
    double freq;
    double sampleRate;
    int noiseCount;
    double freqFactor;

    int filterSize;
    double cutoffFreq;

    int R, N, M;
};

#endif // MAINWINDOW_H
