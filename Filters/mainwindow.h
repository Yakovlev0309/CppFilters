
#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDebug>

#include "signalgen.h"
#include "fir.h"
#include "cic.h"
#include "fft.hpp"
#include "fft_filter.h"
#include "plotter.h"

enum class FilterType
{
    LOW_PASS,
    HIGH_PASS,
    BAND_PASS
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

    void on_bpFirFilter_clicked();

private:
    void updateValues();

    void clearGraphs();
    void replot();

    std::vector<double> absComplex(const std::vector<std::complex<double>> &samples);

    void firFilter(const FilterType &filterType);

private:
    Ui::MainWindow *ui;

    Plotter plt;

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

    int samplesCount;
    double ampl;
    double freq;
    double sampleRate;
    int noiseCount;
    double freqFactor;
    double amplFactor;

    int filterSize;
    double lowCutoffFreq;
    double highCutoffFreq;

    int R, N, M;
};

#endif // MAINWINDOW_H
