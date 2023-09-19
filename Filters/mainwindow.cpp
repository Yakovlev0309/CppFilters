
#include "mainwindow.h"
#include "ui_mainwindow.h"


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

//   ui->plot->setOpenGl(true);

    comparePlot = ui->comparePlot;
    signalPlot = ui->signalPlot;
    filterPlot = ui->filterPlot;
    filteredPlot = ui->filteredPlot;
    fSignalPlot = ui->fSignalPlot;
    fFilterPlot = ui->fFilterPlot;
    fFilteredPlot = ui->fFilteredPlot;
    mSignalPlot = ui->mSignalPlot;
    mFilteredPlot = ui->mFilteredPlot;

    plots.push_back(comparePlot);
    plots.push_back(signalPlot);
    plots.push_back(filterPlot);
    plots.push_back(filteredPlot);
    plots.push_back(fSignalPlot);
    plots.push_back(fFilterPlot);
    plots.push_back(fFilteredPlot);
    plots.push_back(mSignalPlot);
    plots.push_back(mFilteredPlot);

    std::string xLabel = "Время", yLabel = "Амплитуда";
    for (int i = 0; i < plots.size(); i++)
    {
        if (i == 4)
            xLabel = "Частота";
        else if (i == 7)
            yLabel = "дБ";
        initPlot(plots[i], xLabel, yLabel);
    }

    signalType = new QButtonGroup(this);
    signalType->addButton(ui->sinRB);
    ui->sinRB->setChecked(true);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::plot(std::vector<std::complex<double>> &signal, QCustomPlot *plot, const QColor &realColor, const QColor &imagColor, int width)
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

void MainWindow::plot(std::vector<double> &signal, QCustomPlot *plot, const QColor &color, int width)
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

void MainWindow::keyPressEvent(QKeyEvent *event)
{
    for (auto plot : qAsConst(plots))
    {
        if(event->key() == Qt::Key_Shift)
        {
            plot->axisRects().at(0)->setRangeZoom(Qt::Horizontal);
        }
        if(event->key() == Qt::Key_Control)
        {
            plot->axisRects().at(0)->setRangeZoom(Qt::Vertical);
        }
    }
}

void MainWindow::keyReleaseEvent(QKeyEvent *event)
{
    for (auto plot : qAsConst(plots))
    {
        if(event->key() == Qt::Key_Shift || event->key() == Qt::Key_Control)
        {
            plot->axisRects().at(0)->setRangeZoom(Qt::Horizontal | Qt::Vertical);
        }
    }
}

void MainWindow::on_clearGraphs_clicked()
{
    for (int i = 0; i < plots.size(); i++)
    {
        plots[i]->clearGraphs();
        plots[i]->replot();
    }
}

void MainWindow::initPlot(QCustomPlot *plot, std::string &xLabel, std::string &yLabel)
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

void MainWindow::clearGraphs()
{
    for (auto plot : qAsConst(plots))
        plot->clearGraphs();
}

void MainWindow::replot()
{
    for (auto plot : qAsConst(plots))
        plot->replot();
}

void MainWindow::on_slidingWindowAverage_clicked()
{
    updateValues();
    clearGraphs();

    std::vector<std::complex<double>> signal; // Сигнал
    switch (signalType->checkedId())
    {
    default:
        signal = SignalGen::getSin(samplesCount, ampl, freq, sampleRate);
        break;
    }

    signal = SignalGen::addSomeNoise(signal, amplFactor, freq, noiseCount, freqFactor);
    plot(signal, signalPlot);

    auto filtered = FIR::slidingWindowAverage(signal, filterSize);
//    auto filtered = CIC::moving_average_filter(signal, filterSize);
    plot(filtered, filteredPlot);

    auto signalSpectrum = absComplex(Fft::fft(signal, true));
    plot(signalSpectrum, fSignalPlot);
    auto filteredSpectrum = absComplex(Fft::fft(filtered, true));
    plot(filteredSpectrum, fFilteredPlot);

    std::vector<double> mSignal; // Логарифмированный спектр сигнала
    for (auto f : qAsConst(signalSpectrum))
    {
        mSignal.push_back(20 * log10(abs(f)));
    }
    plot(mSignal, mSignalPlot);

    plot(signal, comparePlot);
    plot(filtered, comparePlot, Qt::GlobalColor::darkGreen, Qt::GlobalColor::darkRed, 3);

    replot();
}

void MainWindow::on_lpFirFilter_clicked()
{
    updateValues();
    clearGraphs();

    firFilter(FilterType::LOW_PASS);

    replot();
}

void MainWindow::on_hpFirFilter_clicked()
{
    updateValues();
    clearGraphs();

    firFilter(FilterType::HIGH_PASS);

    replot();
}

void MainWindow::on_bpFirFilter_clicked()
{
    updateValues();
    clearGraphs();

    firFilter(FilterType::BAND_PASS);

    replot();
}

void MainWindow::updateValues()
{
    samplesCount = ui->samplesCount->text().toInt();
    ampl = ui->ampl->text().toDouble();
    freq = ui->freq->text().toDouble();
    sampleRate = ui->sampleRate->text().toDouble();
    noiseCount = ui->noiseCount->text().toInt();
    freqFactor = ui->freqFactor->text().toDouble();
    amplFactor = ui->amplFactor->text().toDouble();

    filterSize = ui->filterSize->text().toInt();
    lowCutoffFreq = ui->lowCutoffFreq->text().toDouble();
    highCutoffFreq = ui->highCutoffFreq->text().toDouble();

    R = ui->R->text().toInt();
    N = ui->N->text().toInt();
    M = ui->M->text().toInt();
}

std::vector<double> MainWindow::absComplex(const std::vector<std::complex<double>> &samples)
{
    int size = samples.size();
    std::vector<double> res(size);
    for (int i = 0; i < size; i++)
    {
        res[i] = abs(samples[i]);
    }
    return res;
}

void MainWindow::firFilter(const FilterType &filterType) // КИХ-фильтр
{
    std::vector<std::complex<double>> signal; // Сигнал
    switch (signalType->checkedId())
    {
    default:
        signal = SignalGen::getSin(samplesCount, ampl, freq, sampleRate);
        break;
    }
    signal = SignalGen::addSomeNoise(signal, amplFactor, freq, noiseCount, freqFactor);
    plot(signal, signalPlot);

    std::vector<double> coeffs;
    switch (filterType)
    {
    case FilterType::LOW_PASS:
        coeffs = FIR::getBandPassFilterCoeffs(filterSize, 0, highCutoffFreq, sampleRate);
//        coeffs = FIR::getLowPassFilterCoeffs(filterSize, cutoffFreq, sampleRate);
        break;
    case FilterType::HIGH_PASS:
        coeffs = FIR::getBandPassFilterCoeffs(filterSize, lowCutoffFreq, sampleRate / 2, sampleRate);
//        coeffs = FIR::getHighPassFilterCoeffs(filterSize, cutoffFreq, sampleRate);
        break;
    default:
        coeffs = FIR::getBandPassFilterCoeffs(filterSize, lowCutoffFreq, highCutoffFreq, sampleRate);
        break;
    }

    plot(coeffs, filterPlot);
    
    auto filtered = FIR::applyFilter(signal, coeffs);
//    auto filtered = FIR::filter(signal, coeffs);
    filtered = FIR::compensatePhaseDelay(filtered, filterSize);
    plot(filtered, filteredPlot);

    auto signalSpectrum = absComplex(Fft::fft(signal, true));
//    signalSpectrum = FIR::slidingWindowAverage(signalSpectrum, 13);
    plot(signalSpectrum, fSignalPlot);
    auto filteredSpectrum = absComplex(Fft::fft(filtered, true));
//    filteredSpectrum = FIR::slidingWindowAverage(filteredSpectrum, 13);
    plot(filteredSpectrum, fFilteredPlot);
    auto coeffsSpectrum = absComplex(Fft::fft(coeffs, true));
    plot(coeffsSpectrum, fFilterPlot, Qt::GlobalColor::darkRed);

    std::vector<double> mSignal; // Логарифмированный спектр сигнала
    for (auto f : qAsConst(signalSpectrum))
    {
        mSignal.push_back(20 * log10(abs(f)));
    }
    plot(mSignal, mSignalPlot);

    std::vector<double> mFiltered; // Логарифмированный спектр фильтра
    for (auto f : qAsConst(coeffsSpectrum))
    {
        mFiltered.push_back(20 * log10(abs(f)));
    }
    plot(mFiltered, mFilteredPlot);

    plot(signal, comparePlot);
    plot(filtered, comparePlot, Qt::GlobalColor::darkGreen, Qt::GlobalColor::darkRed, 3);
}

void MainWindow::on_lpFftFilter_clicked()
{
    updateValues();
    clearGraphs();

    std::vector<std::complex<double>> signal; // Сигнал
    switch (signalType->checkedId())
    {
    default:
        signal = SignalGen::getSin(samplesCount, ampl, freq, sampleRate);
        break;
    }
    signal = SignalGen::addSomeNoise(signal, amplFactor, freq, noiseCount, freqFactor);
    plot(signal, signalPlot);

    int windowSize = filterSize;

    std::vector<std::complex<double>> filtered = signal;
    FFT_Filter::filtration(FILTRATION_TYPE::low, signal.size(), signal.data(), filtered.data(), lowCutoffFreq, highCutoffFreq, sampleRate, windowSize);
    plot(filtered, filteredPlot);

    auto signalSpectrum = absComplex(Fft::fft(signal, true));
    plot(signalSpectrum, fSignalPlot);
    auto filteredSpectrum = absComplex(Fft::fft(filtered, true));
    plot(filteredSpectrum, fFilteredPlot);

    std::vector<double> mSignal; // Логарифмированный спектр сигнала
    for (auto f : qAsConst(signalSpectrum))
    {
        mSignal.push_back(20 * log10(abs(f)));
    }
    plot(mSignal, mSignalPlot);

    plot(signal, comparePlot);
    plot(filtered, comparePlot, Qt::GlobalColor::darkGreen, Qt::GlobalColor::darkRed, 3);

    replot();
}

void MainWindow::on_lpBatterworthFilter_clicked()
{
    updateValues();
    clearGraphs();

    std::vector<std::complex<double>> signal; // Сигнал
    switch (signalType->checkedId())
    {
    default:
        signal = SignalGen::getSin(samplesCount, ampl, freq, sampleRate);
        break;
    }
    signal = SignalGen::addSomeNoise(signal, amplFactor, freq, noiseCount, freqFactor);
    plot(signal, signalPlot);

    auto filtered = FFT_Filter::filter(signal, highCutoffFreq, sampleRate);
    plot(filtered, filteredPlot);

    auto signalSpectrum = absComplex(Fft::fft(signal, true));
    plot(signalSpectrum, fSignalPlot);
    auto filteredSpectrum = absComplex(Fft::fft(filtered, true));
    plot(filteredSpectrum, fFilteredPlot);

    std::vector<double> mSignal; // Логарифмированный спектр сигнала
    for (auto f : qAsConst(signalSpectrum))
    {
        mSignal.push_back(20 * log10(abs(f)));
    }
    plot(mSignal, mSignalPlot);

    plot(signal, comparePlot);
    plot(filtered, comparePlot, Qt::GlobalColor::darkGreen, Qt::GlobalColor::darkRed, 3);

    replot();
}

void MainWindow::on_cicFilter_clicked()
{
    updateValues();
    clearGraphs();

    std::vector<std::complex<double>> signal; // Сигнал
    switch (signalType->checkedId())
    {
    default:
        signal = SignalGen::getSin(samplesCount, ampl, freq, sampleRate);
        break;
    }
    signal = SignalGen::addSomeNoise(signal, amplFactor, freq, noiseCount, freqFactor);
    plot(signal, signalPlot);

//    R = sampleRate / cutoffFreq;
//    ui->R->setText(QString::number(R));

    auto coeffs = FIR::getBandPassFilterCoeffs(filterSize, lowCutoffFreq, highCutoffFreq, sampleRate);
    auto filtered = CIC::cic_decimation_filter(signal, R, M, N);
    filtered = FIR::applyFilter(filtered, coeffs);
    filtered = FIR::compensatePhaseDelay(filtered, filterSize);
    plot(filtered, filteredPlot);

    auto signalSpectrum = absComplex(Fft::fft(signal, true));
    plot(signalSpectrum, fSignalPlot);
    auto filteredSpectrum = absComplex(Fft::fft(filtered, true));
    plot(filteredSpectrum, fFilteredPlot);

    std::vector<double> mSignal; // Логарифмированный спектр сигнала
    for (auto f : qAsConst(signalSpectrum))
    {
        mSignal.push_back(20 * log10(abs(f)));
    }
    plot(mSignal, mSignalPlot);

    plot(signal, comparePlot);
    plot(filtered, comparePlot, Qt::GlobalColor::darkGreen, Qt::GlobalColor::darkRed, 3);

    replot();
}

void MainWindow::on_cicDecimator_clicked()
{
    updateValues();
    clearGraphs();

    std::vector<std::complex<double>> signal; // Сигнал
    switch (signalType->checkedId())
    {
    default:
        signal = SignalGen::getSin(samplesCount, ampl, freq, sampleRate);
        break;
    }
    signal = SignalGen::addSomeNoise(signal, amplFactor, freq, noiseCount, freqFactor);
    plot(signal, signalPlot);

//    R = sampleRate / 2.0 / cutoffFreq;
//    ui->R->setText(QString::number(R));

    auto filtered = CIC::cic_decimation_filter(signal, R, M, N);
//    filtered = CIC::moving_average_filter(filtered, filterSize);
    plot(filtered, filteredPlot);

    auto signalSpectrum = absComplex(Fft::fft(signal, true));
    plot(signalSpectrum, fSignalPlot);
    auto filteredSpectrum = absComplex(Fft::fft(filtered, true));
    plot(filteredSpectrum, fFilteredPlot);

    std::vector<double> mSignal; // Логарифмированный спектр сигнала
    for (auto f : qAsConst(signalSpectrum))
    {
        mSignal.push_back(20 * log10(abs(f)));
    }
    plot(mSignal, mSignalPlot);

    plot(signal, comparePlot);
    plot(filtered, comparePlot, Qt::GlobalColor::darkGreen, Qt::GlobalColor::darkRed, 3);

    replot();
}

void MainWindow::on_cicInterpolator_clicked()
{
    updateValues();
    clearGraphs();

    std::vector<std::complex<double>> signal; // Сигнал
    switch (signalType->checkedId())
    {
    default:
        signal = SignalGen::getSin(samplesCount, ampl, freq, sampleRate);
        break;
    }
    signal = SignalGen::addSomeNoise(signal, amplFactor, freq, noiseCount, freqFactor);
    plot(signal, signalPlot);

    auto filtered = CIC::cic_interpolation_filter(signal, R, M, N);
//    filtered = CIC::moving_average_filter(filtered, filterSize);
    plot(filtered, filteredPlot);

    auto signalSpectrum = absComplex(Fft::fft(signal, true));
    plot(signalSpectrum, fSignalPlot);
    auto filteredSpectrum = absComplex(Fft::fft(filtered, true));
    plot(filteredSpectrum, fFilteredPlot);

    std::vector<double> mSignal; // Логарифмированный спектр сигнала
    for (auto f : qAsConst(signalSpectrum))
    {
        mSignal.push_back(20 * log10(abs(f)));
    }
    plot(mSignal, mSignalPlot);

    plot(signal, comparePlot);
    plot(filtered, comparePlot, Qt::GlobalColor::darkGreen, Qt::GlobalColor::darkRed, 3);

    replot();
}
