QT       += core gui widgets printsupport

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    cic.cpp \
    fft_filter.cpp \
    fir.cpp \
    main.cpp \
    mainwindow.cpp \
    plotter.cpp \
    signalgen.cpp \
    qcustomplot.cpp \
    fft.cpp

HEADERS += \
    cic.h \
    fft_filter.h \
    fir.h \
    mainwindow.h \
    plotter.h \
    signalgen.h \
    qcustomplot.h \
    fft.hpp

FORMS += \
    mainwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
