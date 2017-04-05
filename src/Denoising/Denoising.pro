

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = MeshDenoising
TEMPLATE = app

DEFINES += _USE_MATH_DEFINES

SOURCES += main.cpp \
    neighborhood.cpp \
    curvature.cpp \
    color.cpp \
    noise.cpp \
    util.cpp \
    iomesh.cpp \
    metric.cpp \
    denoising.cpp

HEADERS  += \
    mesh.h \
    util.h \
    neighborhood.h \
    curvature.h \
    color.h \
    noise.h \
    iomesh.h \
    metric.h \
    nanoflann.hpp \
    denoising.h \
    custom.h

FORMS    +=

CONFIG += console

include(../LibsInclude.pri)

RESOURCES +=
