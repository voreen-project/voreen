TEMPLATE = app
TARGET = macdeployqt-voreen
DEPENDPATH += .
INCLUDEPATH += .

# Input
SOURCES += main.cpp shared/shared.cpp
CONFIG += qt warn_on
CONFIG -= app_bundle

