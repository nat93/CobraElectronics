TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
#INCLUDEPATH += "/home/andrii/root/include"
INCLUDEPATH += "/home/anatochi/root-6.10.02/include"

#LIBS += -L/home/andrii/root/lib -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread
LIBS += -L/home/anatochi/root-6.10.02/lib -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread

DISTFILES += \
    convertSigToRootFile.bash \
    filelist.dat \
    filelist_RUN_2.dat

SOURCES += \
    convertSigToRootFile.cc \
    src/CobraClass.cc \
    analysis_v1.cc

HEADERS += \
    src/CobraClass.hh \
    include/Constants.hh



