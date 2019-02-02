ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

#ROOTCFLAGS    = $(shell /usr/bin/root-config --cflags)
#ROOTLIBS      = $(shell /usr/bin/root-config --libs)
#ROOTGLIBS     = $(shell /usr/bin/root-config --glibs)

CXX           = g++
CXXFLAGS      = -g -Wall -fPIC -Wno-deprecated

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit

CXXFLAGS      += $(ROOTCFLAGS)
CXX           += -I./	
LIBS           = $(ROOTLIBS) 

GLIBS          = $(filter-out -lNew, $(NGLIBS))

CXX	      += -I./obj/
OUTLIB	      = ./obj/
.SUFFIXES: .C
.PREFIXES: ./obj/

#----------------------------------------------------#

all: convertSigToRootFile analysis_v1

convertSigToRootFile: convertSigToRootFile.cc
	g++ -o convertSigToRootFile convertSigToRootFile.cc `root-config --cflags --glibs`

analysis_v1: analysis_v1.cc obj/CobraClass.o
	$(CXX) $(CXXFLAGS) -o analysis_v1 $(OUTLIB)/*.o $(GLIBS) $<

obj/CobraClass.o: src/CobraClass.cc src/CobraClass.hh
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)CobraClass.o $<
	
clean:
	rm -f convertSigToRootFile
	rm -f analysis_v1
	rm -f *~
	rm -f *.root
	rm -f $(OUTLIB)*.o
