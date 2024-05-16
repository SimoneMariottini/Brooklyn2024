ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)

ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

CXX           = g++
CXXFLAGS      = -g -Wall -fPIC
CXXFLAGS      += $(ROOTCFLAGS)

LD            =  g++
LDFLAGS       = -g
SOFLAGS       = -shared $(ROOTGLIBS) -lSpectrum -dynamiclib

NGLIBB        += libAnalyze.so 
NGLIBB        += $(ROOTGLIBS) 
NGLIBB        += -lMinuit -lSpectrum
GLIBB          = $(filter-out -lNew, $(NGLIBB))

Analyze:  Analyze.o Lib
	$(LD) $(LDFLAGS) -o Analyze Analyze.o $(GLIBB)

Analyze.o: Analyze.cpp
	$(CXX) $(CXXFLAGS) -c Analyze.cpp

Waveform.o: Source/Waveform.cc Header/Waveform.h
	$(CXX) $(CXXFLAGS) -c Source/Waveform.cc

Event.o: Source/Event.cc Header/Event.h
	$(CXX) $(CXXFLAGS) -c Source/Event.cc

AnaTools.o: Source/AnaTools.cc Header/AnaTools.h
	$(CXX) $(CXXFLAGS) -c Source/AnaTools.cc

Lib: Analyze.o Waveform.o Event.o AnaTools.o
	$(CXX) $(SOFLAGS) Analyze.o Waveform.o AnaTools.o Event.o -o libAnalyze.so  

clean:
	rm -f *.o  *.so Analyze *~
