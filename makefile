CC=g++
BOOSTROOT=$(HOME)/computing/boost
OPENCVROOT=$(HOME)/computing/opencv/
CPPINCLUDE=$(HOME)/cpp/include
SDIR=./src
ODIR=./objects
IDIR=./include
CFLAGS=-Wall -I/usr/local/include -I$(CPPINCLUDE) -I$(BOOSTROOT) -I$(IDIR) -std=gnu++14 -c -D_USE_MATH_DEFINES -O3 -fopenmp
LDFLAGS=-L/usr/local/lib -lfftw3 -lm -fopenmp
_SRCS=ScanParams.cpp Pulse.cpp MatResponse.cpp FiberBundle.cpp CalibMat.cpp scan_material.cpp
_HEADS=DataOps.hpp ScanParams.hpp Pulse.hpp MatResponse.hpp Refraction.hpp FiberBundle.hpp CalibMat.hpp scan_material.hpp 

OBJECTS=$(patsubst %,$(ODIR)/%,$(_SRCS:.cpp=.o))
HEADERS=$(patsubst %,$(IDIR)/%,$(_HEADS))
SOURCES=$(patsubst %,$(SDIR)/%,$(_SRCS))
EXECUTABLE ?= scan_material

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

$(ODIR)/%.o: $(SDIR)/%.cpp $(HEADERS) 
	$(CC) $(CFLAGS) -c $< -o $@ 

.PHONY: clean

clean:
	rm $(ODIR)/*.o $(EXECUTABLE)

# DO NOT DELETE
