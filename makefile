CC=g++
BOOSTROOT=/opt/boost/include
OPENCVROOT=/usr/lib/opencv
PYTHONINCLUDE=/usr/include/python3.7m
SDIR=./src
ODIR=./objects
IDIR=./include
CFLAGS= -Wall -I/usr/local/include -I/usr/local -I${OPENCVROOT} -I$(BOOSTROOT) -I$(IDIR) -I${PYTHONINCLUDE} -std=gnu++14 -c -D_USE_MATH_DEFINES -O3 -fopenmp
LDFLAGS=-L/usr/local/lib -L/use/local -L/usr/lib/opencv -lm -lfftw3 -fopenmp 
CVFLAGS= -ggdb `pkg-config --cflags --libs opencv`
_SRCS=ScanParams.cpp Pulse.cpp MatResponse.cpp FiberBundle.cpp CalibMat.cpp scan_material.cpp
_HEADS=Constants.hpp DataOps.hpp ScanParams.hpp Pulse.hpp MatResponse.hpp Refraction.hpp FiberBundle.hpp CalibMat.hpp scan_material.hpp 

OBJECTS=$(patsubst %,$(ODIR)/%,$(_SRCS:.cpp=.o))
HEADERS=$(patsubst %,$(IDIR)/%,$(_HEADS))
SOURCES=$(patsubst %,$(SDIR)/%,$(_SRCS))
EXECUTABLE ?= scan_material

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) $(CVFLAGS) -o $@

$(ODIR)/%.o: $(SDIR)/%.cpp $(HEADERS) 
	$(CC) $(CFLAGS) -c $< -o $@ 

.PHONY: clean

clean:
	rm $(ODIR)/*.o $(EXECUTABLE)

# DO NOT DELETE
