CC=g++
PYTHONINCLUDE=/usr/include/python3.10
SDIR=./src
ODIR=./objects
IDIR=./include

# Commenting HDF5 for sake of laptop in flight
#H5LIBS=-lhdf5 -lhdf5_cpp 
#CFLAGS=-Wall -I/usr/local/include -I/usr/local -I/usr/local/opencv4/ -I${IDIR} -I${PYTHONINCLUDE} ${CVFLAGS} ${FTFLAGS} ${H5LIBS} -std=gnu++17 -c -D_USE_MATH_DEFINES -O3 -fopenmp
#LDFLAGS=-L/usr/local -L/usr/local/lib ${FTLIBS} ${CVLIBS} ${H5LIBS} -ldl -lrt -lpthread -lm  -fopenmp 

FTLIBS=-L/usr/local/lib -lfftw3
FTCFLAGS=-I/usr/local/include
#FTFLAGS=`pkg-config fftw3 --cflags`
#FTLIBS=`pkg-config fftw3 --libs`

CFLAGS=-Wall -I/usr/local/include -I/usr/local -I${IDIR} -I${PYTHONINCLUDE} ${FTFLAGS} -std=gnu++17 -c -D_USE_MATH_DEFINES -O3 -fopenmp
LDFLAGS=-L/usr/local -L/usr/local/lib ${FTLIBS} -ldl -lrt -lpthread -lm  -fopenmp 
INCFLAGS=-I/usr/local/include -I$(IDIR) 

_SRCS=ScanParams.cpp Pulse.cpp MatResponse.cpp FiberBundle.cpp CalibMat.cpp scan_material.cpp
_HEADS=Constants.hpp DataOps.hpp ScanParams.hpp Pulse.hpp MatResponse.hpp Refraction.hpp FiberBundle.hpp CalibMat.hpp scan_material.hpp 

OBJECTS=$(patsubst %,$(ODIR)/%,$(_SRCS:.cpp=.o))
HEADERS=$(patsubst %,$(IDIR)/%,$(_HEADS))
SOURCES=$(patsubst %,$(SDIR)/%,$(_SRCS))
EXECUTABLE ?= scan_material

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

$(ODIR)/%.o: $(SDIR)/%.cpp $(HEADERS) 
	$(CC) $(INCFLAGS) $(LDFLAGS) $(CFLAGS) -c $< -o $@ 

.PHONY: clean

clean:
	rm $(ODIR)/*.o $(EXECUTABLE)

# DO NOT DELETE
