CC=g++
BOOSTROOT=/opt/boost/include
PYTHONINCLUDE=/usr/include/python3.7m
SDIR=./src
ODIR=./objects
IDIR=./include
FTFLAGS=`pkg-config --cflags --libs fftw3`
CVFLAGS=`pkg-config opencv4 --cflags --libs`
CFLAGS= -Wall -I/usr/local/include -I/usr/local -I/opt/opencv/include/opencv/opencv2 -I$(BOOSTROOT) -I$(IDIR) -I${PYTHONINCLUDE} -std=gnu++17 -c -D_USE_MATH_DEFINES -O3 -fopenmp
LDFLAGS=-L/usr/local/lib -L/use/local -L/opt/opencv/lib -L/opt/opencv/lib64 -lopencv_highgui -lopencv_shape -lopencv_calib3d -lopencv_imgcodecs -lopencv_flann -lopencv_imgproc -lopencv_core -lm -lfftw2 -fopenmp 

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
	$(CC) $(CFLAGS) $(LDFLAGS) -c $< -o $@ 

.PHONY: clean

clean:
	rm $(ODIR)/*.o $(EXECUTABLE)

# DO NOT DELETE
