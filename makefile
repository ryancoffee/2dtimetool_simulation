CC=g++
BOOSTROOT=/usr/include/boost/
PYTHONINCLUDE=/usr/include/python3.10
SDIR=./src
ODIR=./objects
IDIR=./include
FTFLAGS=`pkg-config fftw3 --cflags`
CVFLAGS=`pkg-config opencv4 --cflags`
FTLIBS=`pkg-config fftw3 --libs`
CVLIBS=`pkg-config opencv4 --libs`
CFLAGS= -Wall -I/usr/local/include -I/usr/local -I/usr/local/opencv4/ -I$(BOOSTROOT) -I$(IDIR) -I${PYTHONINCLUDE} ${CVFLAGS} ${FTFLAGS}-std=gnu++17 -c -D_USE_MATH_DEFINES -O3 -fopenmp
#CFLAGS= -Wall -I/usr/local/include -I/usr/local -I/opt/opencv/include/opencv4/opencv2 -I/opt/opencv/include/opencv4 -I$(BOOSTROOT) -I$(IDIR) -I${PYTHONINCLUDE} ${CVFLAGS} ${FTFLAGS}-std=gnu++17 -c -D_USE_MATH_DEFINES -O3 -fopenmp
LDFLAGS=-L/usr/local/lib -L/usr/local ${FTLIBS} ${CVLIBS} -lm  -fopenmp 
#LDFLAGS=-L/usr/local/lib -L/use/local -L/opt/opencv/lib -L/opt/opencv/lib64 -lopencv_stitching -lopencv_aruco -lopencv_bgsegm -lopencv_bioinspired -lopencv_ccalib -lopencv_dpm -lopencv_face -lopencv_freetype -lopencv_fuzzy -lopencv_hdf -lopencv_line_descriptor -lopencv_reg -lopencv_rgbd -lopencv_saliency -lopencv_stereo -lopencv_structured_light -lopencv_phase_unwrapping -lopencv_superres -lopencv_optflow -lopencv_surface_matching -lopencv_datasets -lopencv_text -lopencv_highgui -lopencv_plot -lopencv_videostab -lopencv_video -lopencv_videoio -lopencv_viz -lopencv_shape -lopencv_ml -lopencv_ximgproc -lopencv_xobjdetect -lopencv_objdetect -lopencv_calib3d -lopencv_imgcodecs -lopencv_flann -lopencv_xphoto -lopencv_photo -lopencv_imgproc -lopencv_core -lm -lfftw3 -fopenmp 
#LDFLAGS=-L/usr/local/lib -L/use/local -L/opt/opencv/lib -L/opt/opencv/lib64 -lopencv_gapi -lopencv_stitching -lopencv_aruco -lopencv_bgsegm -lopencv_bioinspired -lopencv_ccalib -lopencv_cvv -lopencv_dnn_objdetect -lopencv_dpm -lopencv_face -lopencv_freetype -lopencv_fuzzy -lopencv_hdf -lopencv_hfs -lopencv_img_hash -lopencv_line_descriptor -lopencv_quality -lopencv_reg -lopencv_rgbd -lopencv_saliency -lopencv_stereo -lopencv_structured_light -lopencv_phase_unwrapping -lopencv_superres -lopencv_optflow -lopencv_surface_matching -lopencv_tracking -lopencv_datasets -lopencv_text -lopencv_highgui -lopencv_dnn -lopencv_plot -lopencv_videostab -lopencv_video -lopencv_videoio -lopencv_viz -lopencv_xfeatures2d -lopencv_shape -lopencv_ml -lopencv_ximgproc -lopencv_xobjdetect -lopencv_objdetect -lopencv_calib3d -lopencv_imgcodecs -lopencv_features2d -lopencv_flann -lopencv_xphoto -lopencv_photo -lopencv_imgproc -lopencv_core -lm -lfftw3 -fopenmp 

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
