
include $(HOME)/src/makeinclude/include.mak

OFLAGS  = $(X_OFLAGS)
#OFLAGS  = $(X_OFLAGS_DEBUG) -DSINGLE_PROCESSOR
INCLUDE = $(X_INC) $(LX_INC) -I/home/monopole/hschmitz/include

CXXFLAGS = $(OFLAGS)


LDFLAGS = $(X_LIB) $(LX_LIB)

EXTRA_HDF = ../../src/$(X_OSTYPE)/hdfstream.o
LOADLIBS =  $(EXTRA_HDF) $(X_HDF_LIBS) -lm

HDFJOINLIBS = $(X_HDF_LIBS) -lboost_program_options -lm 

BINARY_PREFIX = $(X_ROOT_DIR)/bin/$(X_OSTYPE)

CXX = $(X_CXX_MPI)

all: h5tognu h5topnm gnutopnm

h5tognu: h5tognu.o
	$(CXX) $^ -o $@ $(OFLAGS) $(LDFLAGS) $(LOADLIBS)
	cp $@ $(HOME)/bin/

h5topnm: h5topnm.o
	$(CXX) $^ -o $@ $(OFLAGS) $(LDFLAGS) $(LOADLIBS)
	cp $@ $(HOME)/bin/

slicetognu: slicetognu.o
	$(CXX) $^ -o $@ $(OFLAGS) $(LDFLAGS) $(LOADLIBS)
	cp $@ $(HOME)/bin/

gnutopnm: gnutopnm.o
	$(CXX) $^ -o $@ $(OFLAGS) $(LDFLAGS)
	cp $@ $(HOME)/bin/
    
hdfjoin: hdfjoin.o
	$(CXX) $^ -o $@ $(OFLAGS) $(LDFLAGS) $(HDFJOINLIBS) 
	cp $@ $(HOME)/bin/

h5fft2d: h5fft2d.o
	$(CXX) $^ -o $@ $(OFLAGS) $(LDFLAGS) -lhdf5 $(LOADLIBS) -lboost_program_options -lfftw3
	cp $@ $(HOME)/bin/

h5crop: h5crop.o
	$(CXX) $^ -o $@ $(OFLAGS) $(LDFLAGS) -lhdf5 $(LOADLIBS) -lboost_program_options -lfftw3
	cp $@ $(HOME)/bin/

gnutoh5: gnutoh5.o
	$(CXX) $^ -o $@ $(OFLAGS) $(LOADLIBS) $(LDFLAGS)
	cp $@ $(HOME)/bin/


%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $<


clean:
	-rm -f *.o core h5tognu h5topnm


