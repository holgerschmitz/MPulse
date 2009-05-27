
include $(HOME)/src/makeinclude/include.mak

OFLAGS  = $(X_OFLAGS) -Wall
#OFLAGS  = $(X_OFLAGS_BRUTE) -Wall
#OFLAGS  = $(X_OFLAGS_DEBUG) -Wall
INCLUDE = $(X_INC) $(X_BOOST_INC) $(X_HDF_INC) $(LX_INC) $(X_FFTW3_INC)
CXX     = $(X_CXX)

CXXFLAGS = $(OFLAGS)

MPULSE_SOURCES = hdfstream.cpp \
	cpml_border.cpp \
	diagnostic.cpp \
	factor.cpp \
	fdtd_plain.cpp \
	fdtd_nl.cpp \
	fdtd_disp.cpp \
	fdtd_plrc.cpp \
	fielddiag.cpp \
	fieldsim.cpp \
	freqdiag.cpp \
	gaussinit.cpp \
	gaussinject.cpp \
	globals.cpp \
	incsource.cpp \
	mpi_bound.cpp \
	mpulse.cpp \
	periodic_bound.cpp \
	plasmacurrent.cpp \
	plasmadensity.cpp \
	process.cpp \
	pulseinit.cpp \
	rebuild.cpp \
	shortpulseinject.cpp \
	shortpulsefunctions.cpp \
	sources.cpp \
	specfunc.cpp \
	storage.cpp \
	waveinit.cpp


MPULSE_OBJECTS = hdfstream.o \
	cpml_border.o \
	diagnostic.o \
	factor.o \
	fdtd_plain.o \
	fdtd_nl.o \
	fdtd_disp.o \
	fdtd_plrc.o \
	fielddiag.o \
	fieldsim.o \
	freqdiag.o \
	gaussinit.o \
	gaussinject.o \
	globals.o \
	incsource.o \
	mpi_bound.o \
	mpulse.o \
	periodic_bound.o \
	plasmacurrent.o \
	plasmadensity.o \
	process.o \
	pulseinit.o \
	rebuild.o \
	shortpulseinject.o \
	shortpulsefunctions.o \
	sources.o \
	specfunc.o \
	storage.o \
	waveinit.o

LDFLAGS = $(X_LIB) $(LX_LIB)

LOADLIBS = $(X_HDF_LIBS) $(X_FFTW3_LIB) -lm

BINARY_PREFIX = $(X_ROOT_DIR)/bin/$(X_OSTYPE)

all: mpulse

mpulse: $(MPULSE_OBJECTS) 
	$(CXX) $^ -o $@ $(OFLAGS) $(LDFLAGS) $(LOADLIBS)
	cp $@ $(HOME)/bin/$(X_OSTYPE)/mpulse


%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $<


clean:
	-rm -f *.o core mpulse


