
TARGET=mpulse

#OFLAGS  = -g -O0 -Wall
OFLAGS  = -O3 -Wall

INCLUDE = -I/usr/local/include -I/usr/include/hdf5/openmpi
#CXX     = $(X_CXX)
CXX     = mpic++

CXXFLAGS = $(OFLAGS)

SOURCES = src/border.cpp \
  src/cpml_border.cpp \
  src/current.cpp \
  src/diagnostic.cpp \
  src/em_fields.cpp \
  src/fdtd_plain.cpp \
  src/fdtd_plrc.cpp \
  src/incsource.cpp \
  src/plasmacurrent.cpp \
  src/shortpulsefunctions.cpp \
  src/shortpulseinject.cpp \
  src/sources.cpp \
  src/specfunc.cpp \
  src/mpulse.cpp


OBJECTS = $(SOURCES:.cpp=.o)

LDFLAGS = -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi

LOADLIBS = -lhdf5 -lschnek -lm
BINDIR = bin
OBJDIR = obj

FULLTARGET = $(BINDIR)/$(TARGET)

all: $(FULLTARGET)

$(FULLTARGET): $(OBJECTS) 
	@mkdir -p $(BINDIR)
	$(CXX) $^ -o $@ $(OFLAGS) $(LDFLAGS) $(LOADLIBS)


%.o: %.cpp
	$(CXX) -o $@ -c $(CXXFLAGS) $(INCLUDE) $<


clean:
	-rm -f $(OBJECTS) core $(FULLTARGET)


