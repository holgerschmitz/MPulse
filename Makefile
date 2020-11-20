
TARGET=mpulse

#OFLAGS  = -g -O0 -Wall
OFLAGS  = -O3 -Wall -std=c++14

INCLUDE = -I/usr/local/include -I/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include
#CXX     = $(X_CXX)
CXX     = mpic++

DEFINES = -DHUERTO_THREE_DIM

CXXFLAGS = $(OFLAGS)

SOURCES = src/border.cpp \
  src/cpml_border.cpp \
  src/diagnostic.cpp \
  src/fdtd_plrc.cpp \
  src/focusedpulseinject.cpp \
  src/focusedpulsefunctions.cpp \
  src/incsource.cpp \
  src/plasmacurrent.cpp \
  src/shortpulsefunctions.cpp \
  src/shortpulseinject.cpp \
  src/sources.cpp \
  src/specfunc.cpp \
  src/mpulse.cpp \
  huerto/electromagnetics/current.cpp \
  huerto/electromagnetics/em_fields.cpp \
  huerto/electromagnetics/fdtd/fdtd_plain.cpp \
  huerto/maths/functions/core.cpp


OBJECTS = $(addprefix obj/,$(SOURCES:.cpp=.o))

LDFLAGS = -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib -Wl,-rpath,/usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib

LOADLIBS = -lhdf5 -lschnek -lfftw3 -lm
BINDIR = bin
OBJDIR = obj

FULLTARGET = $(BINDIR)/$(TARGET)

all: $(FULLTARGET)

$(FULLTARGET): $(OBJECTS) 
	@mkdir -p $(BINDIR)
	$(CXX) $^ -o $@ $(OFLAGS) $(LDFLAGS) $(LOADLIBS)


obj/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) -o $@ -c $(CXXFLAGS) $(INCLUDE) $(DEFINES) $<


clean:
	-rm -f $(OBJECTS) core $(FULLTARGET)


