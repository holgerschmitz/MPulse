
include $(HOME)/src/makeinclude/include.mak

#OFLAGS  = $(X_OFLAGS) -DNDEBUG
OFLAGS  = $(X_OFLAGS) 
INCLUDE = $(X_INC)



LDFLAGS = $(X_LIB) -L/usr/local/lib

LOADLIBS = $(X_HDF_LIBS) -lschnek -lm

BINARY_PREFIX = $(X_ROOT_DIR)/bin/$(X_OSTYPE)

all: testhdf

testhdf: testhdf.o
	$(CXX) $^ ../src/hdfstream.o -o $@ $(OFLAGS) $(LDFLAGS) $(LOADLIBS)


%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $<


clean:
	-rm -f *.o core testhdf


