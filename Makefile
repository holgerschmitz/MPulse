
TARGET_BASE = mpulse

DIMENSIONS = 3

#OFLAGS  = -g -O0 -Wall
OFLAGS  = -O3 -Wall -std=c++14

INCLUDE = -I/usr/local/include -I/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include
CXX     = mpic++
LINK	 = mpic++

CXXFLAGS = $(OFLAGS)

SOURCES = $(wildcard src/*.cpp) \
  huerto/electromagnetics/current.cpp \
  huerto/electromagnetics/em_fields.cpp \
  huerto/electromagnetics/fdtd/fdtd_plain.cpp \
  huerto/maths/functions/core.cpp

BUILD_DIR = build
BIN_DIR = bin

OBJECTS = $(addprefix obj/,$(SOURCES:.cpp=.o))

LDFLAGS = -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib -Wl,-rpath,/usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib

LOADLIBS = -lhdf5 -lschnek -lfftw3 -lm

DIM1_FLAGS = -DHUERTO_ONE_DIM
DIM2_FLAGS = -DHUERTO_TWO_DIM
DIM3_FLAGS = -DHUERTO_THREE_DIM

FULLTARGET = $(foreach dimension,$(DIMENSIONS),$(BIN_DIR)/$(TARGET_BASE)$(dimension)d)

all: $(FULLTARGET)

define PROGRAM_template =
 TARGET$(1)D_OBJS = $(addprefix $(BUILD_DIR)/$(1)d/,$(patsubst %.cpp,%.o,$(SOURCES)))
 $(BIN_DIR)/$(TARGET_BASE)$(1)d: $$(TARGET$(1)D_OBJS)
	@mkdir -p $(BIN_DIR)
	$(LINK) $$^ -o $$@ $(OFLAGS) $(LDFLAGS) $(LOADLIBS)
 $$(TARGET$(1)D_OBJS): $(BUILD_DIR)/$(1)d/%.o: %.cpp
	@mkdir -p $$(dir $$@)
	$(CXX) -o $$@ -c $(CXXFLAGS) $(INCLUDE) $$(DIM$(1)_FLAGS) $$<
 ALL_OBJS   += $$(TARGET$(1)D_OBJS)
endef

$(foreach dimension,$(DIMENSIONS),$(eval $(call PROGRAM_template,$(dimension))))



clean:
	-rm -f $(OBJECTS) core $(FULLTARGET)


