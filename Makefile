
TARGET_BASE = mpulse

DIMENSIONS = 1 2 3

#OFLAGS  = -g -O0 -Wall -std=c++14
OFLAGS  = -O3 -Wall -std=c++14

INCLUDE = -I/usr/local/include $(HDF_INCLUDE)

CXX     = mpic++
LINK	 = mpic++

CXXFLAGS = $(OFLAGS)

SOURCES = $(wildcard src/*.cpp) \
  huerto/electromagnetics/current.cpp \
  huerto/electromagnetics/em_fields.cpp \
  huerto/electromagnetics/fdtd/fdtd_plain.cpp \
  huerto/electromagnetics/source/border.cpp \
  huerto/electromagnetics/source/incsource.cpp \
  huerto/electromagnetics/source/plane_wave.cpp \
  huerto/electromagnetics/source/beam.cpp \
  huerto/electromagnetics/pml/cpml_border.cpp \
  huerto/maths/functions/core.cpp \
  huerto/maths/random.cpp \
  huerto/simulation/task.cpp

BUILD_DIR = build
BIN_DIR = bin

LDFLAGS = $(HDF_LDFLAGS)

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
	-rm -f $(ALL_OBJS) core $(FULLTARGET)


