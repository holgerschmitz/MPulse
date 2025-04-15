export LD_LIBRARY_PATH := /home/terence411/resources/lib:$(LD_LIBRARY_PATH) # resolves the libsz.so issue

# export LD_LIBRARY_PATH=/home/terence411/resources/lib:$LD_LIBRARY_PATH (in terminal)


TARGET_BASE = mpulse

DIMENSIONS = 1 2 3

#OFLAGS  = -g -O0 -Wall -std=c++17
OFLAGS  = -O3 -Wall -std=c++17


INCLUDE = -I/usr/local/include \
          -I/home/terence411/resources/include \
          -I/home/terence411/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-12.3.0/boost-1.82.0-3zvrwkhbsxoaivfmmy2gonv4qwdn36fb/include \
          -I/home/terence411/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-12.3.0/hdf5-1.14.5-6ibf3iy452splddw7bjwtnmu3ayj6q5t/include \
          -I/home/terence411/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-12.3.0/kokkos-4.4.01-pytybmwfndyyu6x2f5uqo5ucjgfojxzp/include

CXX     = mpiCC
LINK	  = mpiCC

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
  huerto/maths/random.cpp \
  huerto/maths/functions/core.cpp \
  huerto/simulation/task.cpp

BUILD_DIR = build
BIN_DIR = bin

LDFLAGS = $(HDF_LDFLAGS) -L/usr/local/lib -Wl,-rpath,/usr/local/lib \
                         -L/home/terence411/resources/lib -Wl,-rpath,/home/terence411/resources/lib \
                         -L/home/terence411/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-12.3.0/hdf5-1.14.5-6ibf3iy452splddw7bjwtnmu3ayj6q5t/lib -Wl,-rpath,/home/terence411/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-12.3.0/hdf5-1.14.5-6ibf3iy452splddw7bjwtnmu3ayj6q5t/lib \
                         -L/home/terence411/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-12.3.0/kokkos-4.4.01-pytybmwfndyyu6x2f5uqo5ucjgfojxzp/lib

LOADLIBS = -lhdf5 -lschnek -lfftw3 -lm -lkokkoscore -lkokkoscontainers -lkokkossimd

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