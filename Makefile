
OFLAGS  = -O3 -Wall -std=c++17

INCLUDE = -I/usr/local/include -I/home/vol07/scarf237/arch/amd/include

CXX     = mpiCC
LINK	 = mpiCC

CXXFLAGS = $(OFLAGS)

SOURCES = $(wildcard src/*.cpp)
BUILD_DIR = build
BIN_DIR = bin

LDFLAGS = -L/home/vol07/scarf237/arch/amd/lib -Wl,-rpath,/home/vol07/scarf237/arch/amd/lib

LOADLIBS = -lhdf5 -lschnek -lm

FULLTARGET = 

all: bin/mpulse

TARGET_OBJS = $(addprefix $(BUILD_DIR)/,$(patsubst %.cpp,%.o,$(SOURCES)))

bin/mpulse: $(TARGET_OBJS)
	@mkdir -p $(BIN_DIR)
	$(LINK) $^ -o $@ $(OFLAGS) $(LDFLAGS) $(LOADLIBS)

$(TARGET_OBJS): $(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) -o $@ -c $(CXXFLAGS) $(INCLUDE) $<



clean:
	-rm -f $(TARGET_OBJS) core $(FULLTARGET)


