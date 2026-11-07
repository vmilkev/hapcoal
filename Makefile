# hapcoal — SMC / SMC' / MaCS coalescent haplotype simulator (static library)
#
# Targets:
#   make                 lib/libhapcoal.a, bin/hapcoal (CLI) and bin/hapcoal_example (mock client)
#   make lib             only the static library
#   make test            build and run the unit / statistical tests
#   make regression      check that all outputs equal the recorded reference (needs `make data`)
#   make doc             HTML API documentation with Doxygen (if installed) -> doc/html
#   make install         install headers and library under PREFIX (default /usr/local)
#   make print-flags     compiler / linker flags a client program needs
#   make data            download the human chr22 test data into data/
#   make chr22           human chr22 (CEU) validation test case
#   make bench           whole-genome benchmark, 10,000 haplotypes
#   make clean           remove build files (build/, bin/, lib/, doc/)
#   make rclean          remove all results (out/); downloaded data/ is kept
#
# Options (make VAR=value):
#   CXX=g++-15           compiler (default: c++)
#   OPENMP=0|1           enable OpenMP (default 1)
#   ZLIB=0|1             read .vcf.gz / .gz maps through zlib (default 1)
#   NATIVE=0|1           add -march=native (default 0; binaries then not portable)
#   DEBUG=0|1            debug build with assertions and tree/ARG validation after every event
#   PREFIX=/path         installation prefix for `make install`

CXX     ?= c++
AR      ?= ar
OPENMP  ?= 1
ZLIB    ?= 1
NATIVE  ?= 0
DEBUG   ?= 0
PREFIX  ?= /usr/local

# ---------------------------------------------------------------- platform ----
UNAME_S := $(shell uname -s 2>/dev/null || echo Windows)
EXE :=
ifeq ($(OS),Windows_NT)
  EXE := .exe
endif
ifneq (,$(findstring MINGW,$(UNAME_S)))
  EXE := .exe
endif
ifneq (,$(findstring MSYS,$(UNAME_S)))
  EXE := .exe
endif

IS_CLANG := $(shell $(CXX) --version 2>/dev/null | grep -c -i clang)

CXXFLAGS := -std=c++17 -Wall -Wextra -Wno-unused-parameter -MMD -MP
LDFLAGS  :=
LDLIBS   :=

ifeq ($(DEBUG),1)
  CXXFLAGS += -O1 -g -DHAPCOAL_DEBUG
else
  CXXFLAGS += -O3 -DNDEBUG
endif
ifeq ($(NATIVE),1)
  CXXFLAGS += -march=native
endif

# GCC on macOS: point it at the active SDK (Homebrew GCC may be built for an older SDK).
ifeq ($(UNAME_S),Darwin)
  ifeq ($(IS_CLANG),0)
    SDKROOT ?= $(shell xcrun --show-sdk-path 2>/dev/null)
    ifneq ($(SDKROOT),)
      CXXFLAGS += -isysroot $(SDKROOT)
      LDFLAGS  += -isysroot $(SDKROOT)
    endif
  endif
endif

CXXFLAGS += $(EXTRA_CXXFLAGS)

# OpenMP: Apple clang needs libomp (Homebrew) and -Xpreprocessor; GCC / LLVM clang use -fopenmp.
ifeq ($(OPENMP),1)
  ifeq ($(UNAME_S),Darwin)
    ifneq ($(IS_CLANG),0)
      LIBOMP_PREFIX ?= $(shell brew --prefix libomp 2>/dev/null || echo /opt/homebrew/opt/libomp)
      CXXFLAGS += -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
      LDFLAGS  += -L$(LIBOMP_PREFIX)/lib
      LDLIBS   += -lomp
    else
      CXXFLAGS += -fopenmp
      LDFLAGS  += -fopenmp
    endif
  else
    CXXFLAGS += -fopenmp
    LDFLAGS  += -fopenmp
  endif
endif

ifeq ($(ZLIB),1)
  CXXFLAGS += -DHAPCOAL_HAVE_ZLIB
  LDLIBS   += -lz
endif

# ------------------------------------------------------------------ sources ---
# Library: public headers in include/hapcoal, implementation and private headers in src/.
LIB_SRC  := $(wildcard src/*.cpp)
LIB_OBJ  := $(patsubst src/%.cpp,build/lib/%.o,$(LIB_SRC))
LIB      := lib/libhapcoal.a

CLI_OBJ  := build/tools/hapcoal_cli.o
EX_OBJ   := build/examples/main.o
TEST_OBJ := build/tests/test_main.o

BIN      := bin/hapcoal$(EXE)
EXAMPLE  := bin/hapcoal_example$(EXE)
TESTBIN  := bin/hapcoal_tests$(EXE)

.PHONY: all lib test regression doc install print-flags data chr22 bench clean rclean
all: $(LIB) $(BIN) $(EXAMPLE)
lib: $(LIB)

# static library
$(LIB): $(LIB_OBJ)
	@mkdir -p $(dir $@)
	rm -f $@
	$(AR) rcs $@ $^

# library objects may use private headers (src/)
build/lib/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -Iinclude -Isrc -c $< -o $@

# clients (CLI, example) see only the public headers
build/tools/%.o: tools/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -Iinclude -c $< -o $@
build/examples/%.o: examples/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -Iinclude -c $< -o $@
# tests also exercise internals
build/tests/%.o: tests/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -Iinclude -Isrc -c $< -o $@

$(BIN): $(CLI_OBJ) $(LIB)
	@mkdir -p $(dir $@)
	$(CXX) $(LDFLAGS) -o $@ $(CLI_OBJ) $(LIB) $(LDLIBS)
$(EXAMPLE): $(EX_OBJ) $(LIB)
	@mkdir -p $(dir $@)
	$(CXX) $(LDFLAGS) -o $@ $(EX_OBJ) $(LIB) $(LDLIBS)
$(TESTBIN): $(TEST_OBJ) $(LIB)
	@mkdir -p $(dir $@)
	$(CXX) $(LDFLAGS) -o $@ $(TEST_OBJ) $(LIB) $(LDLIBS)

test: $(TESTBIN)
	./$(TESTBIN)

regression: $(BIN)
	sh tests/regression.sh check

doc:
	@command -v doxygen >/dev/null 2>&1 || { echo "doxygen not found (e.g. brew install doxygen / apt install doxygen)"; exit 1; }
	doxygen Doxyfile
	@echo "API documentation: doc/html/index.html"

install: $(LIB)
	mkdir -p $(PREFIX)/include/hapcoal $(PREFIX)/lib
	cp include/hapcoal/*.hpp $(PREFIX)/include/hapcoal/
	cp $(LIB) $(PREFIX)/lib/

# Flags for a client program that links lib/libhapcoal.a
print-flags:
	@echo "CXXFLAGS: -std=c++17 -I$(CURDIR)/include"
	@echo "LIBS:     $(CURDIR)/$(LIB) $(strip $(LDFLAGS) $(LDLIBS))"

data:
	sh scripts/fetch_chr22.sh

# Human chr22 / 1000 Genomes CEU validation (downloads data on first use).
chr22: $(BIN)
	sh scripts/run_chr22.sh

bench: $(BIN)
	sh scripts/bench.sh

clean:
	rm -rf build bin lib doc

# Results only: everything under out/ (default output location of all runs).
rclean:
	rm -rf out
	mkdir -p out

-include $(LIB_OBJ:.o=.d) $(CLI_OBJ:.o=.d) $(EX_OBJ:.o=.d) $(TEST_OBJ:.o=.d)
