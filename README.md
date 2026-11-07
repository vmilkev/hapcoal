# hapcoal

hapcoal is a C++17 static library, with a command-line tool, that simulates haplotypes using the SMC, SMC′ and MaCS approximations of the coalescent with recombination. It supports population size histories, genetic maps and individuals of any ploidy. It returns packed haplotypes, genetic maps (true, realised and LD-inferred) and LD decay curves or matrices, either in memory, per chromosome, or as files. It can also fit a population size history so that the simulated LD matches real data.

**Full documentation** (methods, implementation, architecture, usage, file formats, validation and references): [docs/hapcoal.md](docs/hapcoal.md). API reference: `make doc`, which writes `doc/html/index.html`.

## Build

```sh
make                 # lib/libhapcoal.a, bin/hapcoal (CLI), bin/hapcoal_example (mock client)
make test            # unit + statistical tests
make regression      # all outputs identical to the recorded reference (needs `make data`)
make data            # download 1000 Genomes chr22 + HapMap II map into data/
make chr22           # human chr22 validation test case
make doc             # Doxygen HTML API documentation
make print-flags     # flags to compile / link a client against the library
make install PREFIX=/path
make clean           # remove build files (build/, bin/, lib/, doc/)
make rclean          # remove all results in out/ (downloaded data/ is kept)
```

Build options: `CXX=g++-15`, `OPENMP=0`, `ZLIB=0`, `NATIVE=1`, `DEBUG=1`. Requirements: a C++17 compiler; optionally OpenMP (on macOS: `brew install libomp`), zlib and Doxygen.

## Using the library

```cpp
#include <hapcoal/hapcoal.hpp>

hapcoal::SimParams p;
p.model = hapcoal::Model::MaCS;
p.nInd = 1000;  p.ploidy = 2;
p.chromosomes = {{"1", 0, 50e6, 1.2e-8}};
hapcoal::SimulationResult r = hapcoal::CoalescentSimulator(p).run();
int a = r.chromosomes[0].haplotypes.allele(0, 0);
```

`examples/main.cpp` is a complete mock client: a host-project class that uses in-memory results, per-chromosome streaming and file output.

## Command-line tool

```sh
bin/hapcoal simulate --model smcp --n-ind 500 --chr 1:50000000 --hotspots --out out/run1
bin/hapcoal ld --vcf <phased.vcf.gz> --panel <panel> --pop CEU --out out/obs
bin/hapcoal chr22
python3 scripts/plot_ld.py out/chr22      # LD decay and genetic-map plots
```

## Layout

```
include/hapcoal/  public headers (the library interface)
src/              implementation and internal headers
tools/            command-line tool
examples/         mock client (main.cpp)
tests/            unit tests, output regression test (golden checksums)
scripts/          data download, chr22 run, plotting, benchmark
docs/             project documentation, implementation plan
legacy/           previous implementation (reference only)
```

All results are written under `out/` by default.
