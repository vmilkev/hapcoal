# hapcoal — Sequentially Markov coalescent simulation of haplotypes

**Version 1.0.0** · C++17 static library with OpenMP · SMC, SMC′ and MaCS

hapcoal simulates haplotypes of a sample from a population with a given size history and genetic map. It uses three sequentially Markov approximations of the coalescent with recombination: SMC (McVean & Cardin 2005), SMC′ (Marjoram & Wall 2006) and MaCS (Chen, Marjoram & Wall 2009).

It is a library meant to be called from the classes of a larger program. The `CoalescentSimulator` class returns the following, in memory, per chromosome through a callback, and/or as files:
- packed haplotypes of individuals of any ploidy;
- the true genetic map, a map re-estimated from the simulated recombination events, and a map inferred from linkage disequilibrium (LD);
- LD decay curves and LD matrices.

The library also reads phased VCF files, so observed LD can be computed with the same code. It can fit a population size history to an observed LD decay curve, so that the simulated haplotypes match real data.

---

## Contents

1. Scope and features
2. Quick start
3. Methods
4. Implementation
5. Software architecture
6. Library usage
7. Command-line tool
8. File formats
9. Verification and validation
10. Performance
11. Assumptions and limitations
12. Build, portability and dependencies
13. Version history
14. References

---

## 1. Scope and features

| Area | Capability |
|---|---|
| Models | SMC, SMC′, MaCS (retained-ARG window *h*, given in bp) |
| Population | Reference size N0, piecewise-constant N(t) (file, epoch list or built-in Tennessen 2012 CEU preset), any ploidy p |
| Genome | Any number of chromosomes or regions; uniform rate, simulated hotspot landscape, or user genetic map (HapMap, PLINK, 1000 Genomes/IMPUTE formats, plain or gzip) |
| Output | In-memory haplotype matrices, a per-chromosome callback (bounded memory), and/or files: packed binary haplotypes with index files, optional 0/1 text |
| Maps | True, realised (from simulated recombination events) and LD-inferred (Hill–Weir fit) |
| LD | Haplotype or genotype r², LD decay curve, sparse and dense LD matrices |
| Real data | Phased VCF reader, observed LD decay, fit of a demography to an observed LD decay |
| Performance | O(log n) tree operations, bit-parallel genotypes and LD, OpenMP parallelism, deterministic results for any thread count |

---

## 2. Quick start

```sh
make                  # lib/libhapcoal.a, bin/hapcoal (CLI), bin/hapcoal_example (mock client)
make test             # unit and statistical tests
./bin/hapcoal_example  # library usage demonstration
make print-flags      # flags for compiling / linking a client
```

```cpp
#include <hapcoal/hapcoal.hpp>

hapcoal::SimParams p;
p.model  = hapcoal::Model::MaCS;
p.nInd   = 1000;                         // 1000 diploid individuals
p.ploidy = 2;
p.Ne     = 10000;  p.mu = 1.25e-8;
p.mapFile = "genetic_map_GRCh37_chr22.txt";   // simulate the mapped range of chr22
hapcoal::SimulationResult r = hapcoal::CoalescentSimulator(p).run();
int allele = r.chromosomes[0].haplotypes.allele(/*haplotype*/ 0, /*SNP*/ 0);
```

---

## 3. Methods

### 3.1 Notation and scaling

| Symbol | Meaning |
|---|---|
| n | number of sampled haplotypes; for N_ind individuals of ploidy p, n = p·N_ind |
| N0 | reference effective population size (`SimParams::Ne`) |
| N(g) | population size g generations ago; λ(t) = N/N0 |
| t | time backwards in **coalescent units** of p·N0 generations |
| r(x) | recombination rate (Morgans per bp per generation) at position x |
| G(x) | genetic map, G(x) = ∫ r(y) dy (Morgans) |
| μ | mutation rate per bp per generation |
| τ_x | local genealogy (marginal tree) at position x |
| T(x) | total branch length of τ_x (coalescent units) |
| K(t) | number of branches of τ_x crossing time t |

The population has p·N(g) gene copies. Under the Wright–Fisher model and its diffusion limit, a pair of lineages coalesces at rate 1/(p·N) per generation. In coalescent units this is 1/λ(t), so k lineages coalesce at total rate C(k,2)/λ(t) (Kingman 1982). One coalescent unit of branch length corresponds to p·N0 generations. Hence a branch of length ℓ experiences crossovers at rate p·N0·ℓ per Morgan and mutations at rate μ·p·N0·ℓ per bp. For p = 2 this gives the usual scaled rates θ = 4N0μ and ρ = 4N0r per bp.

### 3.2 Population size history

N(t) is piecewise constant with epochs [s_j, s_{j+1}) and relative sizes λ_j. For a process with hazard c/λ(t) starting at t0, the waiting time T satisfies

  P(T > t) = exp(−H(t0, t)),  H(t0, t) = ∫_{t0}^{t} c/λ(s) ds,

where H is piecewise linear. T is sampled by inversion: draw E ~ Exp(1) and solve H(t0, T) = E, one epoch at a time. Exponential growth phases (presets and fitted models) are discretised into epochs, using the size at each epoch's midpoint.

### 3.3 The coalescent with recombination and its sequential approximations

The ancestry of a recombining sequence is described exactly by the ancestral recombination graph (ARG; Hudson 1983; Griffiths & Marjoram 1997). Wiuf & Hein (1999) showed that the sequence of local trees τ_x can be generated by moving along the chromosome, but the process is not Markov in x.

The sequentially Markov coalescent (SMC) family approximates it by a Markov process along the sequence, using the following steps.

1. **Initial tree.** τ_{x0} is a Kingman coalescent tree under N(t): with k lineages, wait for hazard C(k,2)/λ(t), then merge a uniformly chosen pair.
2. **Distance to the next recombination.** Given τ_x, crossovers occur along the sequence at intensity p·N0·T(x) per Morgan. The genetic distance to the next event is Δg ~ Exp(p·N0·T(x)), and the next breakpoint is x′ = G⁻¹(G(x) + Δg). Sampling on the genetic scale is exact for any piecewise-constant r(x), including narrow hotspots; zero-rate segments are skipped.
3. **Mutations** (infinite sites). On [x, x′) the number of mutations is Poisson(μ·p·N0·T(x)·(x′ − x)). Positions are uniform on the interval. Each mutation falls on a branch chosen with probability proportional to its length, and the haplotypes below that branch carry the derived allele.
4. **Recombination point.** A point is chosen uniformly on τ_x: branch b with probability ℓ_b/T, and time t_r ~ U(t_b, t_{pa(b)}). The lineage above t_r is detached and floats upwards from t_r. Let op = pa(b) be the old parent and s the sibling of b.
5. **Re-coalescence** of the floating lineage, as defined by the model (below), gives τ_{x′}.

**SMC** (McVean & Cardin 2005). The old branch of b above t_r is removed and op disappears, so s is joined directly to the old grandparent. This leaves the pruned tree τ⁻. The floating lineage coalesces with the k⁻(t) branches of τ⁻ crossing t, at hazard k⁻(t)/λ(t), and attaches to one of them chosen uniformly. In terms of the original tree,

  k⁻(t) = K(t) − 1{t < t_op},

so the hazard needs no pruned copy: start with K(t_r) − 1 and skip the decrement at node op. Above the root k⁻ = 1. If the chosen branch is op's branch (t ≥ t_op), the lineage attaches to s's branch, which continues through op's former position in τ⁻.

**SMC′** (Marjoram & Wall 2006). As SMC, but the floating lineage may also coalesce with its own old branch, from t_r up to t_op. The number of available lineages is therefore

  k′(t) = k⁻(t) + 1{t < t_op} = K(t),

i.e. exactly the branches of the original tree crossing t. Coalescing with its own branch leaves the local tree unchanged; this is an *invisible* recombination. Such events are counted, and they lengthen tree spans exactly as under the ARG. SMC′ reproduces many properties of the ARG closely (Wilton et al. 2015); see §9.2 for the two-locus covariance.

**MaCS** (Chen, Marjoram & Wall 2009). MaCS stores the part of the ARG used by the local trees of the last h bp. Edges that have not been ancestral at any position in (x − h, x] are discarded. Recombination points are placed on the local tree only, as above. The floating lineage, however, coalesces with any retained edge: at time t the hazard is k_G(t)/λ(t), where k_G(t) is the number of retained edges crossing t, and the target edge is uniform among them.

If the target edge is not ancestral at x′, the transmitted material follows that edge upwards through the retained graph. At a crossover node it takes the parent edge created most recently, which carries the material of positions to the right. The material joins the local tree at the first node where it meets ancestral material of the other sampled lineages. The new local tree is the remaining tree with the subtree of b attached at that meeting point.

When h = 0, the retained graph at an event is the local tree plus the old branch of b, so MaCS equals SMC′. hapcoal reproduces this bit for bit. As h grows, MaCS approaches the full ARG.

### 3.4 Genetic maps

A map is a set of breakpoints x_i with cumulative positions G_i, linearly interpolated, so r is piecewise constant. The input formats are HapMap II (`chr pos rate cM`), PLINK (`chr id cM pos`) and 1000 Genomes/IMPUTE (`pos rate cM`).

Without an input map, either a uniform map or a **hotspot landscape** is generated. Hotspot centres are uniform on the chromosome outside an optional centromere, at a density of 20 per Mb by default. The relative rate is

  r(x) ∝ 1 + Σ_h (A_h − 1)·exp(−(x − c_h)²/(2 s_h²)),

with peak multipliers A_h ~ U(10, 20) and widths s_h ~ U(800, 1500) bp, typical of human hotspots (Myers et al. 2005). The rate in the centromere is multiplied by 0.1. Finally the map is rescaled so that its mean rate equals the requested rate.

### 3.5 Linkage disequilibrium

**Haplotype r²** (Hill & Robertson 1968). For SNPs i, j with derived-allele frequencies p_i, p_j and joint frequency p_ij,

  D_ij = p_ij − p_i p_j,  r²_ij = D_ij² / (p_i(1 − p_i) p_j(1 − p_j)).

**Genotype r².** The squared Pearson correlation between diploid dosages g = h_{2k} + h_{2k+1} over individuals k. It does not depend on phase, and is used to rule out phasing (switch) errors in real data as a source of discrepancy.

**LD decay.** SNPs with minor-allele frequency ≥ 0.05 (default) are kept. If necessary they are thinned at random to at most 200,000 per chromosome. Pairs within 500 kb are averaged in 1-kb distance bins. The same code processes simulated haplotypes and VCF data.

Under the coalescent, E[r²] ≈ 1/(1 + C) for C = 4N_e·c (Sved 1971). The more accurate finite-sample expectation of Hill & Weir (1988) is

  E[r²] ≈ (10 + C)/((2 + C)(11 + C)) · [1 + (3 + C)(12 + 12C + C²)/(n(2 + C)(11 + C))],

which is the basis of the LD-based map (§3.6).

### 3.6 Genetic maps from simulated data

*Realised map.* Conditional on the sequence of local trees, recombination breakpoints form a Poisson process with intensity r(x)·p·N0·T(x) per bp. For a window W with constant r, the estimator

  r̂_W = N_W / (p·N0·∫_W T(x) dx),

where N_W is the number of events in W, is unbiased, and its variance is that of a Poisson count. Invisible events are included, since they are events of the same process.

*LD-based map.* This map is inferred from the haplotypes alone. In each window (1 Mb by default), pairs up to 100 kb apart are binned by distance (2-kb bins). The population-scaled rate ρ = 2pN_e·c per bp is fitted by weighted least squares of the binned mean r² against the Hill–Weir expectation with C = ρd, using a grid search on log₁₀ρ followed by golden-section refinement. c is converted to cM/Mb with the supplied N_e. The shape of the map is therefore estimated without bias, but the absolute scale is only correct when N_e describes the recent population size.

### 3.7 Fitting a population size history to observed LD

LD at distance d mainly reflects coalescence times of order 1/(2c(d)) generations (Hayes et al. 2003). Short-range LD therefore informs old population sizes, and long-range LD recent ones. hapcoal fits a six-parameter **bottleneck + growth** model:

  N(g) = N0·(N1/N0)^{g/T1} for g < T1;  N(g) = N_B for T1 ≤ g < T1 + ΔT;  N(g) = N_A afterwards.

The model is fitted by minimising

  f(θ) = [ (1/K) Σ_k log²( r̄²_sim,k(θ) / r̄²_obs,k ) ]^{1/2}

over K = 9 distance classes (edges at 0, 1, 2, 5, 10, 20, 50, 100, 200 and 500 kb), so every class carries the same weight. Minimisation uses the Nelder–Mead simplex method on log θ (Nelder & Mead 1965). Each evaluation is a full simulation of the observed region with the observed sample size, using a fixed seed (common random numbers), which makes f smooth enough for a derivative-free optimiser. The simplex starts near the Tennessen et al. (2012) European model.

---

## 4. Implementation

### 4.1 Local tree (`src/LocalTree.*`)

- **Storage.** Nodes are kept in structure-of-arrays form (parent, left, right, time). Samples are 0…n−1 and internal nodes n…2n−2. A prune-and-regraft reuses the removed node's id, so the tree always has exactly 2n − 1 nodes.
- **Fenwick tree** (binary indexed tree; Fenwick 1994) over branch lengths. Choosing a point proportional to length (for recombination and mutation) is an O(log n) descent. Updates are O(log n), and the tree is rebuilt periodically to remove floating-point drift.
- **TimeIndex** (`src/TimeIndex.hpp`). An ordered multiset of node times with their changes to the lineage count (+1 per sample and crossover, −1 per coalescence), stored in sorted blocks of about 256 with block sums. K(t) is a prefix sum, and the coalescence hazard is integrated by walking forwards from t_r through node times and epoch boundaries.
- **Branch crossing t.** If many branches cross t (K·8 ≥ 2n − 1), a branch is found by rejection sampling over all nodes, O(1) expected. Otherwise the K crossing branches are enumerated by descending from the root, O(K).
- **Leaf-set cache.** Each node can carry a bitset of the samples below it, as a lazy cache invalidated along root paths at each regraft. The genotype row of a mutation is then a copy of ⌈n/64⌉ words instead of a subtree traversal. The cache is enabled when it fits in 512 MB (n ≲ 45,000), and falls back to traversal otherwise.

### 4.2 MaCS retained ARG (`src/MacsHistory.*`)

- **Nodes and edges.**
  - Nodes are samples, coalescences (two child edges, one parent edge) and crossovers (one child edge, two parent edges: older and newer).
  - The open-ended edge above the grand most-recent common ancestor (MRCA) has no top node.
  - Each edge stores its *owner*, the local-tree branch it currently belongs to (−1 if non-ancestral), and the last position at which it was ancestral.
- **Recombination.** A crossover node is inserted on b's path at t_r. b's former path up to op becomes non-ancestral, and op's path is reassigned to s. The coalescence time is sampled with K_ARG(t) from the ARG's TimeIndex.
- **Target edge.** The target is chosen uniformly by stratified sampling. With probability K_tree(t)/K_ARG(t) a local-tree branch is drawn by the fast LocalTree sampler and mapped to its ARG edge. Otherwise an edge is drawn from a dense list of non-ancestral edges, excluding b's replaced path, which is already represented by b's branch.
- **Following material.** From a non-ancestral target the material follows the newest-edge path upwards until it meets another owner. Ownership is then updated along the affected paths only.
- **Pruning.** Non-ancestral edges enter a min-heap keyed by last use. Edges older than x − h are removed. Nodes left with one child and one parent edge are fused, and dead branches are removed in cascade.
- **Validation.** In `DEBUG=1` builds the full consistency of the ARG and the local tree is checked after every event.

### 4.3 Random numbers and reproducibility

The standard `std::*_distribution` classes produce implementation-defined sequences, so they are not used. hapcoal draws raw numbers from `std::mt19937_64` (Matsumoto & Nishimura 1998) and implements its own uniform, exponential (inversion), Poisson (Hörmann 1993 PTRS) and unbiased integer distributions. Chromosome i uses an independent stream seeded with `std::seed_seq{seed, i + 1}`. Results therefore depend only on the parameters and seed, not on the number of threads, platform or compiler.

### 4.4 Haplotype storage and output (`src/HaplotypeStore.*`)

During simulation every mutation is appended as a SNP-major row of n bits to a `SnpBuffer`. The buffer keeps rows in memory, or streams them to a temporary file in 16 MB blocks to keep memory bounded. Afterwards the rows are transposed in 64 × 64-bit tiles (Warren 2012) into haplotype-major form:
- in memory (`HaplotypeMatrix` per chromosome), or
- into `PREFIX.hap.bin`, in chunks of about 2 GB, with parallel scattered writes into a pre-allocated file.

For LD calculations, only the selected (common) rows are read back.

### 4.5 Parallelism and vectorisation

- **Chromosomes** are simulated in parallel with OpenMP, longest first.
- **LD** is computed in parallel over the first SNP of each pair, with per-thread bins.
- **Transposition and output** are parallel over haplotype blocks, with one file handle per thread.
- **Inner loops vectorise:** AND + popcount over 64-bit words in LD, OR over words in the leaf-set cache, and word copies.

### 4.6 Complexity per recombination event (n haplotypes)

| Operation | Cost |
|---|---|
| Choose a point on the tree | O(log n) |
| Sample the coalescence time | O(number of node times or epochs passed) + O(n/B + B) for K(t_r) |
| Choose the target lineage | O(1) expected, or O(K(t)) |
| Regraft | O(log n) Fenwick updates + O(B) time-index updates + O(depth) cache invalidation |
| MaCS extra work | O(path lengths) ownership updates, amortised O(1) pruning |
| Per mutation | O(log n) + O(n/64) (genotype row) |

---

## 5. Software architecture

### 5.1 Layout

```
include/hapcoal/        PUBLIC API (installed; clients include <hapcoal/hapcoal.hpp>)
  hapcoal.hpp             umbrella header, version
  CoalescentSimulator.hpp  SimParams, ChromosomeSpec, Model, results, CoalescentSimulator
  Demography.hpp         piecewise-constant N(t), Tennessen 2012 preset
  GeneticMap.hpp         maps: files, hotspots, bp <-> M
  Haplotypes.hpp         BitRows, HaplotypeMatrix, .hap.bin reader
  LDCalculator.hpp       r2, LD decay, matrices, decay utilities
  MapEstimator.hpp       realised and LD-based maps
  VcfReader.hpp          phased VCF, observed LD
  DemographyFit.hpp      fit of a demography to an observed LD decay
src/                   IMPLEMENTATION (private headers are not installed)
  LocalTree.*  MacsHistory.*  TimeIndex.hpp  Rng.hpp  HaplotypeStore.*  + *.cpp of the public API
tools/hapcoal_cli.cpp   command-line tool (public API only)
examples/main.cpp      mock client class of a host project (public API only)
tests/                 unit/statistical tests, output regression test
scripts/               data download, chr22 validation, plotting, benchmark
docs/                  this document, implementation plan
```

### 5.2 Components and dependencies

```
              host project / tools / examples
                          │  (#include <hapcoal/hapcoal.hpp>, link libhapcoal.a)
                          ▼
  ┌──────────────────── CoalescentSimulator ───────────────────┐
  │ SimParams → run(): maps → simulate (parallel) → outputs     │
  │            simulateChromosome(): SMC / SMC' / MaCS           │
  └──┬──────────┬───────────┬────────────┬───────────┬─────────┘
     │          │           │            │           │
  LocalTree  MacsHistory  Demography  GeneticMap  HaplotypeStore ── Haplotypes
     │          │                                    │
  TimeIndex, Rng                                LDCalculator ── MapEstimator
                                                     ▲
                       VcfReader (observed LD) ──────┤
                       DemographyFit ── uses CoalescentSimulator + LDCalculator
```

Dependency rules:
- Public headers include only public headers.
- Internal types (`LocalTree`, `SnpBuffer`, `Rng`) appear in public headers only as forward declarations, on the advanced `simulateChromosome()` interface and on `Demography::sampleWaiting()`.
- The tool and the example compile with `-Iinclude` only, which is checked by the build.

### 5.3 Class responsibilities

| Class / module | Responsibility |
|---|---|
| `CoalescentSimulator` | The single class implementing the three models (`simulateChromosome`) and the library interface (`run`, results, files) |
| `SimParams`, `ChromosomeSpec` | All inputs, with defaults; plain data |
| `SimulationResult`, `ChromosomeResult`, `ChromStats` | All outputs; plain data plus small helpers (`genomeRow`, `positionCM`) |
| `LocalTree` | Local genealogy and its O(log n) operations |
| `MacsHistory` | Retained ARG of MaCS |
| `Demography` | N(t), waiting times, presets |
| `GeneticMap` | Maps, hotspot landscape, conversion |
| `SnpBuffer`, `assembleHaplotypes` | Row storage, transposition, file writing |
| `LDCalculator`, `MapEstimator` | r², decay, matrices; realised and LD-based maps |
| `VcfReader`, `observedLD` | Real data |
| `fitBottleneckGrowth` | Demography inference from LD |

### 5.4 Data flow of `run()`

1. Validate parameters and load the demography and maps (constructor).
2. Build the genetic map of each chromosome.
3. Simulate chromosomes in parallel into `SnpBuffer`s.
4. If `outPrefix` is set, write `hap.bin`, the index files and the true and realised maps.
5. For each chromosome in input order:
   1. compute LD (decay, matrices, LD-based map);
   2. build the `HaplotypeMatrix`;
   3. release the buffer;
   4. call the callback;
   5. append the chromosome to the `SimulationResult`.
6. Write the run log and remove temporary files.

### 5.5 Error handling, threading and logging

- **Errors** are reported by exceptions: `std::invalid_argument` for parameters, `std::runtime_error` for I/O, and `std::logic_error` for internal inconsistencies (never expected). Exceptions thrown in parallel regions are collected and re-thrown on the calling thread.
- **Threading.** `run()` must not be called concurrently on the *same* object. Separate `CoalescentSimulator` objects are independent. The callback is invoked on the calling thread, in chromosome order.
- **Logging.** The library is silent unless `SimParams::log` points to a stream.

### 5.6 Architecture review for the library version

The restructuring into a library was limited to changes that were needed:

| Change | Reason |
|---|---|
| Public `include/hapcoal/` vs private `src/` | Small, stable interface for the host project |
| Command-line parsing moved out of the simulator into `tools/` | A library should not parse command lines |
| Demography fit, Tennessen preset and observed LD moved into the library | They were previously only in `main.cpp`, but are reusable features |
| `run()` returns results, with an optional callback; files optional | Integration into other classes |
| `SnpBuffer` (memory or temporary file) | In-memory use without files; bounded memory for large runs |
| Library silent by default | No unsolicited output |

Not changed, deliberately:
- The single simulator class for the three models (an earlier requirement).
- The separation of the LD and map modules (distinct responsibilities).
- Class names (stability).

All outputs were verified to be bit-identical before and after the restructuring (§9.3).

---

## 6. Library usage

### 6.1 Linking

```sh
make lib
make print-flags
# CXXFLAGS: -std=c++17 -I<hapcoal>/include
# LIBS:     <hapcoal>/lib/libhapcoal.a -L/opt/homebrew/opt/libomp/lib -lomp -lz   (macOS; Linux: -fopenmp -lz)
make install PREFIX=$HOME/.local     # optional: headers + library
```

### 6.2 Receiving results

```cpp
hapcoal::SimParams p;   // see CoalescentSimulator.hpp for every field
p.nInd = 500; p.ploidy = 2; p.model = hapcoal::Model::SMCPrime;
p.chromosomes = {{"1", 0, 50e6, 1.2e-8}, {"2", 0, 40e6, 1.2e-8}};
p.epochs = hapcoal::Demography::tennessen2012CEU();   // or {{0, 20000}, {500, 2000}, ...}

// (a) everything in memory
hapcoal::SimulationResult res = hapcoal::CoalescentSimulator(p).run();
const hapcoal::ChromosomeResult& c1 = res.chromosomes[0];
c1.haplotypes.allele(h, j);   // allele of haplotype h at SNP j
c1.positions[j];              // bp
c1.positionCM(j);             // cM
res.genomeRow(h);             // all chromosomes of haplotype h, consecutive

// (b) streaming: one chromosome at a time, rows in temporary files (bounded memory)
p.tempDir = "/scratch/hapcoal_tmp";
hapcoal::CoalescentSimulator(p).run([&](hapcoal::ChromosomeResult& c) {
    myGenome.add(c.name, std::move(c.haplotypes), c.positions);
});

// (c) files as written by the CLI
p.outPrefix = "out/run1";  p.keepHaplotypes = false;
hapcoal::CoalescentSimulator(p).run();
```

### 6.3 Matching real data

```cpp
hapcoal::LDOptions ld;  ld.threads = 8;
auto ceu = hapcoal::readPanelSamples("integrated_call_samples_v3.20130502.ALL.panel", "CEU");
hapcoal::ObservedLD obs = hapcoal::observedLD("ALL.chr22....vcf.gz", ceu, 16051347, 51229805, ld);

hapcoal::SimParams base;  base.model = hapcoal::Model::SMCPrime;  base.nInd = int(obs.nHap / 2);
base.mapFile = "genetic_map_GRCh37_chr22.txt";  base.mu = 2.36e-8;  base.seed = 2026;
hapcoal::DemographyFitResult fit = hapcoal::fitBottleneckGrowth(obs.decay, base, 120, &std::cerr);
base.epochs = fit.epochs;     // simulate haplotypes whose LD matches the observations
```

`examples/main.cpp` is a complete mock client: a class `FounderPopulation` of a host project that uses all three ways of receiving results.

---

## 7. Command-line tool

```sh
bin/hapcoal simulate --model macs --n-ind 5000 --chromosomes data/grch37_autosomes.tsv --out out/genome
bin/hapcoal simulate --model smcp --n-ind 99 --map data/genetic_map_GRCh37_chr22.txt \
                    --demography out/chr22/fitted_demography.tsv --mu 2.36e-8 --estimate-map --out out/c22
bin/hapcoal ld --vcf <vcf.gz> --panel <panel> --pop CEU --out out/obs_CEU
bin/hapcoal chr22                      # validation test case (make chr22)
python3 scripts/plot_ld.py out/chr22  # plots: LD decay, genetic maps
```

`bin/hapcoal simulate --help` lists all options. They map one-to-one onto `SimParams`.

---

## 8. File formats

| File | Content |
|---|---|
| `PREFIX.hap.bin` | 64-byte header (`HAPCOAL1`, version, n_hap, ploidy, n_chrom, n_snps, words_per_row, header_bytes; little-endian uint64), then n_hap rows of words_per_row uint64. Each row holds all chromosomes consecutively, each padded to 64 bits. Bit j of a chromosome is bit j mod 64 of word `row_word_offset + ⌊j/64⌋`. Haplotypes of individual i are rows i·p … i·p + p − 1. |
| `PREFIX.chrom.tsv` | chrom, start_bp, end_bp, length_cM, n_snps, first_snp_index, row_word_offset, n_words |
| `PREFIX.snps.tsv` | index, chrom, pos_bp (1-based), pos_cM, alt_freq |
| `PREFIX.samples.tsv` | hap_index, individual, haplotype |
| `PREFIX.hap.txt` | optional; one 0/1 line per haplotype |
| `PREFIX.map.true.tsv` | chrom, pos_bp, rate_cM_per_Mb, pos_cM |
| `PREFIX.map.realised.tsv` | chrom, start_bp, end_bp, events, rate_cM_per_Mb, pos_cM |
| `PREFIX.map.ld.tsv` | chrom, start_bp, end_bp, pairs, rho_per_kb, rate_cM_per_Mb, pos_cM |
| `PREFIX.ld_decay.tsv` | chrom, bin_start_bp, bin_mid_bp, mean_r2, n_pairs (`all` = pooled) |
| `PREFIX.ld_matrix.tsv` / `.ld_region.tsv` | sparse pairs (chrom, snp_i, snp_j, pos_i, pos_j, r2) / dense matrix |
| `PREFIX.log` | parameters, per-chromosome statistics, timings |
| demography file | lines `generations_ago  N` (`#` comments) |

---

## 9. Verification and validation

### 9.1 Unit and statistical tests (`make test`; 105 checks)

| Component | Test |
|---|---|
| RNG, demography | Exponential, Poisson and integer distributions; P(T₂ > t) under a bottleneck equals exp(−H) |
| Genetic maps | Interpolation and inversion, plateaus, the three file formats, hotspot rescaling, clipping |
| Initial tree | E[T_MRCA] = 2(1 − 1/n); E[T_total] = 2H_{n−1}; E[ξ_i] = 2/i (site frequency spectrum) |
| SMC / SMC′ / MaCS | Sequence-averaged T_total and T_MRCA equal the Kingman values (stationarity); E[the number of events] = r·p·N0·L·E[T_total]; validity of every tree |
| Two-locus covariance | n = 2: Cov(T(0), T(d)) against an exact chain for SMC and SMC′, and against the ARG value (ρ + 18)/(ρ² + 13ρ + 18) (Griffiths 1981) for MaCS with large h |
| MaCS ARG | Consistency of ARG and local tree after every one of 20,000 events; h = 0 graph equals the local tree |
| Leaf-set cache | Equals subtree traversal after 3,000 regrafts |
| Storage | `hap.bin`, text and in-memory matrices reproduce the input (file and memory buffers) |
| Library API | In-memory, file, streaming and 1-vs-4-thread results are identical; temporary files removed |
| LD | Haplotype and genotype r² equal a naive double-precision calculation |
| Maps | Realised map recovers a 1 → 5 cM/Mb step; the LD-based map detects it |
| VCF | SNP filtering, sample subsets, allele coding |

### 9.2 Two-locus covariance of pairwise coalescence times (n = 2, ρ = 4N0·r·d)

| ρ | SMC exact | SMC hapcoal | SMC′ exact | SMC′ hapcoal | ARG | MaCS, h = 1 Mb |
|---|---|---|---|---|---|---|
| 0.5 | 0.667 | 0.665 | 0.741 | 0.734 | 0.747 | 0.741 |
| 2 | 0.333 | 0.337 | 0.393 | 0.393 | 0.417 | 0.418 |
| 8 | 0.110 | 0.109 | 0.123 | 0.128 | 0.140 | 0.153 |

The SMC′ covariance is close to the ARG value but slightly lower. MaCS with a long history recovers the ARG value.

### 9.3 Output regression (`make regression`)

`tests/regression.sh` runs a fixed set of commands. They cover:
- all three models and two chromosomes;
- hotspot and real maps, demographies and ploidy 3;
- the MaCS history, genotype LD, sparse and dense LD matrices and text output;
- the VCF reader and the chr22 demography fit.

The SHA-256 checksums of all 52 outputs are compared with `tests/golden.sha256`. Any change of results, intended or not, is therefore detected. The library restructuring (version 1.0.0) passed this test with all outputs bit-identical.

The reference was recorded with Apple clang 21. With GCC 15, haplotypes, maps, the observed LD and the fitted demography are identical. Three text outputs differ only by floating-point rounding: GCC contracts multiply-adds by default (`-ffp-contract=fast`), so a SNP position can differ in its last bit. As a result, one SNP pair out of about 2·10⁸ falls into the neighbouring 1-kb LD bin, and one cM value differs in the 6th decimal. Regression checks should therefore be run with the compiler that recorded the reference.

### 9.4 Human chromosome 22 (`make chr22`)

Data: 1000 Genomes phase 3 CEU (99 individuals, 85,630 SNPs with MAF ≥ 0.05) and the HapMap II GRCh37 map (16.05–51.23 Mb, 74.1 cM). Simulations used the same sample size, map and MAF filter. The table gives the RMSE of mean r² over 1-kb bins from 0 to 500 kb:

| Demography | SMC | SMC′ | MaCS |
|---|---|---|---|
| Tennessen et al. 2012 (fitted to the SFS) | 0.044 | 0.047 | 0.037 |
| LD-fitted bottleneck + growth | 0.006 | 0.002 | 0.006 |

The Tennessen model produces about twice the observed LD at 20–200 kb. Genotype r² gives the same picture, so phasing errors in the data are not the cause. The fitted model has:
- ancestral N_A ≈ 11,900;
- a bottleneck N_B ≈ 3,000 between about 2,234 and 2,618 generations ago;
- growth from about 2,400 to 47,000.

With an independent seed, its RMSE is 0.004–0.008. Below about 2 kb the simulations slightly underestimate LD, because the HapMap map is linearly interpolated at that scale and smooths out narrow hotspots.

---

## 10. Performance

Measured on an Apple M1 Pro (10 cores) with Apple clang 21, `-O3` and OpenMP.

| Run | Simulation | Total | Peak memory |
|---|---|---|---|
| chr22 (35.2 Mb), 198 haplotypes | 0.2 s (SMC, SMC′), 0.7–2 s (MaCS) | < 2 s | < 50 MB |
| chr22, 10,000 haplotypes, SMC / SMC′ / MaCS | 1.7 / 1.7 / 3.9 s | 2.7 / 2.5 / 4.7 s | 0.56 GB |
| 22 autosomes (2.88 Gb), 10,000 haplotypes, MaCS (14.1M SNPs, 17.6 GB file) | 22 s | 124 s (I/O bound) | 2.8 GB |

---

## 11. Assumptions and limitations

- **Model assumptions.** Neutral, panmictic, Wright–Fisher population with piecewise-constant size. Sequential Markov approximations (§3.3). Infinite-sites mutation with biallelic SNPs. Crossovers only (no gene conversion). One population (no structure or migration).
- **SNP positions** are continuous internally. Files report ⌊x⌋ + 1, so two SNPs can rarely share a bp position.
- **LD-based map:** the absolute scale depends on the supplied N_e (§3.6).
- **Short-range LD** depends on the resolution of the input map (§9.4).
- **Memory.** In-memory results for genome-scale samples can be very large (10,000 haplotypes × 22 autosomes ≈ 17.6 GB). Use the callback with `tempDir`, or file output.
- **Leaf-set cache** is used up to n ≈ 45,000. Larger samples use subtree traversal, which is slower per mutation.

---

## 12. Build, portability and dependencies

| Item | Requirement |
|---|---|
| Language | ISO C++17. Compiler extensions are used only behind guards (`__builtin_popcountll`, `__restrict`). |
| Compilers tested | Apple clang 21 (macOS arm64) with Homebrew libomp; GCC 15 |
| Optional | OpenMP (`OPENMP=0` to disable), zlib for `.gz` input (`ZLIB=0`), Doxygen (`make doc`) |
| Build system | GNU make (also make 3.81 as shipped with macOS); `make DEBUG=1` adds assertions and full tree/ARG validation |
| Platforms | macOS, Linux; Windows via MSYS2/MinGW (`.exe` suffix, 64-bit file offsets via `_fseeki64`) |

---

## 13. Version history

| Version | Changes |
|---|---|
| 1.0.0 | Static library with public API (`include/hapcoal`) and in-memory or streaming results. Command-line tool and mock client built on the library. Demography fit, Tennessen preset and observed LD in the library. Doxygen documentation. Outputs bit-identical to 0.9. |
| 0.9 | Corrected rewrite of the legacy simulator: exact SMC/SMC′ transitions, MaCS with retained ARG, demography, maps, LD, VCF input, chr22 validation, performance work |

---

## 14. References

- 1000 Genomes Project Consortium (2015). A global reference for human genetic variation. *Nature* 526:68–74. doi:10.1038/nature15393
- Adrion JR et al. (2020). A community-maintained standard library of population genetic models. *eLife* 9:e54967. doi:10.7554/eLife.54967
- Chen GK, Marjoram P, Wall JD (2009). Fast and flexible simulation of DNA sequence data. *Genome Research* 19:136–142. doi:10.1101/gr.083634.108
- Fenwick PM (1994). A new data structure for cumulative frequency tables. *Software: Practice and Experience* 24:327–336. doi:10.1002/spe.4380240306
- Griffiths RC (1981). Neutral two-locus multiple allele models with recombination. *Theoretical Population Biology* 19:169–186. doi:10.1016/0040-5809(81)90016-2
- Griffiths RC, Marjoram P (1997). An ancestral recombination graph. In: Donnelly P, Tavaré S (eds), *Progress in Population Genetics and Human Evolution*, IMA Volumes in Mathematics and its Applications 87:257–270. Springer.
- Hayes BJ, Visscher PM, McPartlan HC, Goddard ME (2003). Novel multilocus measure of linkage disequilibrium to estimate past effective population size. *Genome Research* 13:635–643. doi:10.1101/gr.387103
- Hill WG, Robertson A (1968). Linkage disequilibrium in finite populations. *Theoretical and Applied Genetics* 38:226–231. doi:10.1007/BF01245622
- Hill WG, Weir BS (1988). Variances and covariances of squared linkage disequilibria in finite populations. *Theoretical Population Biology* 33:54–78. doi:10.1016/0040-5809(88)90004-4
- Hörmann W (1993). The transformed rejection method for generating Poisson random variables. *Insurance: Mathematics and Economics* 12:39–45. doi:10.1016/0167-6687(93)90997-4
- Hudson RR (1983). Properties of a neutral allele model with intragenic recombination. *Theoretical Population Biology* 23:183–201. doi:10.1016/0040-5809(83)90013-8
- International HapMap Consortium (2007). A second generation human haplotype map of over 3.1 million SNPs. *Nature* 449:851–861. doi:10.1038/nature06258
- Kingman JFC (1982). The coalescent. *Stochastic Processes and their Applications* 13:235–248. doi:10.1016/0304-4149(82)90011-4
- Marjoram P, Wall JD (2006). Fast "coalescent" simulation. *BMC Genetics* 7:16. doi:10.1186/1471-2156-7-16
- Matsumoto M, Nishimura T (1998). Mersenne twister: a 623-dimensionally equidistributed uniform pseudo-random number generator. *ACM Transactions on Modeling and Computer Simulation* 8:3–30. doi:10.1145/272991.272995
- McVean GAT, Cardin NJ (2005). Approximating the coalescent with recombination. *Philosophical Transactions of the Royal Society B* 360:1387–1393. doi:10.1098/rstb.2005.1673
- Myers S, Bottolo L, Freeman C, McVean G, Donnelly P (2005). A fine-scale map of recombination rates and hotspots across the human genome. *Science* 310:321–324. doi:10.1126/science.1117196
- Nelder JA, Mead R (1965). A simplex method for function minimization. *The Computer Journal* 7:308–313. doi:10.1093/comjnl/7.4.308
- Sved JA (1971). Linkage disequilibrium and homozygosity of chromosome segments in finite populations. *Theoretical Population Biology* 2:125–141. doi:10.1016/0040-5809(71)90011-6
- Tennessen JA et al. (2012). Evolution and functional impact of rare coding variation from deep sequencing of human exomes. *Science* 337:64–69. doi:10.1126/science.1219240
- Warren HS (2012). *Hacker's Delight*, 2nd ed. Addison-Wesley.
- Wilton PR, Carmi S, Hickerson MJ (2015). The SMC′ is a highly accurate approximation to the ancestral recombination graph. *Genetics* 200:343–355. doi:10.1534/genetics.114.173898
- Wiuf C, Hein J (1999). Recombination as a point process along sequences. *Theoretical Population Biology* 55:248–259. doi:10.1006/tpbi.1998.1403
