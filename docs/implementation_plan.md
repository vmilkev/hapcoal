# Plan: Coalescent haplotype simulator (SMC / SMC′ / MaCS)

## Context
`macs_simulator.cpp` is meant to simulate haplotypes along chromosomes using sequentially Markov coalescent models. Today it doesn't compile: the LD-binning refactor at L591-598 is half-done. Its SMC, SMC′ and MaCS logic is also incorrect (audit below). The goal is a correct, fast, memory-efficient C++17 (+OpenMP) simulator, with:
- one public class for the methods and the interface
- a test-case `main.cpp`
- a cross-platform Makefile
- validation against the observed LD decay of human chr22 (1000 Genomes CEU)

**User decisions:**
- Map output: true + realised + LD-based estimate.
- Demography: piecewise-constant N(t).
- Data: download 1000G chr22 + HapMap map.
- Haplotype output: packed binary + index files, with optional text 0/1.
- Target scale: 10k haplotypes over the whole genome.
- Observed LD: computed inside the C++ program (VCF reader).
- Old files: moved to `legacy/`.

**Environment:** macOS arm64, Apple clang 21, Homebrew libomp, g++-15, zlib in the SDK, 10 cores, 16 GB RAM. Python has no matplotlib or msprime yet.

---

## Stage 2/3 — Audit of the current implementation (problems and fixes)

### A. Model correctness (blocking)
1. **Coalescence time sampling is wrong in every model** (L465-491). It draws `Exp(k)` with k fixed at the recombination time, then walks the target up its ancestors. The number of lineages k(t) drops at every node of the remaining tree, and the rate also depends on N(t). **Fix:** step through the intervals between sorted node times (and demographic epochs), with a constant hazard k/λ(t) in each. Above the root, k = 1.
2. **The floating lineage can coalesce into its own subtree** (L433). `find_active_lineages_iterative(broken_node, …)` adds lineages from inside the detached subtree. `old_parent.time` has already been reset to 0 (L427), so the traversal goes down into it. This creates invalid or cyclic trees. **Fix:** only branches of the pruned tree are targets.
3. **SMC and SMC′ are identical in the code, and both are wrong** (L436-461).
   - **SMC (McVean & Cardin 2005):** remove the old branch above the recombination point. Coalesce with any branch of the pruned tree, including the merged sister + old-parent branch.
   - **SMC′ (Marjoram & Wall 2006):** same as SMC, but the lineage's own old branch [t_r, t_oldparent) also stays a target, which allows invisible recombinations.
   - The code instead removes the sister branch.
4. **This is not MaCS** (L445-484). MaCS (Chen, Marjoram & Wall 2009) keeps the *union of branches from recent local trees*, within a history window *along the sequence*. The new lineage may coalesce with any retained branch. The code instead caps the coalescence *time* at `t_r + 3.5`, which truncates the time distribution with no basis in MaCS. **Fix:** retained-branch store (Step 4). With h = 0 it reduces exactly to SMC′.
5. **Recombination positions ignore the genetic map** (L333-350). The waiting distance is drawn in bp, using the rate at the interval start only, so narrow hotspots are skipped. **Fix:** draw the waiting distance in genetic units, then convert with `genetic_to_physical` (L177). That function exists but is unused.

### B. Biological realism
6. N is constant; there is no size history (needed to fit CEU LD).
7. There are no individuals or ploidy, and time is scaled as if ploidy were always 2.
8. `generate_hotspot_map` uses a 40 kb grid on 200 Mb, so hotspots about 1 kb wide are missed. Rates are printed with `setprecision(5)`, which shows 1e-8 as `0.00000`.
9. `load_genetic_map` is a stub, so a user-supplied map is silently ignored.
10. `delta_x` has a 1 bp minimum (L352), which biases high-rate regions.
11. `std::*_distribution` output depends on the implementation, so results with the same seed differ across platforms.

### C. Performance and memory
12. `select_branch_by_cdf` scans all branches in O(n) for every mutation and recombination. **Fix:** Fenwick tree, O(log n).
13. Per-event vector allocations and O(n) traversals.
14. The full SNP×haplotype matrix stays in memory, and the text output uses 1 byte per allele: 1.1 GB for 1k haplotypes on 400 Mb, and about 15 GB for chr1 at 10k haplotypes. **Fix:** stream to disk.
15. Scalar popcount loop and no OpenMP. **Fix:** word-wise AND + popcount with OpenMP.
16. `plot_ld.py` reads columns the simulator no longer writes.

---

## Stage 8.1 — Design

### Scaling and units
- Haploid population size is p·N(t), where p = ploidy.
- Time is measured in units of p·N0 generations.
- λ(t) = N(t)/N0. Coalescence hazard is C(k,2)/λ(t) for k lineages; for the floating lineage it is k/λ(t).
- Recombination per bp per unit time is r(x)·p·N0. Mutation is μ·p·N0.
- Demography is given as epochs (generations, N), converted to coalescent units. The waiting time is drawn by inverting the cumulative hazard across epoch boundaries.

### Architecture (`src/`)
| File | Content |
|---|---|
| `CoalescentSimulator.hpp/.cpp` | **The single public class:** `Params` (model, n_ind, ploidy or n_hap, N0, demography, mu, rho, chromosome list, map path, MaCS history h, seed, threads, output options), `run()`, and the CLI parsing `fromArgs()`. It runs the SMC, SMC′ and MaCS transitions. |
| `LocalTree.hpp/.cpp` | SoA node arrays (parent, child[2], time); Fenwick tree over branch lengths; sorted internal-node times (k(t) = n − #internal nodes ≤ t); prune/regraft; initial tree under N(t) |
| `MacsHistory.hpp/.cpp` | Retained branches, each tagged with the last local-tree index where it existed; pruning at window h |
| `GeneticMap.hpp/.cpp` | HapMap format (`Chromosome Position(bp) Rate(cM/Mb) Map(cM)`) and the existing simple format; piecewise-constant rate segments; hotspot generator with explicit hotspot intervals (no fixed grid); bp↔cM conversion (reuses the logic from legacy L177-207) |
| `Demography.hpp` | Piecewise λ(t); integrated-hazard inversion |
| `Rng.hpp` | `mt19937_64` plus our own uniform, exponential and Poisson sampling, for identical results on every platform; one independent stream per chromosome via `seed_seq{seed, chr}` |
| `HaplotypeStore.hpp/.cpp` | SNP-major packed bitsets streamed to a temp file per chromosome; blocked 64×64 bit transpose into haplotype-major rows; binary and text writers; index files |
| `LDCalculator.hpp/.cpp` | r² by AND + popcount over 64-bit words (portable popcount wrapper, written so the compiler can auto-vectorise); sliding window over SNP-major data; OpenMP with per-thread bins; decay curve; optional sparse or dense LD matrix |
| `MapEstimator.hpp/.cpp` | Realised map from breakpoint density: r̂ = count / (p·N0·T_total·bp), per window; LD-based estimate by fitting 4N·c per window to r² vs distance with the Hill–Weir expectation |
| `VcfReader.hpp/.cpp` | Phased VCF reader (gz via zlib, enabled by `ZLIB=1`, otherwise plain text); biallelic SNPs; sample subset from the 1000G panel file |
| `main.cpp` | Subcommands `simulate` and `ld`; default test case is chr22 / CEU |

### Algorithms per event (target: amortised O(log n))
- **Recombination point:** cumulative-genetic-distance draw. The branch is chosen by Fenwick search; the time is uniform on that branch.
- **Floating lineage:** step through the time intervals above t_r using the sorted internal times and the epoch boundaries.
  - The target branch crossing t_c is picked uniformly.
  - If k ≥ (2n−2)/8, use rejection sampling over all branches.
  - Otherwise, enumerate downward from the root (O(k)).
- **SMC′:** counts +1 lineage on [t_r, t_op) for the lineage's own old branch. Hitting it means no topology change.
- **MaCS:** the targets are local branches plus retained non-local branches that cross t_c.
  - On hitting a non-local branch, follow its retained parent chain upward until it meets a local branch, which gives the attachment point.
  - If the chain reaches the old root first, resume the floating process from that time.
  - Window h is a number of previous local trees. Before implementing, I'll check this against the MaCS paper and the `-h` semantics in the source.
- **Mutations:** each tree span gets Poisson(μ·p·N0·T_total·span) mutations. Each is placed on a Fenwick-chosen branch; the subtree leaves are written into a 64-bit-word bitset row. Rows go to a per-chromosome SNP-major buffer and are flushed to a temp file in blocks.

### Output (prefix `out/<name>`)
- `.hap.bin`
  - Header: magic, version, n_hap, ploidy, n_chr, total SNPs, words per row.
  - Then one row per haplotype, with **all chromosomes consecutive** (requirement 6). Each chromosome block is padded to a 64-bit word.
  - Haplotypes of individual i are rows i·p … i·p+p−1.
- `.hap.txt`: optional 0/1 text version of the same rows.
- `.chrom.tsv`: chr, length_bp, length_cM, n_snps, first_snp_idx, row_word_offset, n_words. This is what's needed to extract one chromosome.
- `.snps.tsv`: global_idx, chr, pos_bp, pos_cM, alt allele frequency.
- `.samples.tsv`: individual, haplotype index.
- `.map.true.tsv`, `.map.realised.tsv`, `.map.ldest.tsv` (only with `--estimate-map`).
- `.ld_decay.tsv` (optional), `.ld_matrix.*` (optional; sparse within `--ld-max-dist`, or dense for `--ld-region`).

### Memory and parallelism at 10k haplotypes, whole genome
- **Memory per chromosome:** tree O(n), MaCS history O(n+h), SNP buffer about 64 MB per block, LD window about (window SNPs × 1.25 KB).
- **Expected peak:** under 2–3 GB with several chromosomes in parallel.
- **OpenMP:**
  - chromosomes in parallel (dynamic schedule, largest first)
  - LD over SNP i
  - transpose over haplotype blocks
- Results are deterministic regardless of thread count, because each chromosome has its own RNG stream.
- Large outputs go to `out/`, which can be redirected with `--out`. The project folder is OneDrive-synced, so I'll recommend a non-synced path for whole-genome runs.

### CLI (sketch)
```
hapcoal simulate --model macs|smc|smcp --n-ind 5000 --ploidy 2 | --n-hap 10000
                --Ne 10000 --mu 1.25e-8 [--rho 1e-8] [--demography epochs.tsv]
                --chromosomes chroms.tsv | --chr 22:51304566 [--map hapmap_chr22.txt]
                [--history 1000] [--seed 42] [--threads 10] --out out/sim
                [--text] [--ld-decay] [--ld-matrix] [--ld-max-dist 500000] [--maf 0.05] [--estimate-map]
hapcoal ld --vcf chr22.vcf.gz --samples CEU.txt [--maf 0.05] --out out/obs
```

### Makefile
- Variables: `CXX ?= c++`, `OPENMP ?= 1`, `ZLIB ?= 1`, `NATIVE ?= 0` (adds `-march=native`), `DEBUG ?= 0`.
- OS detection with `uname`:
  - **Darwin + Apple clang:** `-Xpreprocessor -fopenmp -lomp`, using `$(brew --prefix libomp)`.
  - **Linux / GCC:** `-fopenmp`.
  - **MSYS2 / MinGW:** `-fopenmp`, `.exe` suffix.
- Flags: `-std=c++17 -O3`.
- Targets: `all`, `test`, `chr22`, `bench`, `clean`.

---

## Stage 8.2 — Implementation steps (each ends with a short report: what was done, test output, numbers)
0. **Housekeeping.**
   - Create `legacy/` and move in `macs_simulator.cpp`, the `ok_code_*.txt` snapshots, `plot_ld.py`, `output.ms`, the old CSVs, PNG and map files. Ask before moving the 1.1 GB `haplotypes.csv`, or delete it.
   - Add the `src/` skeleton, Makefile and `.gitignore`; `make` builds an empty binary on macOS.
1. **Rng, Demography, GeneticMap.** Unit tests:
   - map round-trip bp↔cM
   - HapMap parse
   - hotspot map integrates to the expected cM
   - hazard inversion vs. analytic constant-N result
2. **LocalTree with no recombination.** Check against theory with 10⁴ replicates:
   - E[TMRCA] = 2(1−1/n)
   - E[T_total] = 2·H_{n−1}
   - site frequency spectrum ∝ 1/i
   - under a bottleneck, the pairwise coalescence time distribution matches the analytic result
3. **SMC and SMC′ transitions.** Check:
   - breakpoint count = ρL·E[T_total]/2
   - tree validity checker (debug build)
   - n = 2 SMC′ two-locus covariance of T against the analytic formula
   - r² decay compared with msprime `smc` / `smc_prime` (optional Python venv in scratchpad; network install, I'll ask first)
4. **MaCS history.**
   - h = 0 reproduces SMC′ exactly (same seed, identical output).
   - Increasing h moves LD decay toward the full Hudson coalescent (msprime `hudson`).
   - Measure the cost of h.
5. **Mutations + HaplotypeStore.** Streaming, transpose, binary/text writers, index files; round-trip test against a naive dense matrix on a small case.
6. **LDCalculator.** Check against a naive double-precision loop; OpenMP speedup; sparse/dense matrix output.
7. **MapEstimator.** On a simulated hotspot map, the realised map and LD estimate should track the true map (report correlation per window).
8. **VcfReader + `ld` mode.**
   - Download into `data/`: 1000G phase 3 chr22 VCF (about 200 MB), the sample panel, and the HapMap II GRCh37 chr22 map.
   - Compute the observed CEU LD decay: biallelic SNPs, MAF ≥ 0.05.
9. **chr22 test case (`main.cpp`, `make chr22`).**
   - Simulate 99 CEU-sized individuals × 2 haplotypes with the HapMap map, μ = 1.25e-8, and a published CEU piecewise N(t) (simplified Tennessen/Gravel-style epochs), for all three models.
   - Compare the simulated and observed decay curves (RMSE per distance bin).
   - Updated `scripts/plot_ld.py` (needs matplotlib; I'll ask before installing).
10. **Performance.**
    - Benchmark 10k haplotypes on chr22, then all 22 autosomes (GRCh37 lengths, HapMap maps).
    - Measure wall time, peak RSS (`/usr/bin/time -l`) and thread scaling.
    - Profile and optimise hotspots.
    - Report the final numbers.

## Verification (end-to-end)
- `make clean && make && make test` passes, with both Apple clang + libomp and g++-15.
- `make chr22` produces every output file. The simulated LD decay for CEU-like settings lies close to the observed 1000G CEU curve (RMSE reported per model). MaCS and SMC′ should match better than SMC.
- Extracting chr22 from `.hap.bin` using `.chrom.tsv` gives the same bits as the text output.
- The whole-genome 10k-haplotype run finishes within the memory budget (under 16 GB, target under 4 GB), with reported runtime.
