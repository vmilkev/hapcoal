/**
 * @file CoalescentSimulator.hpp
 * @brief Public interface of hapcoal: parameters, results and the simulator class.
 *
 * One class, CoalescentSimulator, implements the three sequentially Markov
 * approximations of the coalescent with recombination and the library
 * interface:
 *
 * | Model | Reference | Floating lineage may coalesce with |
 * |---|---|---|
 * | SMC   | McVean & Cardin (2005)      | branches of the remaining tree (its own old branch removed) |
 * | SMC′  | Marjoram & Wall (2006)      | the remaining tree **and** its own old branch (invisible recombination) |
 * | MaCS  | Chen, Marjoram & Wall (2009)| any edge of the ARG retained over the last `history` bp |
 *
 * Scaling: haploid population size p·N(t) (p = ploidy); time in units of
 * p·N0 generations; per-bp rates per generation are multiplied by p·N0.
 * See docs/hapcoal.md for the mathematical description.
 *
 * Minimal use:
 * @code
 *   hapcoal::SimParams p;
 *   p.model = hapcoal::Model::MaCS;
 *   p.nInd = 100;  p.ploidy = 2;
 *   p.chromosomes = {{"1", 0, 5e6, -1}};
 *   hapcoal::CoalescentSimulator sim(p);
 *   hapcoal::SimulationResult res = sim.run();
 *   int a = res.chromosomes[0].haplotypes.allele(0, 0);
 * @endcode
 */
#pragma once
#include <cstdint>
#include <functional>
#include <iosfwd>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "hapcoal/Demography.hpp"
#include "hapcoal/GeneticMap.hpp"
#include "hapcoal/Haplotypes.hpp"
#include "hapcoal/LDCalculator.hpp"
#include "hapcoal/MapEstimator.hpp"

namespace hapcoal {

class LocalTree;  // internal (src/LocalTree.hpp)
class SnpBuffer;  // internal (src/HaplotypeStore.hpp)

/// Coalescent approximation.
enum class Model {
    SMC,       ///< sequentially Markov coalescent (McVean & Cardin 2005)
    SMCPrime,  ///< SMC′ (Marjoram & Wall 2006)
    MaCS,      ///< Markovian coalescent simulator (Chen, Marjoram & Wall 2009)
};
/// "SMC", "SMC'" or "MaCS".
std::string modelName(Model m);
/// Parse "smc", "smcp"/"smc'"/"smcprime" or "macs" (case-insensitive).
Model parseModel(const std::string& s);

/// One chromosome (or region) to simulate.
struct ChromosomeSpec {
    std::string name;      ///< chromosome name ("chr" prefix is ignored when matching maps)
    double start = 0.0;    ///< region start (bp, 0-based, inclusive)
    double end = 0.0;      ///< region end (bp, exclusive); 0 = take the range of the genetic map
    double rate = -1.0;    ///< mean recombination rate per bp per generation; < 0 = SimParams::rho
};

/// All parameters of a simulation run.
struct SimParams {
    // ---------------------------------------------------- population genetics
    Model model = Model::MaCS;  ///< coalescent approximation
    int nInd = 0;               ///< number of individuals (number of haplotypes = nInd·ploidy)
    int ploidy = 2;             ///< ploidy p (haplotypes per individual; scales time as p·N)
    int nHap = 0;               ///< alternative to nInd: number of haplotypes
    double Ne = 10000.0;        ///< reference effective population size N0
    double mu = 1.25e-8;        ///< mutation rate per bp per generation
    double rho = 1e-8;          ///< recombination rate per bp per generation (if no map)
    std::vector<std::pair<double, double>> epochs;  ///< demography: (generations ago, N); empty = constant N0
    std::string demographyFile; ///< alternatively a demography file ("generations_ago N" per line)
    double history = 1e5;       ///< MaCS history window h (bp); 0 reproduces SMC′ exactly

    // ----------------------------------------------------------------- genome
    std::vector<ChromosomeSpec> chromosomes;  ///< chromosomes; empty = all chromosomes of mapFile
    std::string mapFile;        ///< optional genetic map file (any supported format, several chromosomes)
    bool hotspots = false;      ///< without a map: simulate a hotspot map (mean rate = rate/rho) instead of uniform
    GeneticMap::HotspotParams hotspotParams;  ///< hotspot landscape parameters

    // -------------------------------------------------------------------- run
    uint64_t seed = 42;         ///< master seed; chromosome i uses an independent stream (seed, i+1)
    int threads = 0;            ///< OpenMP threads (0 = all available); results do not depend on it
    std::ostream* log = nullptr;  ///< progress messages (null = silent)

    // ---------------------------------------------------------------- results
    bool keepHaplotypes = true; ///< fill ChromosomeResult::haplotypes (haplotype-major, in memory)
    std::string outPrefix;      ///< write output files PREFIX.* (empty = no files)
    bool writeText = false;     ///< also write PREFIX.hap.txt (0/1 text)
    std::string tempDir;        ///< SNP rows during simulation: empty = memory (or PREFIX_tmp if writing
                                ///< files); otherwise temporary files in this directory (bounded memory)
    bool keepTemp = false;      ///< keep temporary files

    // ---------------------------------------------------------------- analyses
    bool ldDecay = true;        ///< compute the LD decay curve
    bool ldMatrix = false;      ///< write a sparse LD matrix (files only)
    std::string ldRegion;       ///< "chr:start-end" (1-based): dense LD matrix of a region (files only)
    LDOptions ld;               ///< LD options (MAF filter, distances, thinning, genotype r²)
    bool estimateMap = false;   ///< infer a genetic map from the simulated haplotypes (LD-based)
    double realisedWindow = 1e5;///< window length of the realised map (bp)

    /// Number of haplotypes (nInd·ploidy if nInd > 0, else nHap).
    int numHaplotypes() const { return nInd > 0 ? nInd * ploidy : nHap; }
};

/// Statistics of one simulated chromosome.
struct ChromStats {
    std::string name;
    double length = 0;            ///< region length (bp)
    double lengthCM = 0;          ///< map length (cM)
    uint64_t recombinations = 0;  ///< recombination events (including invisible ones)
    uint64_t invisible = 0;       ///< events that left the local tree unchanged (SMC′, MaCS)
    uint64_t nonLocal = 0;        ///< MaCS: coalescences with non-ancestral edges
    uint64_t snps = 0;            ///< segregating sites
    size_t maxArgEdges = 0;       ///< MaCS: maximum number of retained ARG edges
    double seconds = 0;           ///< simulation time of this chromosome
};

/// Everything produced for one chromosome.
struct ChromosomeResult {
    std::string name;
    double start = 0, end = 0;             ///< simulated region [start, end) (bp, 0-based)
    std::vector<double> positions;         ///< SNP positions (bp, 0-based, continuous; file output uses floor+1)
    std::vector<uint32_t> altCounts;       ///< derived-allele count per SNP
    HaplotypeMatrix haplotypes;            ///< haplotype-major alleles (if SimParams::keepHaplotypes)
    GeneticMap map;                        ///< true map used for simulation
    RealisedMap realisedMap;               ///< realised map (events / branch length per window)
    std::vector<LDMapWindow> ldMap;        ///< LD-based map (if SimParams::estimateMap)
    LDDecay ldDecay;                       ///< LD decay (if SimParams::ldDecay)
    ChromStats stats;                      ///< simulation statistics

    /// Genetic position (cM, relative to the region start) of SNP @p j.
    double positionCM(size_t j) const { return (map.toMorgans(positions[j]) - map.toMorgans(start)) * 100.0; }
};

/// Result of CoalescentSimulator::run().
struct SimulationResult {
    size_t nHap = 0;                               ///< haplotypes
    int ploidy = 2;                                ///< ploidy (haplotypes of individual i: i·p … i·p+p−1)
    std::vector<ChromosomeResult> chromosomes;     ///< per chromosome, in input order
    LDDecay ldDecay;                               ///< LD decay over all chromosomes
    double seconds = 0;                            ///< wall time of run()

    /// Total number of SNPs.
    size_t numSnps() const;
    /**
     * @brief All chromosomes of haplotype @p h as one consecutive bit row.
     * Chromosome c starts at word offset Σ_{c'<c} ceil(nSnps_{c'}/64) (same
     * layout as the rows of `PREFIX.hap.bin`). Requires keepHaplotypes and no callback.
     */
    std::vector<uint64_t> genomeRow(size_t h) const;
    /// Word offset of chromosome @p c in genomeRow().
    size_t rowWordOffset(size_t c) const;
};

/// Called once per chromosome, in input order, from the thread that called run().
/// The callee may move data (e.g. haplotypes) out of the result.
using ChromosomeCallback = std::function<void(ChromosomeResult&)>;

/// Coalescent haplotype simulator (SMC, SMC′, MaCS) and library interface.
class CoalescentSimulator {
public:
    /**
     * @brief Validate parameters, build the demography and load the genetic map.
     * @throws std::invalid_argument for inconsistent parameters
     * @throws std::runtime_error if a file cannot be read
     */
    explicit CoalescentSimulator(SimParams p);

    /**
     * @brief Simulate all chromosomes and compute the requested outputs.
     *
     * Chromosomes are simulated in parallel (OpenMP); results are identical for
     * any number of threads. If @p onChromosome is given, each ChromosomeResult
     * is passed to it as soon as it is complete and its haplotypes are not kept
     * in the returned result (bounded memory: use together with SimParams::tempDir);
     * all other per-chromosome data are kept.
     */
    SimulationResult run(const ChromosomeCallback& onChromosome = nullptr);

    /// Parameters (after defaults were filled in by the constructor).
    const SimParams& params() const { return p_; }
    /// LD decay over all chromosomes of the last run().
    const LDDecay& decay() const { return decay_; }
    /// Statistics of the last run(), in input order.
    const std::vector<ChromStats>& stats() const { return stats_; }

    // ---------------------------------------------------------------- advanced
    /// Called for every local-tree span during simulateChromosome() (testing / custom statistics).
    struct TreeSpanCallback {
        virtual ~TreeSpanCallback() = default;
        /// The local tree @p tree is valid on [a, b) (bp).
        virtual void span(double a, double b, const LocalTree& tree) = 0;
    };

    /**
     * @brief Simulate one chromosome (the core of the three models).
     * @param spec      region; @param map  its genetic map
     * @param stream    random stream index (chromosome i of run() uses i+1)
     * @param snps      receives one SNP-major row per mutation
     * @param realised  receives events / branch length per window (may be null)
     * @param cb        called for every local tree (may be null)
     */
    ChromStats simulateChromosome(const ChromosomeSpec& spec, const GeneticMap& map, uint64_t stream,
                                  SnpBuffer& snps, RealisedMap* realised, TreeSpanCallback* cb = nullptr) const;

    /// Genetic map of a chromosome: user map (clipped to the region), hotspot map or uniform map.
    /// Fills spec.start/end from the map if spec.end == 0.
    GeneticMap buildMap(ChromosomeSpec& spec, uint64_t stream) const;

private:
    SimParams p_;
    Demography demo_;
    std::map<std::string, GeneticMap> userMaps_;
    std::vector<ChromStats> stats_;
    LDDecay decay_;
};

}  // namespace hapcoal
