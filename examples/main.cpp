/**
 * @file main.cpp
 * @brief Mock client: how a class of a larger project uses the hapcoal library.
 *
 * `FounderPopulation` stands for a class of the host project (e.g. a breeding or
 * genomic-prediction simulator) that needs founder haplotypes. It only includes
 * the public header <hapcoal/hapcoal.hpp> and links lib/libhapcoal.a.
 *
 * It demonstrates the three ways of receiving results:
 *  1. in memory  — run() returns all chromosomes (haplotypes, maps, LD);
 *  2. streaming  — a callback receives each chromosome as soon as it is done,
 *                  while SNP rows are kept in temporary files (bounded memory);
 *  3. files      — the standard output files are written (outPrefix).
 *
 * Build and run:  make && ./bin/hapcoal_example
 */
#include <hapcoal/hapcoal.hpp>

#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace host {

/// Toy class of a host project: a founder population of diploid individuals.
class FounderPopulation {
public:
    /// @param nIndividuals  number of diploid founders
    /// @param seed          seed of the coalescent simulation
    FounderPopulation(int nIndividuals, uint64_t seed) {
        params_.model = hapcoal::Model::MaCS;      // MaCS with a 100 kb history window
        params_.history = 1e5;
        params_.nInd = nIndividuals;
        params_.ploidy = 2;
        params_.Ne = 10000;
        params_.mu = 1.25e-8;
        params_.seed = seed;
        // a recent bottleneck: N = 20000 now, 2000 between 200 and 1000 generations ago, 10000 before
        params_.epochs = {{0, 20000}, {200, 2000}, {1000, 10000}};
        // two chromosomes with simulated recombination hotspots
        params_.chromosomes = {{"1", 0, 4e6, 1.2e-8}, {"2", 0, 3e6, 1.0e-8}};
        params_.hotspots = true;
        params_.estimateMap = true;               // also infer a map from the haplotypes
    }

    /// 1. In-memory use: simulate and keep all haplotypes.
    void simulate() {
        hapcoal::CoalescentSimulator sim(params_);
        founders_ = sim.run();
    }

    /// Print a summary computed from the in-memory results.
    void summarise(std::ostream& os) const {
        os << "founders: " << founders_.nHap / founders_.ploidy << " individuals, " << founders_.nHap
           << " haplotypes, " << founders_.numSnps() << " SNPs, " << std::fixed << std::setprecision(2)
           << founders_.seconds << " s\n";
        for (const auto& c : founders_.chromosomes) {
            const hapcoal::HaplotypeMatrix& H = c.haplotypes;
            // heterozygosity of individual 0 = sites where its two haplotypes differ
            size_t het = 0;
            for (size_t w = 0; w < H.words; ++w) het += hapcoal::popcount64(H.row(0)[w] ^ H.row(1)[w]);
            // mean derived-allele frequency
            const double meanFreq =
                std::accumulate(c.altCounts.begin(), c.altCounts.end(), 0.0) / (c.altCounts.size() * double(H.nHap));
            os << "  chr" << c.name << ": " << c.positions.size() << " SNPs, " << std::setprecision(1)
               << c.map.totalMorgans() * 100 << " cM (realised " << realisedLengthCM(c) << " cM), "
               << c.stats.recombinations << " recombinations; individual 0 heterozygous at " << het
               << " sites; mean derived frequency " << std::setprecision(3) << meanFreq << "\n";
            os << "        LD r2: " << std::setprecision(3) << c.ldDecay.meanR2(0, 1e3) << " (<1 kb), "
               << c.ldDecay.meanR2(1e4, 2e4) << " (10-20 kb), " << c.ldDecay.meanR2(1e5, 2e5) << " (100-200 kb)\n";
            os << "        first SNP: " << std::setprecision(0) << c.positions[0] << " bp = " << std::setprecision(4)
               << c.positionCM(0) << " cM, alleles of haplotypes 0-7: ";
            for (size_t h = 0; h < 8; ++h) os << H.allele(h, 0);
            os << "\n";
        }
        // whole-genome row of haplotype 0 (chromosomes consecutive, as in .hap.bin)
        const std::vector<uint64_t> g = founders_.genomeRow(0);
        os << "  genome row of haplotype 0: " << g.size() << " words; chromosome 2 starts at word "
           << founders_.rowWordOffset(1) << "\n";
    }

    /// 2. Streaming use: each chromosome is handed over when finished; SNP rows
    ///    live in temporary files, so memory stays bounded for large genomes.
    void streamTo(std::ostream& os, const std::string& tempDir) const {
        hapcoal::SimParams p = params_;
        p.tempDir = tempDir;
        hapcoal::CoalescentSimulator sim(p);
        sim.run([&](hapcoal::ChromosomeResult& c) {
            // here a host project would copy the haplotypes into its own data structures
            os << "  received chr" << c.name << ": " << c.haplotypes.nHap << " x " << c.haplotypes.nSnps
               << " haplotype matrix\n";
        });
    }

    /// 3. File output (the same files as the `hapcoal simulate` command).
    void writeFiles(const std::string& prefix) const {
        hapcoal::SimParams p = params_;
        p.outPrefix = prefix;
        p.keepHaplotypes = false;
        hapcoal::CoalescentSimulator(p).run();
    }

private:
    /// Length of the realised map (cM) from the simulated recombination events.
    static double realisedLengthCM(const hapcoal::ChromosomeResult& c) {
        const double scale = 2.0 * 10000;  // p * N0
        double M = 0;
        const auto& rm = c.realisedMap;
        for (size_t k = 0; k < rm.events.size(); ++k)
            M += rm.rate(k, scale) * (std::min(rm.end, rm.start + (k + 1) * rm.window) - (rm.start + k * rm.window));
        return M * 100;
    }

    hapcoal::SimParams params_;
    hapcoal::SimulationResult founders_;
};

}  // namespace host

int main() {
    try {
        std::cout << "hapcoal " HAPCOAL_VERSION_STRING " library example\n\n[1] in-memory results\n";
        host::FounderPopulation pop(250, 7);
        pop.simulate();
        pop.summarise(std::cout);

        std::cout << "\n[2] streaming with a per-chromosome callback\n";
        pop.streamTo(std::cout, "out/example/tmp");

        std::cout << "\n[3] file output\n";
        pop.writeFiles("out/example/founders");
        std::cout << "  written out/example/founders.{hap.bin,chrom.tsv,snps.tsv,samples.tsv,map.*.tsv,ld_decay.tsv}\n";
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
