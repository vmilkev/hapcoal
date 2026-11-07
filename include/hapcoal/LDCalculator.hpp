/**
 * @file LDCalculator.hpp
 * @brief Linkage disequilibrium (r²) from packed haplotypes: decay curves and matrices.
 *
 * **Haplotype r²** (Hill & Robertson 1968). For SNPs i, j with derived-allele
 * frequencies p_i, p_j and joint frequency p_ij among n haplotypes,
 *   D = p_ij − p_i p_j,   r² = D² / (p_i(1−p_i) p_j(1−p_j)).
 * p_ij = popcount(row_i AND row_j)/n is computed over 64-bit words with a
 * branch-free loop that compilers vectorise (NEON `cnt`, AVX2/AVX-512).
 *
 * **Genotype r²** (LDOptions::genotype): squared Pearson correlation of the
 * diploid dosages g = h(2k) + h(2k+1) of consecutive haplotype pairs; it does
 * not depend on phase and is therefore robust to phasing (switch) errors in
 * real data.
 *
 * Pairs within LDOptions::maxDist are averaged in distance bins (LD decay
 * curve); optionally per-window bins are collected for LD-based map estimation
 * (MapEstimator.hpp) and a sparse matrix is written. OpenMP parallelises over
 * the first SNP of each pair with per-thread accumulators.
 */
#pragma once
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

#include "hapcoal/Haplotypes.hpp"

namespace hapcoal {

/// Options of the LD calculations.
struct LDOptions {
    double maf = 0.05;           ///< minor-allele-frequency filter (SNPs with MAF < maf are ignored)
    double maxDist = 500000.0;   ///< maximum pair distance for the decay curve (bp)
    double binSize = 1000.0;     ///< width of the decay distance bins (bp)
    size_t maxSnps = 200000;     ///< random thinning to at most this many SNPs per chromosome (0 = all)
    uint64_t seed = 1;           ///< seed of the thinning
    int threads = 0;             ///< OpenMP threads (0 = all available)
    bool genotype = false;       ///< genotype (dosage) r² instead of haplotype r² (diploid pairing)
    // sparse matrix output
    bool matrix = false;             ///< write pairs to the matrix stream
    double matrixMaxDist = 100000.0; ///< maximum pair distance in the matrix (bp)
    double matrixMinR2 = 0.0;        ///< only pairs with r² >= this value
    // per-window bins (LD-based genetic map estimation)
    bool windows = false;            ///< collect per-window distance bins
    double windowSize = 1e6;         ///< window length (bp); a pair belongs to the window of its first SNP
    double windowMaxDist = 100000.0; ///< maximum pair distance used in windows (bp)
    double windowBin = 2000.0;       ///< distance bin width inside windows (bp)
};

/// LD decay curve: sum of r² and number of pairs per distance bin [b·binSize, (b+1)·binSize).
struct LDDecay {
    double binSize = 0;          ///< bin width (bp)
    std::vector<double> sum;     ///< Σ r² per bin
    std::vector<uint64_t> n;     ///< number of pairs per bin
    /// Add the bins of @p o (bin sizes must agree).
    void add(const LDDecay& o);
    /// Mean r² over all pairs with distance in [a, b) bp (NaN if no pairs).
    double meanR2(double a, double b) const;
};

/// Per-window distance bins (input of estimateLDMap()).
struct LDWindows {
    double start = 0, windowSize = 0, maxDist = 0, bin = 0;
    size_t nWin = 0, nBin = 0;
    std::vector<double> sum;   ///< nWin × nBin sums of r²
    std::vector<uint64_t> n;   ///< nWin × nBin pair counts
};

/**
 * @brief SNPs passing the MAF filter, randomly thinned to opts.maxSnps.
 * @param counts derived-allele count per SNP
 * @param nHap   number of haplotypes
 * @param opts   MAF filter, maximum number of SNPs and thinning seed
 * @return sorted indices of the selected SNPs
 */
std::vector<size_t> selectSnps(const std::vector<uint32_t>& counts, size_t nHap, const LDOptions& opts);

/**
 * @brief LD of all SNP pairs within the maximum distance on one chromosome.
 * @param rows        selected SNP rows (sorted by position)
 * @param chrom       chromosome name (matrix output only)
 * @param opts        LD options
 * @param regionStart, regionEnd  region (bp) used to lay out the per-window bins
 * @param decay       [in,out] decay bins (accumulated)
 * @param windows     per-window bins (may be null; requires opts.windows)
 * @param matrixOut   receives lines "chrom snp_i snp_j pos_i pos_j r2" (may be null)
 */
void computeLD(const BitRows& rows, const std::string& chrom, const LDOptions& opts, double regionStart,
               double regionEnd, LDDecay& decay, LDWindows* windows, std::ostream* matrixOut);

/// Dense tab-separated r² matrix of all rows with position in [a, b] (bp, 0-based).
void writeDenseMatrix(const BitRows& rows, double a, double b, const std::string& path);

/// Haplotype r² of two rows (@p cx, @p cy are their derived-allele counts).
double r2Pair(const uint64_t* x, const uint64_t* y, size_t words, uint32_t cx, uint32_t cy, size_t nHap);

/// Write decay bins as lines "label  bin_start_bp  bin_mid_bp  mean_r2  n_pairs" (bins with pairs only).
void writeDecay(std::ostream& out, const std::string& label, const LDDecay& d);

/**
 * @brief Root-mean-square difference of the mean r² of two decay curves.
 * Bins with pairs in both curves and start ≤ @p maxDist are compared.
 * @param a, b      the two decay curves (same bin size)
 * @param maxDist   largest bin start included (bp)
 * @param meanDiff  optional: mean of (a − b) over the compared bins
 */
double decayRmse(const LDDecay& a, const LDDecay& b, double maxDist, double* meanDiff = nullptr);

}  // namespace hapcoal
