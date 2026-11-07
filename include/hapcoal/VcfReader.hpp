/**
 * @file VcfReader.hpp
 * @brief Phased VCF input and observed LD of real data.
 *
 * Reads biallelic SNPs (single-base REF and ALT, not flagged MULTI_ALLELIC)
 * from a phased VCF (`.vcf`, or `.vcf.gz` when built with zlib) for a subset
 * of samples and packs the haplotypes (2 per diploid sample, in sample order)
 * into SNP-major BitRows — the representation used for simulated data — so that
 * observed and simulated LD are computed by the same code.
 */
#pragma once
#include <string>
#include <vector>

#include "hapcoal/Haplotypes.hpp"
#include "hapcoal/LDCalculator.hpp"

namespace hapcoal {

/// Haplotypes read from a VCF.
struct VcfData {
    std::string chrom;                 ///< chromosome name (normalised, "22")
    BitRows rows;                      ///< SNP-major haplotypes of the kept SNPs
    std::vector<std::string> samples;  ///< sample IDs in haplotype order (2 haplotypes each)
    size_t linesRead = 0;              ///< data lines read
    size_t snpsKept = 0;               ///< SNPs kept (biallelic, in region, MAF filter)
    size_t unphased = 0;               ///< genotypes with '/' separator (treated as phased)
    size_t missing = 0;                ///< genotypes with '.' (treated as 0)
};

/// Sample IDs of population @p pop (panel column 2) or super-population (column 3)
/// from a 1000 Genomes panel file ("sample pop super_pop gender").
std::vector<std::string> readPanelSamples(const std::string& panelPath, const std::string& pop);

/**
 * @brief Read a phased VCF.
 * @param path         VCF file (.vcf or .vcf.gz)
 * @param samples      sample IDs to keep (empty = all)
 * @param minMaf       minimum minor-allele frequency of kept SNPs
 * @param regionStart, regionEnd  1-based inclusive position range
 * @throws std::runtime_error on missing file / header or if no requested sample is present
 */
VcfData readPhasedVcf(const std::string& path, const std::vector<std::string>& samples, double minMaf,
                      double regionStart = 0, double regionEnd = 1e18);

/// Result of observedLD().
struct ObservedLD {
    std::string chrom;     ///< chromosome
    LDDecay decay;         ///< LD decay of the observed haplotypes
    size_t nHap = 0;       ///< haplotypes used
    size_t snpsUsed = 0;   ///< SNPs entering the LD calculation (after MAF filter and thinning)
    VcfData stats;         ///< reading statistics (rows are moved out; only counters are kept)
};

/**
 * @brief LD decay of real (phased) data, computed exactly like simulated LD.
 * @param vcfPath      phased VCF (one chromosome)
 * @param samples      sample IDs (empty = all; see readPanelSamples())
 * @param regionStart, regionEnd  1-based inclusive range
 * @param opts         LD options (MAF filter, distances, threads, genotype r²)
 */
ObservedLD observedLD(const std::string& vcfPath, const std::vector<std::string>& samples, double regionStart,
                      double regionEnd, const LDOptions& opts);

}  // namespace hapcoal
