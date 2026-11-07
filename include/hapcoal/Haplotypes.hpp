/**
 * @file Haplotypes.hpp
 * @brief Packed-bit haplotype containers and the reader for `.hap.bin` files.
 *
 * Alleles are biallelic (0 = ancestral, 1 = derived) and stored one bit per
 * allele in 64-bit words, least significant bit first. Two layouts are used:
 *
 *  - **SNP-major** (BitRows): one row per SNP, bit *h* = allele of haplotype *h*.
 *    This layout is what the simulator produces and what LD calculations need
 *    (r² of two SNPs = popcount of the AND of two rows).
 *  - **Haplotype-major** (HaplotypeMatrix): one row per haplotype, bit *j* =
 *    allele at SNP *j*. This is the natural layout for downstream users
 *    (genomes of individuals) and the layout of `.hap.bin` files.
 */
#pragma once
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace hapcoal {

/// Number of 64-bit words needed to store @p bits bits.
inline size_t wordsFor(size_t bits) { return (bits + 63) / 64; }

/// Population count of a 64-bit word (compiler intrinsic where available).
inline int popcount64(uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(x);
#else
    x = x - ((x >> 1) & 0x5555555555555555ULL);
    x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
    x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
    return static_cast<int>((x * 0x0101010101010101ULL) >> 56);
#endif
}

/**
 * @brief SNP-major packed alleles of a set of SNPs (used for LD).
 *
 * Row *i* occupies `data[i*words .. i*words+words)`; bit *h* of the row is the
 * allele of haplotype *h*. Bits beyond @ref nHap in the last word are zero.
 */
struct BitRows {
    size_t nHap = 0;                ///< number of haplotypes (bits per row)
    size_t words = 0;               ///< 64-bit words per row = wordsFor(nHap)
    std::vector<uint64_t> data;     ///< size() * words packed words
    std::vector<double> pos;        ///< SNP positions (bp, 0-based, may be fractional)
    std::vector<uint32_t> count;    ///< derived-allele count per SNP
    std::vector<int64_t> globalIdx; ///< index in the genome-wide SNP list (-1 if none)

    /// Number of SNPs (rows).
    size_t size() const { return pos.size(); }
    /// Pointer to the first word of row @p i.
    const uint64_t* row(size_t i) const { return data.data() + i * words; }
};

/**
 * @brief Haplotype-major packed alleles of one chromosome.
 *
 * Row *h* (haplotype *h*) occupies `bits[h*words .. h*words+words)`; bit *j*
 * is the allele at SNP *j*. Haplotypes of individual *i* (ploidy *p*) are rows
 * *i·p … i·p+p−1*.
 */
struct HaplotypeMatrix {
    size_t nHap = 0;             ///< number of haplotypes (rows)
    size_t nSnps = 0;            ///< number of SNPs (bits per row)
    size_t words = 0;            ///< 64-bit words per row = wordsFor(nSnps)
    std::vector<uint64_t> bits;  ///< nHap * words packed words

    /// True if the matrix holds no data (e.g. haplotypes were not requested).
    bool empty() const { return bits.empty(); }
    /// Allele (0/1) of haplotype @p h at SNP @p j.
    int allele(size_t h, size_t j) const { return static_cast<int>((bits[h * words + j / 64] >> (j % 64)) & 1ULL); }
    /// Pointer to the packed row of haplotype @p h.
    const uint64_t* row(size_t h) const { return bits.data() + h * words; }
};

/**
 * @brief Read one chromosome back from a `.hap.bin` file.
 *
 * @param binPath       path of the `PREFIX.hap.bin` file
 * @param rowWordOffset `row_word_offset` of the chromosome (from `PREFIX.chrom.tsv`)
 * @param nSnps         `n_snps` of the chromosome (from `PREFIX.chrom.tsv`)
 * @return the chromosome as a HaplotypeMatrix
 * @throws std::runtime_error if the file is missing or not a hapcoal file
 */
HaplotypeMatrix readChromosomeFromBin(const std::string& binPath, size_t rowWordOffset, size_t nSnps);

}  // namespace hapcoal
