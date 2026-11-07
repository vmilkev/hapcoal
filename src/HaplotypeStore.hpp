/**
 * @file HaplotypeStore.hpp
 * @brief (internal) Storage of simulated SNP rows and assembly of the output files.
 *
 * During simulation every new SNP is appended as a SNP-major row of nHap bits
 * to a SnpBuffer, which keeps the rows in memory or streams them to a temporary
 * file (bounded memory for genome-scale runs). After simulation the rows are
 *  - transposed in 64×64-bit tiles into the haplotype-major `.hap.bin` file,
 *    with all chromosomes of a haplotype stored consecutively, plus index files;
 *  - and/or transposed into an in-memory HaplotypeMatrix per chromosome;
 *  - and selectively re-read (common SNPs only) for LD calculations.
 *
 * File layout of `PREFIX.hap.bin` (all integers little-endian uint64):
 *   header (64 bytes): "HAPCOAL1", version=1, nHap, ploidy, nChrom, nSnpsTotal,
 *                      wordsPerRow, headerBytes=64
 *   then nHap rows of wordsPerRow words; chromosome c occupies words
 *   [rowWordOffset_c, rowWordOffset_c + ceil(nSnps_c/64)) of every row.
 */
#pragma once
#include <cstdint>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "hapcoal/GeneticMap.hpp"
#include "hapcoal/Haplotypes.hpp"

namespace hapcoal {

/**
 * @brief Append-only SNP-major row store, in memory or in a temporary file.
 *
 * Usage: `uint64_t* r = beginRow(); ...set bits...; commitRow(pos, count);`
 * repeated, then close(); afterwards rows can be read with readRows().
 */
class SnpBuffer {
public:
    /// @param nHap  haplotypes per row
    /// @param path  temporary file path; empty = keep rows in memory
    SnpBuffer(size_t nHap, const std::string& path = "");
    ~SnpBuffer();
    SnpBuffer(const SnpBuffer&) = delete;
    SnpBuffer& operator=(const SnpBuffer&) = delete;

    /// Zeroed row for the next SNP; valid until commitRow().
    uint64_t* beginRow();
    /// Commit the row returned by beginRow() at position @p pos with derived count @p altCount.
    void commitRow(double pos, uint32_t altCount);
    /// Finish writing (flushes the file buffer). Further reads are allowed.
    void close();
    /// Copy rows [first, first+n) into @p dst (n * words() words).
    void readRows(size_t first, size_t n, uint64_t* dst) const;
    /// Free memory / delete the temporary file (unless keepFile() was called).
    void release();
    /// Keep the temporary file on release/destruction.
    void keepFile() { keep_ = true; }

    size_t numSnps() const { return nRows_; }  ///< rows committed (independent of moved-out positions)
    size_t words() const { return words_; }
    bool inMemory() const { return path_.empty(); }
    std::vector<double>& positions() { return positions_; }
    std::vector<uint32_t>& counts() { return counts_; }

private:
    std::string path_;              // empty: memory mode
    bool keep_ = false;             // keep the temporary file
    size_t words_;
    std::FILE* f_ = nullptr;        // write handle (file mode)
    mutable std::FILE* rf_ = nullptr;  // read handle (file mode)
    std::vector<uint64_t> buf_;     // file mode: write buffer; memory mode: all rows
    size_t rowsInBuf_ = 0, maxRows_ = 0;
    size_t nRows_ = 0;              // rows committed
    std::vector<double> positions_;
    std::vector<uint32_t> counts_;
    void flush();
};

/// Everything known about one simulated chromosome when assembling outputs.
struct ChromData {
    std::string name;
    double start = 0, end = 0;           ///< simulated region [start, end) in bp
    std::unique_ptr<SnpBuffer> snps;     ///< SNP-major rows
    std::vector<double> pos;             ///< SNP positions (bp)
    std::vector<uint32_t> count;         ///< derived-allele counts
    GeneticMap map;                      ///< map used for simulation ("true" map)
    size_t firstSnp = 0;                 ///< global index of the first SNP (set by assembleHaplotypes)
    size_t rowWordOffset = 0;            ///< word offset in .hap.bin rows (set by assembleHaplotypes)
};

/// Read the SNP rows with indices @p selected (sorted ascending) into a BitRows.
BitRows loadRows(const SnpBuffer& snps, size_t nHap, const std::vector<size_t>& selected,
                 const std::vector<double>& pos, const std::vector<uint32_t>& count, size_t firstGlobal);

/// Transpose all rows of @p snps into a haplotype-major matrix.
HaplotypeMatrix toHaplotypeMatrix(const SnpBuffer& snps, size_t nHap, int threads);

/// Write PREFIX.hap.bin (+ .hap.txt if @p writeText), .chrom.tsv, .snps.tsv and .samples.tsv.
void assembleHaplotypes(const std::string& prefix, std::vector<ChromData>& chroms, size_t nHap, int ploidy,
                        bool writeText, int threads);

}  // namespace hapcoal
