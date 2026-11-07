/**
 * @file HaplotypeStore.cpp
 * @brief (internal) SnpBuffer, bit transposition and output-file assembly.
 */
#include "HaplotypeStore.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace hapcoal {

namespace {
/// 64-bit file seek (fseek takes a long, which is 32-bit on Windows).
int seek64(std::FILE* f, uint64_t off) {
#if defined(_WIN32)
    return _fseeki64(f, static_cast<__int64>(off), SEEK_SET);
#else
    return fseeko(f, static_cast<off_t>(off), SEEK_SET);
#endif
}

constexpr uint64_t kHeaderBytes = 64;  ///< size of the .hap.bin header

/**
 * In-place 64×64 bit-matrix transpose: bit i of a[j] becomes bit j of a[i].
 * Recursive block swap (Warren, Hacker's Delight §7-3) in 6 rounds of 32
 * word operations, i.e. ~192 operations for 4096 bits.
 */
void transpose64(uint64_t a[64]) {
    uint64_t m = 0x00000000FFFFFFFFULL;
    for (int j = 32; j != 0; j >>= 1, m ^= (m << j)) {
        for (int k = 0; k < 64; k = ((k | j) + 1) & ~j) {
            const uint64_t t = ((a[k] >> j) ^ a[k | j]) & m;
            a[k] ^= t << j;
            a[k | j] ^= t;
        }
    }
}

/**
 * Transpose a chunk of @p nb SNP-major rows (W words each) into haplotype-major
 * row segments: out[h * outStride + jb] holds SNPs jb*64..jb*64+63 of haplotype h.
 * Parallel over 64-haplotype blocks.
 */
void transposeChunk(const uint64_t* in, size_t nb, size_t W, size_t nHap, uint64_t* out, size_t outStride,
                    int threads) {
    const size_t snpBlocks = (nb + 63) / 64;
#pragma omp parallel for schedule(static) num_threads(threads > 0 ? threads : 1) if (threads != 1)
    for (long long hb = 0; hb < static_cast<long long>(W); ++hb) {
        uint64_t tile[64];
        for (size_t jb = 0; jb < snpBlocks; ++jb) {
            for (size_t r = 0; r < 64; ++r) {
                const size_t j = jb * 64 + r;
                tile[r] = j < nb ? in[j * W + static_cast<size_t>(hb)] : 0;
            }
            transpose64(tile);
            for (size_t r = 0; r < 64; ++r) {
                const size_t h = static_cast<size_t>(hb) * 64 + r;
                if (h < nHap) out[h * outStride + jb] = tile[r];
            }
        }
    }
}

/// SNPs per transposition chunk: ~2 GB budget for input + output, multiple of 64.
size_t chunkSnps(size_t W, size_t nS) {
    size_t B = (2048ull << 20) / (16 * W);
    B = std::max<size_t>(64, std::min<size_t>(B, 1u << 22));
    B = std::min(B, (nS + 63) / 64 * 64);  // never larger than the chromosome
    B -= B % 64;
    return B;
}
}  // namespace

// ============================================================== SnpBuffer ====
SnpBuffer::SnpBuffer(size_t nHap, const std::string& path) : path_(path), words_(wordsFor(nHap)) {
    if (!path_.empty()) {
        f_ = std::fopen(path_.c_str(), "wb");
        if (!f_) throw std::runtime_error("Cannot create temporary SNP file: " + path_);
        maxRows_ = std::max<size_t>(64, (16u << 20) / (words_ * 8));  // 16 MB write buffer
        buf_.assign(maxRows_ * words_, 0);
    }
}

SnpBuffer::~SnpBuffer() {
    try {
        close();
        release();
    } catch (...) {
    }
}

uint64_t* SnpBuffer::beginRow() {
    if (path_.empty()) {  // memory mode: append a zero row
        buf_.resize(buf_.size() + words_, 0);
        return buf_.data() + buf_.size() - words_;
    }
    if (rowsInBuf_ == maxRows_) flush();
    uint64_t* r = buf_.data() + rowsInBuf_ * words_;
    std::memset(r, 0, words_ * sizeof(uint64_t));
    return r;
}

void SnpBuffer::commitRow(double pos, uint32_t altCount) {
    if (!path_.empty()) ++rowsInBuf_;
    ++nRows_;
    positions_.push_back(pos);
    counts_.push_back(altCount);
}

void SnpBuffer::flush() {
    if (rowsInBuf_ == 0 || !f_) return;
    if (std::fwrite(buf_.data(), sizeof(uint64_t), rowsInBuf_ * words_, f_) != rowsInBuf_ * words_)
        throw std::runtime_error("Write error on " + path_);
    rowsInBuf_ = 0;
}

void SnpBuffer::close() {
    if (!f_) return;
    flush();
    std::fclose(f_);
    f_ = nullptr;
    buf_.clear();
    buf_.shrink_to_fit();
}

void SnpBuffer::readRows(size_t first, size_t n, uint64_t* dst) const {
    if (n == 0) return;
    if (path_.empty()) {
        std::memcpy(dst, buf_.data() + first * words_, n * words_ * sizeof(uint64_t));
        return;
    }
    if (!rf_) {
        rf_ = std::fopen(path_.c_str(), "rb");
        if (!rf_) throw std::runtime_error("Cannot open " + path_);
    }
    if (seek64(rf_, static_cast<uint64_t>(first) * words_ * 8) != 0 ||
        std::fread(dst, sizeof(uint64_t) * words_, n, rf_) != n)
        throw std::runtime_error("Temporary SNP file truncated: " + path_);
}

void SnpBuffer::release() {
    if (rf_) {
        std::fclose(rf_);
        rf_ = nullptr;
    }
    if (!path_.empty() && !f_ && !keep_) std::remove(path_.c_str());
    buf_.clear();
    buf_.shrink_to_fit();
}

// =============================================================== loadRows ====
BitRows loadRows(const SnpBuffer& snps, size_t nHap, const std::vector<size_t>& selected,
                 const std::vector<double>& pos, const std::vector<uint32_t>& count, size_t firstGlobal) {
    BitRows R;
    R.nHap = nHap;
    R.words = wordsFor(nHap);
    R.data.resize(selected.size() * R.words);
    R.pos.reserve(selected.size());
    R.count.reserve(selected.size());
    R.globalIdx.reserve(selected.size());
    // Read sequentially in 16 MB chunks and keep the selected rows.
    const size_t total = pos.size();
    const size_t chunkRows = std::max<size_t>(64, (16u << 20) / (R.words * 8));
    std::vector<uint64_t> chunk(chunkRows * R.words);
    size_t nextRow = 0, si = 0;
    while (si < selected.size() && nextRow < total) {
        const size_t got = std::min(chunkRows, total - nextRow);
        snps.readRows(nextRow, got, chunk.data());
        while (si < selected.size() && selected[si] < nextRow + got) {
            const size_t r = selected[si] - nextRow;
            std::memcpy(R.data.data() + R.pos.size() * R.words, chunk.data() + r * R.words, R.words * 8);
            R.pos.push_back(pos[selected[si]]);
            R.count.push_back(count[selected[si]]);
            R.globalIdx.push_back(static_cast<int64_t>(firstGlobal + selected[si]));
            ++si;
        }
        nextRow += got;
    }
    if (si != selected.size()) throw std::runtime_error("SNP buffer has fewer rows than selected");
    return R;
}

// ====================================================== toHaplotypeMatrix ====
HaplotypeMatrix toHaplotypeMatrix(const SnpBuffer& snps, size_t nHap, int threads) {
    HaplotypeMatrix M;
    const size_t nS = snps.numSnps(), W = wordsFor(nHap);
    M.nHap = nHap;
    M.nSnps = nS;
    M.words = wordsFor(nS);
    M.bits.assign(nHap * M.words, 0);
    if (nS == 0) return M;
    const size_t B = chunkSnps(W, nS);
    std::vector<uint64_t> inBuf(B * W);
    const size_t outStride = B / 64;
    std::vector<uint64_t> outBuf(64 * W * outStride);
    for (size_t j0 = 0; j0 < nS; j0 += B) {
        const size_t nb = std::min(B, nS - j0);
        snps.readRows(j0, nb, inBuf.data());
        std::fill(outBuf.begin(), outBuf.end(), 0);
        transposeChunk(inBuf.data(), nb, W, nHap, outBuf.data(), outStride, threads);
        const size_t blocks = (nb + 63) / 64;
        for (size_t h = 0; h < nHap; ++h)
            std::memcpy(M.bits.data() + h * M.words + j0 / 64, outBuf.data() + h * outStride, blocks * 8);
    }
    return M;
}

// ===================================================== assembleHaplotypes ====
void assembleHaplotypes(const std::string& prefix, std::vector<ChromData>& chroms, size_t nHap, int ploidy,
                        bool writeText, int threads) {
    const size_t W = wordsFor(nHap);
    // ---- global SNP indices and word offsets (each chromosome padded to 64 bits)
    size_t totalWords = 0, totalSnps = 0;
    for (auto& c : chroms) {
        c.firstSnp = totalSnps;
        c.rowWordOffset = totalWords;
        totalSnps += c.pos.size();
        totalWords += wordsFor(c.pos.size());
    }

    // ---- index files
    {
        std::ofstream ct(prefix + ".chrom.tsv");
        ct << "chrom\tstart_bp\tend_bp\tlength_cM\tn_snps\tfirst_snp_index\trow_word_offset\tn_words\n";
        for (const auto& c : chroms)
            ct << c.name << '\t' << std::fixed << std::setprecision(0) << c.start << '\t' << c.end << '\t'
               << std::setprecision(4) << c.map.totalMorgans() * 100.0 << '\t' << c.pos.size() << '\t'
               << c.firstSnp << '\t' << c.rowWordOffset << '\t' << wordsFor(c.pos.size()) << '\n';
    }
    {
        std::ofstream st(prefix + ".snps.tsv");
        st << "index\tchrom\tpos_bp\tpos_cM\talt_freq\n";
        for (const auto& c : chroms)
            for (size_t j = 0; j < c.pos.size(); ++j) {
                const double cM = (c.map.toMorgans(c.pos[j]) - c.map.toMorgans(c.start)) * 100.0;
                st << (c.firstSnp + j) << '\t' << c.name << '\t'
                   << static_cast<long long>(std::floor(c.pos[j])) + 1 << '\t' << std::fixed
                   << std::setprecision(6) << cM << '\t' << std::setprecision(5)
                   << static_cast<double>(c.count[j]) / static_cast<double>(nHap) << '\n';
            }
    }
    {
        std::ofstream sm(prefix + ".samples.tsv");
        sm << "hap_index\tindividual\thaplotype\n";
        for (size_t h = 0; h < nHap; ++h)
            sm << h << "\tind" << (h / static_cast<size_t>(ploidy)) << '\t' << (h % static_cast<size_t>(ploidy)) << '\n';
    }

    // ---- binary haplotype-major file: header
    const std::string binPath = prefix + ".hap.bin";
    std::FILE* out = std::fopen(binPath.c_str(), "wb+");
    if (!out) throw std::runtime_error("Cannot create " + binPath);
    {
        uint64_t hdr[8] = {0};
        std::memcpy(hdr, "HAPCOAL1", 8);
        hdr[1] = 1;  // version
        hdr[2] = nHap;
        hdr[3] = static_cast<uint64_t>(ploidy);
        hdr[4] = chroms.size();
        hdr[5] = totalSnps;
        hdr[6] = totalWords;
        hdr[7] = kHeaderBytes;
        std::fwrite(hdr, 8, 8, out);
    }
    // Pre-allocate the file by writing zeros sequentially: the scattered row
    // segment writes below then overwrite allocated blocks instead of filling a
    // sparse file hole by hole (which is very slow on some file systems).
    if (totalWords > 0 && nHap > 0) {
        std::vector<char> zeros(64u << 20, 0);
        uint64_t left = static_cast<uint64_t>(nHap) * totalWords * 8;
        while (left > 0) {
            const size_t n = static_cast<size_t>(std::min<uint64_t>(left, zeros.size()));
            if (std::fwrite(zeros.data(), 1, n, out) != n) throw std::runtime_error("Write error on " + binPath);
            left -= n;
        }
        std::fflush(out);
    }

    // ---- per chromosome: transpose chunks of B SNPs and write row segments
    for (const auto& c : chroms) {
        const size_t nS = c.pos.size();
        if (nS == 0) continue;
        const size_t B = chunkSnps(W, nS);
        std::vector<uint64_t> inBuf(B * W);
        const size_t outWordsPerHapMax = B / 64;
        std::vector<uint64_t> outBuf(static_cast<size_t>(64) * W * outWordsPerHapMax);
        for (size_t j0 = 0; j0 < nS; j0 += B) {
            const size_t nb = std::min(B, nS - j0);
            c.snps->readRows(j0, nb, inBuf.data());
            const size_t snpBlocks = (nb + 63) / 64;
            std::fill(outBuf.begin(), outBuf.end(), 0);
            transposeChunk(inBuf.data(), nb, W, nHap, outBuf.data(), outWordsPerHapMax, threads);
            // Scattered row-segment writes: one unbuffered handle per thread, each
            // writing a contiguous range of haplotypes (parallel I/O, no stdio copies).
            const size_t wordStart = c.rowWordOffset + j0 / 64;
            const int nt = threads > 0 ? threads : 1;
            std::string ioError;
#pragma omp parallel num_threads(nt) if (nt > 1)
            {
                int tid = 0, tn = 1;
#ifdef _OPENMP
                tid = omp_get_thread_num();
                tn = omp_get_num_threads();
#endif
                const size_t h0 = nHap * static_cast<size_t>(tid) / static_cast<size_t>(tn);
                const size_t h1 = nHap * static_cast<size_t>(tid + 1) / static_cast<size_t>(tn);
                std::FILE* fo = std::fopen(binPath.c_str(), "r+b");
                if (!fo) {
#pragma omp critical(ioerr)
                    ioError = "Cannot reopen " + binPath;
                } else {
                    std::setvbuf(fo, nullptr, _IONBF, 0);
                    for (size_t h = h0; h < h1; ++h) {
                        const uint64_t off = kHeaderBytes + (static_cast<uint64_t>(h) * totalWords + wordStart) * 8;
                        if (seek64(fo, off) != 0 ||
                            std::fwrite(outBuf.data() + h * outWordsPerHapMax, 8, snpBlocks, fo) != snpBlocks) {
#pragma omp critical(ioerr)
                            ioError = "Write error on " + binPath;
                            break;
                        }
                    }
                    std::fclose(fo);
                }
            }
            if (!ioError.empty()) throw std::runtime_error(ioError);
        }
    }
    std::fclose(out);

    // ---- optional text output: one 0/1 line per haplotype, chromosomes consecutive
    if (writeText) {
        std::FILE* bin = std::fopen(binPath.c_str(), "rb");
        std::FILE* txt = std::fopen((prefix + ".hap.txt").c_str(), "wb");
        if (!bin || !txt) throw std::runtime_error("Cannot write text haplotypes");
        std::vector<uint64_t> row(totalWords);
        std::string line;
        line.reserve(totalSnps + 1);
        seek64(bin, kHeaderBytes);
        for (size_t h = 0; h < nHap; ++h) {
            if (std::fread(row.data(), 8, totalWords, bin) != totalWords) throw std::runtime_error("Short read on " + binPath);
            line.clear();
            for (const auto& c : chroms)
                for (size_t j = 0; j < c.pos.size(); ++j)
                    line.push_back(((row[c.rowWordOffset + j / 64] >> (j % 64)) & 1ULL) ? '1' : '0');
            line.push_back('\n');
            std::fwrite(line.data(), 1, line.size(), txt);
        }
        std::fclose(bin);
        std::fclose(txt);
    }
}

// ================================================== readChromosomeFromBin ====
HaplotypeMatrix readChromosomeFromBin(const std::string& binPath, size_t rowWordOffset, size_t nSnps) {
    std::FILE* f = std::fopen(binPath.c_str(), "rb");
    if (!f) throw std::runtime_error("Cannot open " + binPath);
    uint64_t hdr[8];
    if (std::fread(hdr, 8, 8, f) != 8 || std::memcmp(hdr, "HAPCOAL1", 8) != 0) {
        std::fclose(f);
        throw std::runtime_error("Not a hapcoal haplotype file: " + binPath);
    }
    HaplotypeMatrix M;
    M.nHap = hdr[2];
    M.nSnps = nSnps;
    M.words = wordsFor(nSnps);
    M.bits.assign(M.nHap * M.words, 0);
    const size_t totalWords = hdr[6];
    for (size_t h = 0; h < M.nHap; ++h) {
        if (seek64(f, hdr[7] + (static_cast<uint64_t>(h) * totalWords + rowWordOffset) * 8) != 0 ||
            std::fread(M.bits.data() + h * M.words, 8, M.words, f) != M.words) {
            std::fclose(f);
            throw std::runtime_error("Short read on " + binPath);
        }
    }
    std::fclose(f);
    return M;
}

}  // namespace hapcoal
