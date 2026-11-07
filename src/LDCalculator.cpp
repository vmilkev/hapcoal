/**
 * @file LDCalculator.cpp
 * @brief r² computation (haplotype and genotype), decay bins, matrices and decay utilities.
 */
#include "hapcoal/LDCalculator.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Rng.hpp"

namespace hapcoal {

void LDDecay::add(const LDDecay& o) {
    if (sum.size() < o.sum.size()) {
        sum.resize(o.sum.size(), 0.0);
        n.resize(o.n.size(), 0);
    }
    if (o.binSize > 0) binSize = o.binSize;
    for (size_t i = 0; i < o.sum.size(); ++i) {
        sum[i] += o.sum[i];
        n[i] += o.n[i];
    }
}

double LDDecay::meanR2(double a, double b) const {
    double s = 0;
    uint64_t cnt = 0;
    for (size_t i = 0; i < sum.size(); ++i) {
        const double x = i * binSize;
        if (x >= a && x < b) {
            s += sum[i];
            cnt += n[i];
        }
    }
    return cnt ? s / static_cast<double>(cnt) : NAN;
}

void writeDecay(std::ostream& out, const std::string& label, const LDDecay& d) {
    for (size_t b = 0; b < d.sum.size(); ++b) {
        if (d.n[b] == 0) continue;
        out << label << '\t' << static_cast<long long>(b * d.binSize) << '\t'
            << static_cast<long long>((b + 0.5) * d.binSize) << '\t' << std::setprecision(6)
            << d.sum[b] / static_cast<double>(d.n[b]) << '\t' << d.n[b] << '\n';
    }
}

double decayRmse(const LDDecay& a, const LDDecay& b, double maxDist, double* meanDiff) {
    double s = 0, sd = 0;
    size_t n = 0;
    const size_t nb = std::min(a.sum.size(), b.sum.size());
    for (size_t i = 0; i < nb; ++i) {
        if (i * a.binSize > maxDist || a.n[i] == 0 || b.n[i] == 0) continue;
        const double d = a.sum[i] / a.n[i] - b.sum[i] / b.n[i];
        s += d * d;
        sd += d;
        ++n;
    }
    if (meanDiff) *meanDiff = n ? sd / n : 0;
    return n ? std::sqrt(s / n) : NAN;
}

/// Number of haplotypes carrying the derived allele at both SNPs (popcount of AND).
static inline uint32_t andCount(const uint64_t* __restrict x, const uint64_t* __restrict y, size_t W) {
    uint32_t c = 0;
    for (size_t w = 0; w < W; ++w) c += static_cast<uint32_t>(popcount64(x[w] & y[w]));
    return c;
}

/// Sum over diploid individuals of g_i * g_j with g = h(2k) + h(2k+1): even bits are the
/// first, odd bits the second haplotype of an individual, so the product expands into
/// four popcounts of ANDs of the even/odd bit planes.
static constexpr uint64_t kEven = 0x5555555555555555ULL;
static inline uint32_t dosageProduct(const uint64_t* __restrict x, const uint64_t* __restrict y, size_t W) {
    uint32_t c = 0;
    for (size_t w = 0; w < W; ++w) {
        const uint64_t a1 = x[w] & kEven, b1 = (x[w] >> 1) & kEven;
        const uint64_t a2 = y[w] & kEven, b2 = (y[w] >> 1) & kEven;
        c += static_cast<uint32_t>(popcount64(a1 & a2) + popcount64(a1 & b2) + popcount64(b1 & a2) +
                                   popcount64(b1 & b2));
    }
    return c;
}

double r2Pair(const uint64_t* x, const uint64_t* y, size_t words, uint32_t cx, uint32_t cy, size_t nHap) {
    const double n = static_cast<double>(nHap);
    const double pa = cx / n, pb = cy / n;
    const double denom = pa * (1 - pa) * pb * (1 - pb);
    if (denom <= 0) return 0.0;
    const double D = andCount(x, y, words) / n - pa * pb;
    return D * D / denom;
}

std::vector<size_t> selectSnps(const std::vector<uint32_t>& counts, size_t nHap, const LDOptions& opts) {
    std::vector<size_t> sel;
    const double n = static_cast<double>(nHap);
    for (size_t i = 0; i < counts.size(); ++i) {
        const double p = counts[i] / n;
        const double maf = std::min(p, 1.0 - p);
        if (maf >= opts.maf && maf > 0.0) sel.push_back(i);
    }
    if (opts.maxSnps > 0 && sel.size() > opts.maxSnps) {
        // Reservoir-free thinning: random subset of size maxSnps, kept sorted.
        Rng rng(opts.seed, 0x1D1D);
        for (size_t i = 0; i < opts.maxSnps; ++i) {
            const size_t j = i + static_cast<size_t>(rng.below(sel.size() - i));
            std::swap(sel[i], sel[j]);
        }
        sel.resize(opts.maxSnps);
        std::sort(sel.begin(), sel.end());
    }
    return sel;
}

void computeLD(const BitRows& rows, const std::string& chrom, const LDOptions& opts, double regionStart,
               double regionEnd, LDDecay& decay, LDWindows* windows, std::ostream* matrixOut) {
    const size_t S = rows.size();
    const size_t W = rows.words;
    const double nH = static_cast<double>(rows.nHap);
    const size_t nBins = static_cast<size_t>(std::ceil(opts.maxDist / opts.binSize)) + 1;
    decay.binSize = opts.binSize;
    if (decay.sum.size() < nBins) {
        decay.sum.resize(nBins, 0.0);
        decay.n.resize(nBins, 0);
    }
    if (windows) {
        windows->start = regionStart;
        windows->windowSize = opts.windowSize;
        windows->maxDist = opts.windowMaxDist;
        windows->bin = opts.windowBin;
        windows->nWin = static_cast<size_t>(std::ceil((regionEnd - regionStart) / opts.windowSize)) + 1;
        windows->nBin = static_cast<size_t>(std::ceil(opts.windowMaxDist / opts.windowBin)) + 1;
        windows->sum.assign(windows->nWin * windows->nBin, 0.0);
        windows->n.assign(windows->nWin * windows->nBin, 0);
    }
    const double maxD = std::max(opts.maxDist, (opts.matrix ? opts.matrixMaxDist : 0.0));

    // Precompute p q for each SNP (haplotype r2), or dosage sums (genotype r2).
    const bool geno = opts.genotype;
    if (geno && rows.nHap % 2 != 0) throw std::invalid_argument("genotype r2 needs an even number of haplotypes");
    const double nI = nH / 2.0;
    std::vector<double> p(S), pq(S), gs(S);
    for (size_t i = 0; i < S; ++i) {
        p[i] = rows.count[i] / nH;
        pq[i] = p[i] * (1.0 - p[i]);
        if (geno) {
            // var term: n * sum g^2 - (sum g)^2 ; sum g^2 = count + 2 * #(both copies)
            uint32_t both = 0;
            const uint64_t* r = rows.row(i);
            for (size_t w = 0; w < W; ++w) both += static_cast<uint32_t>(popcount64(r[w] & (r[w] >> 1) & kEven));
            const double sg = rows.count[i], sg2 = rows.count[i] + 2.0 * both;
            gs[i] = sg;
            pq[i] = nI * sg2 - sg * sg;  // > 0 unless monomorphic in genotypes
        }
    }

    int nThreads = 1;
#ifdef _OPENMP
    nThreads = opts.threads > 0 ? opts.threads : omp_get_max_threads();
#endif
    std::vector<LDDecay> tDecay(static_cast<size_t>(nThreads));
    std::vector<LDWindows> tWin(windows ? static_cast<size_t>(nThreads) : 0);
    for (auto& d : tDecay) { d.binSize = opts.binSize; d.sum.assign(nBins, 0.0); d.n.assign(nBins, 0); }
    for (auto& w : tWin) { w = *windows; }

#pragma omp parallel num_threads(nThreads)
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        LDDecay& D = tDecay[static_cast<size_t>(tid)];
        LDWindows* Wn = windows ? &tWin[static_cast<size_t>(tid)] : nullptr;
        std::ostringstream mbuf;
#pragma omp for schedule(dynamic, 64)
        for (long long ii = 0; ii < static_cast<long long>(S); ++ii) {
            const size_t i = static_cast<size_t>(ii);
            if (pq[i] <= 0) continue;
            const uint64_t* ri = rows.row(i);
            const size_t win = Wn ? static_cast<size_t>((rows.pos[i] - Wn->start) / Wn->windowSize) : 0;
            for (size_t j = i + 1; j < S; ++j) {
                const double d = rows.pos[j] - rows.pos[i];
                if (d > maxD) break;
                if (pq[j] <= 0) continue;
                double r2;
                if (geno) {
                    const double cov = nI * dosageProduct(ri, rows.row(j), W) - gs[i] * gs[j];
                    r2 = cov * cov / (pq[i] * pq[j]);
                } else {
                    const double pab = andCount(ri, rows.row(j), W) / nH;
                    const double Dij = pab - p[i] * p[j];
                    r2 = Dij * Dij / (pq[i] * pq[j]);
                }
                if (d <= opts.maxDist) {
                    const size_t b = static_cast<size_t>(d / opts.binSize);
                    D.sum[b] += r2;
                    ++D.n[b];
                }
                if (Wn && d <= Wn->maxDist && win < Wn->nWin) {
                    const size_t b = static_cast<size_t>(d / Wn->bin);
                    Wn->sum[win * Wn->nBin + b] += r2;
                    ++Wn->n[win * Wn->nBin + b];
                }
                if (matrixOut && d <= opts.matrixMaxDist && r2 >= opts.matrixMinR2)
                    mbuf << chrom << '\t' << rows.globalIdx[i] << '\t' << rows.globalIdx[j] << '\t'
                         << static_cast<long long>(rows.pos[i]) + 1 << '\t' << static_cast<long long>(rows.pos[j]) + 1
                         << '\t' << std::setprecision(5) << r2 << '\n';
            }
            if (matrixOut && mbuf.tellp() > (1 << 22)) {
#pragma omp critical(ldmatrix)
                *matrixOut << mbuf.str();
                mbuf.str("");
                mbuf.clear();
            }
        }
        if (matrixOut) {
#pragma omp critical(ldmatrix)
            *matrixOut << mbuf.str();
        }
    }
    for (const auto& d : tDecay) decay.add(d);
    if (windows)
        for (const auto& w : tWin)
            for (size_t k = 0; k < w.sum.size(); ++k) {
                windows->sum[k] += w.sum[k];
                windows->n[k] += w.n[k];
            }
}

void writeDenseMatrix(const BitRows& rows, double a, double b, const std::string& path) {
    std::vector<size_t> idx;
    for (size_t i = 0; i < rows.size(); ++i)
        if (rows.pos[i] >= a && rows.pos[i] <= b) idx.push_back(i);
    std::ofstream out(path);
    out << "pos";
    for (size_t i : idx) out << '\t' << static_cast<long long>(rows.pos[i]) + 1;
    out << '\n' << std::setprecision(4);
    for (size_t i : idx) {
        out << static_cast<long long>(rows.pos[i]) + 1;
        for (size_t j : idx)
            out << '\t' << (i == j ? 1.0 : r2Pair(rows.row(i), rows.row(j), rows.words, rows.count[i], rows.count[j], rows.nHap));
        out << '\n';
    }
}

}  // namespace hapcoal
