/**
 * @file VcfReader.cpp
 * @brief Minimal phased-VCF reader (plain or gzip) and observed LD.
 */
#include "hapcoal/VcfReader.hpp"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#ifdef HAPCOAL_HAVE_ZLIB
#include <zlib.h>
#endif

#include "hapcoal/GeneticMap.hpp"

namespace hapcoal {

namespace {
/// Line reader over plain or gzip-compressed files (4 MB read buffer, memchr line split).
class LineReader {
public:
    explicit LineReader(const std::string& path) {
        gz_ = path.size() > 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
        if (gz_) {
#ifdef HAPCOAL_HAVE_ZLIB
            g_ = gzopen(path.c_str(), "rb");
            if (!g_) throw std::runtime_error("Cannot open " + path);
            gzbuffer(g_, 1 << 20);
#else
            throw std::runtime_error("Built without zlib (ZLIB=0); decompress first: " + path);
#endif
        } else {
            f_ = std::fopen(path.c_str(), "rb");
            if (!f_) throw std::runtime_error("Cannot open " + path);
        }
        buf_.resize(1 << 22);
    }
    ~LineReader() {
#ifdef HAPCOAL_HAVE_ZLIB
        if (g_) gzclose(g_);
#endif
        if (f_) std::fclose(f_);
    }
    bool next(std::string& line) {
        line.clear();
        for (;;) {
            if (pos_ >= len_) {
                if (!fill()) return !line.empty();
            }
            const char* s = buf_.data() + pos_;
            const char* nl = static_cast<const char*>(std::memchr(s, '\n', len_ - pos_));
            if (nl) {
                line.append(s, static_cast<size_t>(nl - s));
                pos_ += static_cast<size_t>(nl - s) + 1;
                return true;
            }
            line.append(s, len_ - pos_);
            pos_ = len_;
        }
    }

private:
    bool gz_ = false;
#ifdef HAPCOAL_HAVE_ZLIB
    gzFile g_ = nullptr;
#endif
    std::FILE* f_ = nullptr;
    std::vector<char> buf_;
    size_t pos_ = 0, len_ = 0;
    bool fill() {
        pos_ = 0;
        long n = 0;
#ifdef HAPCOAL_HAVE_ZLIB
        if (gz_) n = gzread(g_, buf_.data(), static_cast<unsigned>(buf_.size()));
        else
#endif
            n = static_cast<long>(std::fread(buf_.data(), 1, buf_.size(), f_));
        len_ = n > 0 ? static_cast<size_t>(n) : 0;
        return len_ > 0;
    }
};
}  // namespace

std::vector<std::string> readPanelSamples(const std::string& panelPath, const std::string& pop) {
    std::ifstream in(panelPath);
    if (!in) throw std::runtime_error("Cannot open panel file " + panelPath);
    std::vector<std::string> out;
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream ss(line);
        std::string id, p, sp;
        if (!(ss >> id >> p >> sp)) continue;
        if (id == "sample") continue;
        if (p == pop || sp == pop) out.push_back(id);
    }
    if (out.empty()) throw std::runtime_error("No samples of population " + pop + " in " + panelPath);
    return out;
}

VcfData readPhasedVcf(const std::string& path, const std::vector<std::string>& samples, double minMaf,
                      double regionStart, double regionEnd) {
    LineReader rd(path);
    VcfData out;
    std::string line;
    std::vector<size_t> cols;  // VCF column indices (0-based) of selected samples
    size_t nHap = 0, W = 0;
    std::vector<uint64_t> row;
    while (rd.next(line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            if (line.compare(0, 6, "#CHROM") == 0) {
                std::vector<std::string> hdr;
                std::istringstream ss(line);
                for (std::string t; std::getline(ss, t, '\t');) hdr.push_back(t);
                std::unordered_map<std::string, size_t> idx;
                for (size_t i = 9; i < hdr.size(); ++i) idx[hdr[i]] = i;
                if (samples.empty()) {
                    for (size_t i = 9; i < hdr.size(); ++i) { cols.push_back(i); out.samples.push_back(hdr[i]); }
                } else {
                    for (const auto& s : samples) {
                        auto it = idx.find(s);
                        if (it != idx.end()) { cols.push_back(it->second); out.samples.push_back(s); }
                    }
                }
                // keep column order ascending for a single left-to-right scan
                std::vector<size_t> order(cols.size());
                for (size_t i = 0; i < order.size(); ++i) order[i] = i;
                std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return cols[a] < cols[b]; });
                std::vector<size_t> c2;
                std::vector<std::string> s2;
                for (size_t i : order) { c2.push_back(cols[i]); s2.push_back(out.samples[i]); }
                cols.swap(c2);
                out.samples.swap(s2);
                if (cols.empty()) throw std::runtime_error("None of the requested samples are in " + path);
                nHap = 2 * cols.size();
                W = wordsFor(nHap);
                out.rows.nHap = nHap;
                out.rows.words = W;
                row.assign(W, 0);
            }
            continue;
        }
        ++out.linesRead;
        if (cols.empty()) throw std::runtime_error("VCF has no #CHROM header line: " + path);
        // Parse fixed columns.
        const char* s = line.c_str();
        const char* f[9];
        size_t col = 0;
        const char* p = s;
        f[0] = p;
        while (col < 8) {
            const char* t = std::strchr(p, '\t');
            if (!t) break;
            p = t + 1;
            f[++col] = p;
        }
        if (col < 8) continue;
        const auto fieldLen = [&](int k) { return static_cast<size_t>(f[k + 1] - f[k] - 1); };
        if (fieldLen(3) != 1 || fieldLen(4) != 1) continue;  // SNPs only, biallelic ALT
        const double pos = std::strtod(f[1], nullptr);
        if (pos < regionStart || pos > regionEnd) continue;
        {
            const std::string info(f[7], fieldLen(7));
            if (info.find("MULTI_ALLELIC") != std::string::npos) continue;
        }
        if (out.chrom.empty()) out.chrom = GeneticMap::normaliseChr(std::string(f[0], fieldLen(0)));
        // Scan sample fields.
        std::fill(row.begin(), row.end(), 0);
        uint32_t cnt = 0;
        size_t curCol = 8;
        const char* q = f[8];
        size_t h = 0;
        for (size_t c : cols) {
            while (curCol < c) {
                const char* t = std::strchr(q, '\t');
                if (!t) { q = nullptr; break; }
                q = t + 1;
                ++curCol;
            }
            if (!q) break;
            const char a = q[0];
            const char sep = q[1];
            const char b = (sep == '|' || sep == '/') ? q[2] : a;  // haploid calls: duplicate
            if (sep == '/') ++out.unphased;
            if (a == '.' || b == '.') ++out.missing;
            if (a == '1') { row[h / 64] |= 1ULL << (h % 64); ++cnt; }
            ++h;
            if (b == '1') { row[h / 64] |= 1ULL << (h % 64); ++cnt; }
            ++h;
        }
        if (h != nHap) continue;  // truncated line
        const double freq = static_cast<double>(cnt) / static_cast<double>(nHap);
        if (std::min(freq, 1.0 - freq) < minMaf || cnt == 0 || cnt == nHap) continue;
        out.rows.data.insert(out.rows.data.end(), row.begin(), row.end());
        out.rows.pos.push_back(pos - 1.0);  // internal coordinates are 0-based like the simulator
        out.rows.count.push_back(cnt);
        out.rows.globalIdx.push_back(static_cast<int64_t>(out.snpsKept));
        ++out.snpsKept;
    }
    return out;
}

ObservedLD observedLD(const std::string& vcfPath, const std::vector<std::string>& samples, double regionStart,
                      double regionEnd, const LDOptions& opts) {
    VcfData v = readPhasedVcf(vcfPath, samples, opts.maf, regionStart, regionEnd);
    // Same SNP selection (MAF filter + thinning) as for simulated chromosomes.
    std::vector<uint32_t> cnt = v.rows.count;
    const auto sel = selectSnps(cnt, v.rows.nHap, opts);
    BitRows rows;
    if (sel.size() == v.rows.size()) {
        rows = std::move(v.rows);
    } else {
        rows.nHap = v.rows.nHap;
        rows.words = v.rows.words;
        for (size_t i : sel) {
            rows.data.insert(rows.data.end(), v.rows.row(i), v.rows.row(i) + v.rows.words);
            rows.pos.push_back(v.rows.pos[i]);
            rows.count.push_back(v.rows.count[i]);
            rows.globalIdx.push_back(v.rows.globalIdx[i]);
        }
    }
    ObservedLD res;
    computeLD(rows, v.chrom, opts, regionStart, regionEnd, res.decay, nullptr, nullptr);
    res.chrom = v.chrom;
    res.nHap = rows.nHap;
    res.snpsUsed = rows.size();
    v.rows = BitRows{};
    res.stats = std::move(v);
    return res;
}

}  // namespace hapcoal
