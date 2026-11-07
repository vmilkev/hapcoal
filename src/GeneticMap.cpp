/**
 * @file GeneticMap.cpp
 * @brief Genetic maps: construction, hotspot landscapes, file formats, bp <-> Morgan conversion.
 */
#include "hapcoal/GeneticMap.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "Rng.hpp"

#ifdef HAPCOAL_HAVE_ZLIB
#include <zlib.h>
#endif

namespace hapcoal {

GeneticMap GeneticMap::uniform(double start, double end, double ratePerBp) {
    if (!(end > start)) throw std::invalid_argument("GeneticMap::uniform: end must be > start");
    GeneticMap m;
    m.pos_ = {start, end};
    m.cum_ = {0.0, ratePerBp * (end - start)};
    return m;
}

GeneticMap GeneticMap::fromPoints(const std::vector<double>& posBp, const std::vector<double>& cM) {
    if (posBp.size() != cM.size() || posBp.size() < 2)
        throw std::invalid_argument("GeneticMap::fromPoints: need >= 2 matching points");
    GeneticMap m;
    m.pos_.reserve(posBp.size());
    m.cum_.reserve(posBp.size());
    for (size_t i = 0; i < posBp.size(); ++i) {
        const double M = cM[i] / 100.0;
        if (!m.pos_.empty()) {
            if (posBp[i] <= m.pos_.back()) continue;          // drop duplicate / unsorted positions
            if (M < m.cum_.back()) throw std::runtime_error("Genetic map is not monotone in cM");
        }
        m.pos_.push_back(posBp[i]);
        m.cum_.push_back(M);
    }
    if (m.pos_.size() < 2) throw std::runtime_error("Genetic map has fewer than two usable points");
    // Shift so that the first point is 0 M (only differences matter).
    const double c0 = m.cum_.front();
    for (double& c : m.cum_) c -= c0;
    return m;
}

GeneticMap GeneticMap::withHotspots(double start, double end, double meanRatePerBp,
                                    const HotspotParams& hp, uint64_t seed, uint64_t stream) {
    Rng rng(seed, stream);
    const double L = end - start;
    if (!(L > 0)) throw std::invalid_argument("withHotspots: empty region");
    const double cs = hp.centromereStart >= 0 ? start + hp.centromereStart * L : -1.0;
    const double ce = hp.centromereEnd >= 0 ? start + hp.centromereEnd * L : -1.0;
    const auto inCentromere = [&](double x) { return cs >= 0 && x >= cs && x <= ce; };

    struct HS { double c, amp, sd; };
    std::vector<HS> hs;
    const int nh = std::max(0, static_cast<int>(std::lround(hp.densityPerMb * L / 1e6)));
    while (static_cast<int>(hs.size()) < nh) {
        const double c = rng.uniform(start, end);
        if (inCentromere(c)) continue;
        hs.push_back({c, rng.uniform(hp.ampMin, hp.ampMax), rng.uniform(hp.widthMin, hp.widthMax)});
    }
    std::sort(hs.begin(), hs.end(), [](const HS& a, const HS& b) { return a.c < b.c; });

    // Breakpoints: region ends, centromere edges, and 13 steps across +-3 sd of each hotspot.
    std::vector<double> bp = {start, end};
    if (cs >= 0) { bp.push_back(cs); bp.push_back(ce); }
    for (const auto& h : hs)
        for (int k = -6; k <= 6; ++k) {
            const double x = h.c + (k + 0.5) * 0.5 * h.sd;
            if (x > start && x < end) bp.push_back(x);
        }
    std::sort(bp.begin(), bp.end());
    bp.erase(std::unique(bp.begin(), bp.end()), bp.end());

    double maxSd = hp.widthMax;
    std::vector<double> rel(bp.size() - 1);
    for (size_t i = 0; i + 1 < bp.size(); ++i) {
        const double mid = 0.5 * (bp[i] + bp[i + 1]);
        double mult = 1.0;
        auto lo = std::lower_bound(hs.begin(), hs.end(), mid - 6 * maxSd,
                                   [](const HS& h, double v) { return h.c < v; });
        for (auto it = lo; it != hs.end() && it->c <= mid + 6 * maxSd; ++it) {
            const double d = (mid - it->c) / it->sd;
            mult += (it->amp - 1.0) * std::exp(-0.5 * d * d);
        }
        if (inCentromere(mid)) mult *= hp.centromereFactor;
        rel[i] = mult;
    }
    // Rescale so that the mean rate over the region equals meanRatePerBp.
    double integral = 0.0;
    for (size_t i = 0; i < rel.size(); ++i) integral += rel[i] * (bp[i + 1] - bp[i]);
    const double scale = meanRatePerBp * L / integral;

    GeneticMap m;
    m.pos_ = bp;
    m.cum_.resize(bp.size());
    m.cum_[0] = 0.0;
    for (size_t i = 0; i < rel.size(); ++i) m.cum_[i + 1] = m.cum_[i] + rel[i] * scale * (bp[i + 1] - bp[i]);
    return m;
}

std::string GeneticMap::normaliseChr(std::string c) {
    if (c.size() > 3 && (c.compare(0, 3, "chr") == 0 || c.compare(0, 3, "Chr") == 0 || c.compare(0, 3, "CHR") == 0))
        c = c.substr(3);
    return c;
}

namespace {
bool isNumber(const std::string& s) {
    if (s.empty()) return false;
    char* end = nullptr;
    std::strtod(s.c_str(), &end);
    return end && *end == '\0';
}

std::vector<std::string> readLines(const std::string& path) {
    std::vector<std::string> lines;
    const bool gz = path.size() > 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
    if (gz) {
#ifdef HAPCOAL_HAVE_ZLIB
        gzFile f = gzopen(path.c_str(), "rb");
        if (!f) throw std::runtime_error("Cannot open map file: " + path);
        std::string cur;
        char buf[1 << 16];
        int n;
        while ((n = gzread(f, buf, sizeof(buf))) > 0)
            for (int i = 0; i < n; ++i) {
                if (buf[i] == '\n') { lines.push_back(cur); cur.clear(); }
                else cur.push_back(buf[i]);
            }
        if (!cur.empty()) lines.push_back(cur);
        gzclose(f);
        return lines;
#else
        throw std::runtime_error("Built without zlib; decompress the map first: " + path);
#endif
    }
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open map file: " + path);
    std::string line;
    while (std::getline(in, line)) lines.push_back(line);
    return lines;
}
}  // namespace

std::map<std::string, GeneticMap> GeneticMap::loadFile(const std::string& path,
                                                       const std::string& defaultChr) {
    std::map<std::string, std::pair<std::vector<double>, std::vector<double>>> raw;
    for (std::string line : readLines(path)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::vector<std::string> tok;
        for (std::string t; ss >> t;) tok.push_back(t);
        std::string chr;
        double pos, cM;
        if (tok.size() >= 4 && isNumber(tok[1]) && isNumber(tok[2]) && isNumber(tok[3])) {
            chr = tok[0]; pos = std::stod(tok[1]); cM = std::stod(tok[3]);      // chr pos rate cM
        } else if (tok.size() >= 4 && !isNumber(tok[1]) && isNumber(tok[2]) && isNumber(tok[3])) {
            chr = tok[0]; cM = std::stod(tok[2]); pos = std::stod(tok[3]);      // PLINK: chr id cM pos
        } else if (tok.size() == 3 && isNumber(tok[0]) && isNumber(tok[1]) && isNumber(tok[2])) {
            chr = defaultChr; pos = std::stod(tok[0]); cM = std::stod(tok[2]);  // pos rate cM
        } else {
            continue;  // header or malformed
        }
        auto& r = raw[normaliseChr(chr)];
        r.first.push_back(pos);
        r.second.push_back(cM);
    }
    std::map<std::string, GeneticMap> out;
    for (auto& kv : raw) {
        const std::string& chr = kv.first;
        const auto& r = kv.second;
        // sort by position (files are normally sorted already)
        std::vector<size_t> idx(r.first.size());
        for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
        const std::vector<double>& rp = r.first;
        std::sort(idx.begin(), idx.end(), [&rp](size_t a, size_t b) { return rp[a] < rp[b]; });
        std::vector<double> p, c;
        for (size_t i : idx) { p.push_back(r.first[i]); c.push_back(r.second[i]); }
        if (p.size() >= 2) out.emplace(chr, fromPoints(p, c));
    }
    if (out.empty()) throw std::runtime_error("No usable genetic map entries in " + path);
    return out;
}

double GeneticMap::toMorgans(double bp) const {
    if (bp <= pos_.front()) return cum_.front();
    if (bp >= pos_.back()) return cum_.back();
    const size_t i = static_cast<size_t>(std::upper_bound(pos_.begin(), pos_.end(), bp) - pos_.begin()) - 1;
    const double f = (bp - pos_[i]) / (pos_[i + 1] - pos_[i]);
    return cum_[i] + f * (cum_[i + 1] - cum_[i]);
}

double GeneticMap::toBp(double M) const {
    if (M <= cum_.front()) return pos_.front();
    if (M >= cum_.back()) return pos_.back();
    // first breakpoint with cum > M; segment i = that-1 has cum_[i] <= M < cum_[i+1] (positive rate)
    const size_t i = static_cast<size_t>(std::upper_bound(cum_.begin(), cum_.end(), M) - cum_.begin()) - 1;
    const double dc = cum_[i + 1] - cum_[i];
    const double f = dc > 0 ? (M - cum_[i]) / dc : 0.0;
    return pos_[i] + f * (pos_[i + 1] - pos_[i]);
}

double GeneticMap::rateAt(double bp) const {
    if (bp < pos_.front() || bp >= pos_.back()) return 0.0;
    const size_t i = static_cast<size_t>(std::upper_bound(pos_.begin(), pos_.end(), bp) - pos_.begin()) - 1;
    return (cum_[i + 1] - cum_[i]) / (pos_[i + 1] - pos_[i]);
}

GeneticMap GeneticMap::clipped(double a, double b) const {
    a = std::max(a, start());
    b = std::min(b, end());
    if (!(b > a)) throw std::invalid_argument("GeneticMap::clipped: empty interval");
    std::vector<double> p{a}, c{toMorgans(a) * 100.0};
    for (size_t i = 0; i < pos_.size(); ++i)
        if (pos_[i] > a && pos_[i] < b) { p.push_back(pos_[i]); c.push_back(cum_[i] * 100.0); }
    p.push_back(b);
    c.push_back(toMorgans(b) * 100.0);
    return fromPoints(p, c);
}

}  // namespace hapcoal
