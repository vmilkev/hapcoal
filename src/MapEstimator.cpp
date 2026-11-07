/**
 * @file MapEstimator.cpp
 * @brief Realised and LD-based genetic maps.
 */
#include "hapcoal/MapEstimator.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>

namespace hapcoal {

void RealisedMap::init(double s, double e, double w) {
    start = s;
    end = e;
    window = w;
    const size_t n = static_cast<size_t>(std::ceil((e - s) / w));
    events.assign(n, 0);
    integral.assign(n, 0.0);
}

void RealisedMap::addSpan(double a, double b, double T) {
    if (events.empty()) return;
    while (a < b) {
        size_t k = static_cast<size_t>((a - start) / window);
        if (k >= integral.size()) k = integral.size() - 1;
        const double wEnd = std::min(b, start + (k + 1) * window);
        const double seg = (k + 1 == integral.size()) ? b - a : wEnd - a;
        integral[k] += T * seg;
        if (k + 1 == integral.size()) break;
        a = wEnd;
    }
}

void RealisedMap::addEvent(double x) {
    if (events.empty()) return;
    size_t k = static_cast<size_t>((x - start) / window);
    if (k >= events.size()) k = events.size() - 1;
    ++events[k];
}

double hillWeirR2(double C, double n) {
    const double a = (10.0 + C) / ((2.0 + C) * (11.0 + C));
    const double b = 1.0 + ((3.0 + C) * (12.0 + 12.0 * C + C * C)) / (n * (2.0 + C) * (11.0 + C));
    return a * b;
}

void writeTrueMap(std::ostream& out, const std::string& chrom, const GeneticMap& map) {
    const auto& p = map.positions();
    const auto& c = map.cumulative();
    for (size_t i = 0; i < p.size(); ++i) {
        const double rate = (i + 1 < p.size()) ? (c[i + 1] - c[i]) / (p[i + 1] - p[i]) * 1e8 : 0.0;
        out << chrom << '\t' << static_cast<long long>(p[i]) << '\t' << std::setprecision(6) << rate << '\t'
            << std::setprecision(8) << c[i] * 100.0 << '\n';
    }
}

void writeRealisedMap(std::ostream& out, const std::string& chrom, const RealisedMap& rm, double scale) {
    double cum = 0.0;
    for (size_t k = 0; k < rm.events.size(); ++k) {
        const double a = rm.start + k * rm.window;
        const double b = std::min(rm.end, a + rm.window);
        const double rate = rm.integral[k] > 0 ? rm.events[k] / (scale * rm.integral[k]) : 0.0;  // per bp per gen
        out << chrom << '\t' << static_cast<long long>(a) << '\t' << static_cast<long long>(b) << '\t'
            << rm.events[k] << '\t' << std::setprecision(6) << rate * 1e8 << '\t' << std::setprecision(8)
            << cum * 100.0 << '\n';
        cum += rate * (b - a);
    }
}

namespace {
/// Weighted sum of squared differences between observed binned r² and the
/// Hill-Weir expectation for rho = 10^logRho in window @p win.
double sse(double logRho, const LDWindows& w, size_t win, double n) {
    const double rho = std::pow(10.0, logRho);
    double s = 0.0;
    for (size_t b = 0; b < w.nBin; ++b) {
        const uint64_t c = w.n[win * w.nBin + b];
        if (c == 0) continue;
        const double d = (b + 0.5) * w.bin;
        const double obs = w.sum[win * w.nBin + b] / static_cast<double>(c);
        const double e = hillWeirR2(rho * d, n) - obs;
        s += static_cast<double>(c) * e * e;
    }
    return s;
}
}  // namespace

std::vector<LDMapWindow> estimateLDMap(const LDWindows& w, double regionEnd, size_t nHap, double Ne, int ploidy) {
    std::vector<LDMapWindow> res;
    const double n = static_cast<double>(nHap);
    for (size_t k = 0; k < w.nWin; ++k) {
        LDMapWindow r{};
        r.start = w.start + k * w.windowSize;
        r.end = std::min(regionEnd, r.start + w.windowSize);
        if (r.start >= regionEnd) break;
        for (size_t b = 0; b < w.nBin; ++b) r.pairs += w.n[k * w.nBin + b];
        if (r.pairs >= 50) {
            // grid search on log10(rho) then golden-section refinement
            double best = -9, bestV = 1e300;
            for (double g = -9.0; g <= 0.0; g += 0.05) {
                const double v = sse(g, w, k, n);
                if (v < bestV) { bestV = v; best = g; }
            }
            double a = best - 0.05, b = best + 0.05;
            const double gr = 0.5 * (std::sqrt(5.0) - 1.0);
            for (int it = 0; it < 40; ++it) {
                const double c = b - gr * (b - a), d = a + gr * (b - a);
                if (sse(c, w, k, n) < sse(d, w, k, n)) b = d;
                else a = c;
            }
            r.rhoPerBp = std::pow(10.0, 0.5 * (a + b));
            r.cMperMb = r.rhoPerBp / (2.0 * ploidy * Ne) * 1e8;
        }
        res.push_back(r);
    }
    return res;
}

void writeLDMap(std::ostream& out, const std::string& chrom, const std::vector<LDMapWindow>& win) {
    double cum = 0.0;
    for (const auto& r : win) {
        out << chrom << '\t' << static_cast<long long>(r.start) << '\t' << static_cast<long long>(r.end) << '\t'
            << r.pairs << '\t' << std::setprecision(6) << r.rhoPerBp * 1e3 << '\t' << r.cMperMb << '\t'
            << std::setprecision(8) << cum << '\n';
        cum += r.cMperMb * (r.end - r.start) / 1e6;
    }
}

}  // namespace hapcoal
