/**
 * @file Demography.cpp
 * @brief Piecewise-constant N(t): construction, presets and waiting-time sampling.
 */
#include "hapcoal/Demography.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "Rng.hpp"

namespace hapcoal {

Demography::Demography() : epochs_{{0.0, 1.0}} {}

Demography Demography::fromGenerations(std::vector<std::pair<double, double>> genN, double N0,
                                       int ploidy) {
    if (N0 <= 0.0 || ploidy <= 0) throw std::invalid_argument("Demography: N0 and ploidy must be > 0");
    std::sort(genN.begin(), genN.end());
    Demography d;
    d.epochs_.clear();
    if (genN.empty() || genN.front().first > 0.0) d.epochs_.push_back({0.0, 1.0});
    const double unit = static_cast<double>(ploidy) * N0;  // generations per coalescent unit
    for (const auto& [g, N] : genN) {
        if (g < 0.0 || N <= 0.0) throw std::invalid_argument("Demography: times must be >= 0 and sizes > 0");
        const double t = g / unit;
        if (!d.epochs_.empty() && d.epochs_.back().start == t)
            d.epochs_.back().lambda = N / N0;
        else
            d.epochs_.push_back({t, N / N0});
    }
    return d;
}

std::vector<std::pair<double, double>> Demography::tennessen2012CEU() {
    std::vector<std::pair<double, double>> e;
    const double N_A = 7310, N_AF = 14474, N_B = 1861, N_EU0 = 1032;
    const double T_AF = 148e3 / 25, T_B = 51e3 / 25, T_EU = 23e3 / 25, T_EG = 5115.0 / 25;  // generations
    const double r1 = 0.00307, r2 = 0.0195;                                               // growth rates
    const double N_EU1 = N_EU0 * std::exp(r1 * (T_EU - T_EG));
    // recent rapid growth: N(g) = N_EU1 * exp(r2 * (T_EG - g)), g in [0, T_EG]; 24 epochs
    const int k2 = 24;
    for (int i = 0; i < k2; ++i) {
        const double g0 = T_EG * i / k2, gm = T_EG * (i + 0.5) / k2;
        e.emplace_back(g0, N_EU1 * std::exp(r2 * (T_EG - gm)));
    }
    // slower growth: N(g) = N_EU0 * exp(r1 * (T_EU - g)), g in [T_EG, T_EU]; 16 epochs
    const int k1 = 16;
    for (int i = 0; i < k1; ++i) {
        const double g0 = T_EG + (T_EU - T_EG) * i / k1, gm = T_EG + (T_EU - T_EG) * (i + 0.5) / k1;
        e.emplace_back(g0, N_EU0 * std::exp(r1 * (T_EU - gm)));
    }
    e.emplace_back(T_EU, N_B);   // out-of-Africa bottleneck
    e.emplace_back(T_B, N_AF);   // African expansion
    e.emplace_back(T_AF, N_A);   // ancestral
    return e;
}

Demography Demography::fromFile(const std::string& path, double N0, int ploidy) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open demography file: " + path);
    std::vector<std::pair<double, double>> genN;
    std::string line;
    while (std::getline(in, line)) {
        auto h = line.find('#');
        if (h != std::string::npos) line.erase(h);
        std::istringstream ss(line);
        double g, N;
        if (ss >> g >> N) genN.emplace_back(g, N);
    }
    if (genN.empty()) throw std::runtime_error("Demography file has no epochs: " + path);
    return fromGenerations(std::move(genN), N0, ploidy);
}

size_t Demography::epochIndex(double t) const {
    // epochs_ is small (tens of entries): upper_bound on start times.
    auto it = std::upper_bound(epochs_.begin(), epochs_.end(), t,
                               [](double v, const Epoch& e) { return v < e.start; });
    return static_cast<size_t>(std::max<std::ptrdiff_t>(0, (it - epochs_.begin()) - 1));
}

double Demography::nextBoundary(double t) const {
    auto it = std::upper_bound(epochs_.begin(), epochs_.end(), t,
                               [](double v, const Epoch& e) { return v < e.start; });
    return it == epochs_.end() ? INF : it->start;
}

// Inverse-transform sampling: draw E ~ Exp(1) and find T with H(t0, T) = E, where
// H is piecewise linear with slope c/lambda_j in epoch j.
double Demography::sampleWaiting(double t0, double c, Rng& rng) const {
    double E = rng.exp1();
    double t = t0;
    size_t j = epochIndex(t);
    for (;;) {
        const double rate = c / epochs_[j].lambda;
        const double end = (j + 1 < epochs_.size()) ? epochs_[j + 1].start : INF;
        const double need = E / rate;
        if (t + need < end) return t + need;
        E -= rate * (end - t);
        t = end;
        ++j;
    }
}

double Demography::integratedHazard(double t0, double t1, double c) const {
    double H = 0.0, t = t0;
    size_t j = epochIndex(t);
    while (t < t1) {
        const double end = std::min(t1, (j + 1 < epochs_.size()) ? epochs_[j + 1].start : INF);
        H += c / epochs_[j].lambda * (end - t);
        t = end;
        ++j;
    }
    return H;
}

std::string Demography::describe(double N0, int ploidy) const {
    std::ostringstream os;
    const double unit = static_cast<double>(ploidy) * N0;
    for (const auto& e : epochs_) os << "  from " << e.start * unit << " gen: N=" << e.lambda * N0 << "\n";
    return os.str();
}

}  // namespace hapcoal
