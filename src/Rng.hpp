/**
 * @file Rng.hpp
 * @brief (internal) portable, reproducible random number generation.
 *
 * std::*_distribution output is implementation-defined, so the same seed gives
 * different streams with libstdc++, libc++ and MSVC. We therefore use only the
 * raw std::mt19937_64 engine (fully specified by the standard) and implement
 * all distributions ourselves.
 */
#pragma once
#include <cmath>
#include <cstdint>
#include <random>

namespace hapcoal {

/**
 * @brief Reproducible random numbers: std::mt19937_64 plus own distributions.
 *
 * Every (seed, stream) pair defines an independent generator (seeded through
 * std::seed_seq), so each chromosome gets its own stream and results do not
 * depend on the number of threads.
 */
class Rng {
public:
    /// @param seed master seed; @param stream stream index (e.g. chromosome number)
    explicit Rng(uint64_t seed = 42, uint64_t stream = 0) { reseed(seed, stream); }

    /// Re-initialise the generator for (seed, stream).
    void reseed(uint64_t seed, uint64_t stream) {
        std::seed_seq ss{static_cast<uint32_t>(seed), static_cast<uint32_t>(seed >> 32),
                         static_cast<uint32_t>(stream), static_cast<uint32_t>(stream >> 32),
                         0x9E3779B9u, 0x7F4A7C15u};
        g_.seed(ss);
    }

    /// Raw 64-bit random number.
    uint64_t next() { return g_(); }

    /// Uniform on [0,1).
    double uniform() { return static_cast<double>(g_() >> 11) * 0x1.0p-53; }
    /// Uniform on (0,1].
    double uniformPos() { return static_cast<double>((g_() >> 11) + 1) * 0x1.0p-53; }
    /// Uniform on [a,b).
    double uniform(double a, double b) { return a + (b - a) * uniform(); }
    /// Standard exponential.
    /// Standard exponential, by inversion: -log(U), U in (0,1].
    double exp1() { return -std::log(uniformPos()); }
    /// Exponential with the given rate.
    double exponential(double rate) { return exp1() / rate; }

    /// Uniform integer in [0, n), unbiased (rejection of the incomplete last block).
    uint64_t below(uint64_t n) {
        const uint64_t threshold = (0 - n) % n;  // = 2^64 mod n
        for (;;) {
            const uint64_t x = g_();
            if (x >= threshold) return x % n;
        }
    }

    /// Poisson variate: multiplication method for small means, PTRS (Hörmann 1993) otherwise.
    int64_t poisson(double mean) {
        if (mean <= 0.0) return 0;
        if (mean < 10.0) {
            const double L = std::exp(-mean);
            int64_t k = 0;
            double p = uniformPos();
            while (p > L) { ++k; p *= uniformPos(); }
            return k;
        }
        const double slam = std::sqrt(mean), loglam = std::log(mean);
        const double b = 0.931 + 2.53 * slam;
        const double a = -0.059 + 0.02483 * b;
        const double invalpha = 1.1239 + 1.1328 / (b - 3.4);
        const double vr = 0.9277 - 3.6224 / (b - 2.0);
        for (;;) {
            const double U = uniform() - 0.5;
            const double V = uniform();
            const double us = 0.5 - std::fabs(U);
            const double kd = std::floor((2.0 * a / us + b) * U + mean + 0.43);
            if (us >= 0.07 && V <= vr) return static_cast<int64_t>(kd);
            if (kd < 0.0 || (us < 0.013 && V > us)) continue;
            if (std::log(V) + std::log(invalpha) - std::log(a / (us * us) + b) <=
                -mean + kd * loglam - std::lgamma(kd + 1.0))
                return static_cast<int64_t>(kd);
        }
    }

private:
    std::mt19937_64 g_;
};

}  // namespace hapcoal
