/**
 * @file Demography.hpp
 * @brief Piecewise-constant population size history N(t).
 *
 * Time is measured backwards from the present in **coalescent units** of
 * p·N0 generations (p = ploidy, N0 = reference size). In epoch j the relative
 * population size is λ_j = N_j / N0 and a set of k lineages coalesces at total
 * rate C(k,2)/λ_j. With a time-varying hazard h(t) = c/λ(t) the waiting time T
 * from t0 satisfies P(T > t) = exp(−∫_{t0}^{t} h(s) ds); it is sampled by
 * inverting the piecewise-linear cumulative hazard (sampleWaiting()).
 */
#pragma once
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace hapcoal {

class Rng;  // internal random number generator (src/Rng.hpp)

/// Piecewise-constant population size history, in coalescent time units.
class Demography {
public:
    /// One epoch: from coalescent time @ref start (inclusive) to the next epoch's start.
    struct Epoch {
        double start;   ///< start of the epoch (coalescent units, backwards in time)
        double lambda;  ///< relative size N(t)/N0 in the epoch
    };

    /// Constant population size (a single epoch with λ = 1).
    Demography();

    /**
     * @brief Build from (generations ago, N) pairs.
     * @param genN    epochs as (start in generations before present, diploid/ploid size N);
     *                unsorted input is sorted; an epoch at time 0 with size N0 is added if missing
     * @param N0      reference effective population size (defines the time unit)
     * @param ploidy  ploidy p; one coalescent unit = p·N0 generations
     * @throws std::invalid_argument for non-positive sizes, negative times, N0 or ploidy <= 0
     */
    static Demography fromGenerations(std::vector<std::pair<double, double>> genN, double N0, int ploidy);

    /// Read epochs from a whitespace-separated file "generations_ago  N" ('#' starts a comment).
    static Demography fromFile(const std::string& path, double N0, int ploidy);

    /**
     * @brief Tennessen et al. (2012) European (CEU) history as (generations, N) epochs.
     *
     * Model "OutOfAfrica_2T12" (as in stdpopsim; generation time 25 y,
     * intended mutation rate 2.36e-8): ancestral N=7310; African expansion to
     * 14474 at 148 kya; out-of-Africa bottleneck N=1861 at 51 kya; European
     * founding N=1032 at 23 kya growing at 0.307 %/generation, then at
     * 1.95 %/generation from 5.115 kya. Exponential growth is discretised into
     * 40 piecewise-constant epochs (size at each epoch midpoint).
     */
    static std::vector<std::pair<double, double>> tennessen2012CEU();

    /// All epochs, sorted by start time; the first starts at 0.
    const std::vector<Epoch>& epochs() const { return epochs_; }

    /// Relative size λ(t) at coalescent time @p t.
    double lambda(double t) const { return epochs_[epochIndex(t)].lambda; }
    /// Index of the epoch containing time @p t.
    size_t epochIndex(double t) const;
    /// First epoch boundary strictly greater than @p t (infinity if none).
    double nextBoundary(double t) const;

    /**
     * @brief Absolute time of the next event of a process with hazard c/λ(s), s ≥ t0.
     * @param t0   start time (coalescent units)
     * @param c    constant multiplicity of the hazard (e.g. C(k,2) for k lineages)
     * @param rng  random number generator (internal type)
     */
    double sampleWaiting(double t0, double c, Rng& rng) const;

    /// Integrated hazard ∫_{t0}^{t1} c/λ(s) ds (used in tests).
    double integratedHazard(double t0, double t1, double c) const;

    /// True for a constant size history (λ ≡ 1).
    bool isConstant() const { return epochs_.size() == 1 && epochs_[0].lambda == 1.0; }
    /// Human-readable list of epochs in generations and absolute sizes.
    std::string describe(double N0, int ploidy) const;

    static constexpr double INF = std::numeric_limits<double>::infinity();

private:
    std::vector<Epoch> epochs_;
};

}  // namespace hapcoal
