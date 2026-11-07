/**
 * @file DemographyFit.hpp
 * @brief Fit a population size history to an observed LD decay curve.
 *
 * LD between SNPs at distance d mainly reflects coalescence times of order
 * 1/(2c(d)) generations (Hayes et al. 2003), so the LD decay curve informs the
 * recent history of N. A parametric **bottleneck + growth** model with six
 * parameters θ = (N0, N1, T1, NB, ΔT, NA) is fitted:
 *
 *   N(g) = N0·(N1/N0)^(g/T1)          for 0 ≤ g < T1   (exponential growth, 16 epochs)
 *   N(g) = NB                         for T1 ≤ g < T1 + ΔT  (bottleneck)
 *   N(g) = NA                         for g ≥ T1 + ΔT  (ancestral)
 *
 * by minimising, with the Nelder–Mead simplex method on log θ,
 *   f(θ) = sqrt( (1/K) Σ_k [ log( r̄²_sim,k(θ) / r̄²_obs,k ) ]² ),
 * where r̄²_k is the mean r² in K = 9 distance classes with edges
 * 0, 1, 2, 5, 10, 20, 50, 100, 200, 500 kb (every class weighs the same).
 * Each evaluation is a full simulation with the same seed (common random
 * numbers), which makes f smooth enough for a derivative-free optimiser.
 */
#pragma once
#include <iosfwd>
#include <utility>
#include <vector>

#include "hapcoal/CoalescentSimulator.hpp"
#include "hapcoal/LDCalculator.hpp"

namespace hapcoal {

/// Epochs (generations ago, N) of the bottleneck + growth model for
/// log-parameters (log N0, log N1, log T1, log NB, log ΔT, log NA).
std::vector<std::pair<double, double>> bottleneckGrowthEpochs(const std::vector<double>& logParams);

/// Result of fitBottleneckGrowth().
struct DemographyFitResult {
    std::vector<std::pair<double, double>> epochs;  ///< best epochs (generations ago, N)
    double objective = 0;   ///< best f(θ) (rms log ratio of class means)
    double rmse = 0;        ///< RMSE of mean r² per 1-kb bin for the best fit
    int evaluations = 0;    ///< number of simulations run
};

/**
 * @brief Nelder–Mead fit of the bottleneck + growth model to an observed LD decay.
 *
 * @param observed  observed LD decay (e.g. observedLD())
 * @param base      simulation parameters (model, sample size, map, mu, seed, LD options);
 *                  demography, file output and haplotype storage are overridden
 * @param maxEvaluations  maximum number of simulations
 * @param log       progress messages for every improvement (may be null)
 *
 * The simplex starts near the Tennessen et al. (2012) European model
 * (N0 = 50000, N1 = 2000, T1 = 920, NB = 1861, ΔT = 1120, NA = 10000) with
 * initial steps of 0.7 on the log scale, and stops when the objective range of
 * the simplex is below 1e−3 or after @p maxEvaluations simulations.
 */
DemographyFitResult fitBottleneckGrowth(const LDDecay& observed, SimParams base, int maxEvaluations,
                                        std::ostream* log = nullptr);

}  // namespace hapcoal
