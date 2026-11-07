/**
 * @file MapEstimator.hpp
 * @brief Genetic maps reported with the simulated haplotypes.
 *
 *  - **True map**: the map used to simulate (input map or generated hotspot map).
 *  - **Realised map**: the rate re-estimated from the simulated ancestry. Along
 *    the sequence, recombination events occur as a Poisson process with
 *    intensity r(x)·p·N0·T(x), where T(x) is the total branch length of the local
 *    tree (coalescent units). Hence, in a window W,
 *        r̂_W = (number of events in W) / (p·N0·∫_W T(x) dx)
 *    is an unbiased estimator of the mean rate in W.
 *  - **LD-based map**: inferred from the simulated haplotypes alone. In each
 *    window the population-scaled rate ρ = 2·p·Ne·c per bp is fitted (least
 *    squares, pairs weighted by count) to binned mean r² versus distance d using
 *    the expectation of Hill & Weir (1988) with C = ρd and sample size n:
 *        E[r²] ≈ (10 + C)/((2 + C)(11 + C)) · [1 + (3 + C)(12 + 12C + C²)/(n(2 + C)(11 + C))].
 *    c is converted to cM/Mb with the supplied Ne, so the absolute scale is only
 *    right if Ne describes the relevant (recent) population size.
 */
#pragma once
#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

#include "hapcoal/GeneticMap.hpp"
#include "hapcoal/LDCalculator.hpp"

namespace hapcoal {

/// Realised recombination map: events and branch-length integral per window.
struct RealisedMap {
    double start = 0, end = 0, window = 1e5;  ///< region (bp) and window length (bp)
    std::vector<uint64_t> events;             ///< recombination events per window
    std::vector<double> integral;             ///< ∫ T(x) dx per window (coalescent units × bp)

    /// Set up empty windows of length @p w over [s, e).
    void init(double s, double e, double w);
    /// Accumulate total branch length @p T over [a, b).
    void addSpan(double a, double b, double T);
    /// Count a recombination event at position @p x.
    void addEvent(double x);
    /// Estimated rate (per bp per generation) of window @p k given the scale p·N0.
    double rate(size_t k, double scale) const { return integral[k] > 0 ? events[k] / (scale * integral[k]) : 0.0; }
};

/// Hill & Weir (1988) expected r² for scaled distance C = ρd and @p n haplotypes.
double hillWeirR2(double C, double n);

/// Write the true map as lines "chrom pos_bp rate_cM_per_Mb pos_cM".
void writeTrueMap(std::ostream& out, const std::string& chrom, const GeneticMap& map);
/// Write the realised map as lines "chrom start_bp end_bp events rate_cM_per_Mb pos_cM".
void writeRealisedMap(std::ostream& out, const std::string& chrom, const RealisedMap& rm, double scale);

/// One window of the LD-based map.
struct LDMapWindow {
    double start, end;   ///< window (bp)
    double rhoPerBp;     ///< fitted ρ = 2·p·Ne·c per bp (0 if fewer than 50 pairs)
    double cMperMb;      ///< c converted to cM/Mb with the given Ne
    uint64_t pairs;      ///< number of SNP pairs used
};

/**
 * @brief Fit ρ per window (grid search on log10 ρ, then golden-section refinement).
 * @param w          per-window r² bins from computeLD()
 * @param regionEnd  end of the region (bp)
 * @param nHap       number of haplotypes (sample size in the Hill–Weir formula)
 * @param Ne, ploidy population size and ploidy used to convert ρ to c
 */
std::vector<LDMapWindow> estimateLDMap(const LDWindows& w, double regionEnd, size_t nHap, double Ne, int ploidy);
/// Write the LD-based map as lines "chrom start_bp end_bp pairs rho_per_kb rate_cM_per_Mb pos_cM".
void writeLDMap(std::ostream& out, const std::string& chrom, const std::vector<LDMapWindow>& win);

}  // namespace hapcoal
