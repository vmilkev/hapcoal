/**
 * @file GeneticMap.hpp
 * @brief Piecewise-linear genetic map (physical bp → genetic Morgans) of one chromosome.
 *
 * The map is a set of breakpoints x_0 < x_1 < … < x_m (bp) with cumulative
 * genetic positions G_0 ≤ G_1 ≤ … ≤ G_m (Morgans). Between breakpoints the
 * recombination rate is constant, r_i = (G_{i+1} − G_i)/(x_{i+1} − x_i) per bp
 * per generation, so G(x) is piecewise linear and invertible on segments with
 * r_i > 0. The simulator draws distances between recombination events on the
 * genetic scale and maps them back with toBp(), which honours hotspots of any
 * width exactly (zero-rate segments are skipped).
 */
#pragma once
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace hapcoal {

/// Genetic map of one chromosome region.
class GeneticMap {
public:
    GeneticMap() = default;

    /// Uniform rate @p ratePerBp (Morgans per bp per generation) on [start, end).
    static GeneticMap uniform(double start, double end, double ratePerBp);

    /**
     * @brief Map from explicit points.
     * @param posBp  increasing positions (bp); duplicates are dropped
     * @param cM     cumulative genetic positions (cM), non-decreasing
     * The first point is shifted to 0 M (only differences matter).
     */
    static GeneticMap fromPoints(const std::vector<double>& posBp, const std::vector<double>& cM);

    /// Parameters of a simulated recombination-hotspot landscape (withHotspots()).
    struct HotspotParams {
        double densityPerMb = 20.0;                  ///< hotspots per Mb
        double ampMin = 10.0, ampMax = 20.0;         ///< peak rate multiplier, uniform in [ampMin, ampMax]
        double widthMin = 800.0, widthMax = 1500.0;  ///< Gaussian s.d. of a hotspot (bp), uniform
        double centromereStart = -1.0;               ///< centromere start as a fraction of the length (<0: none)
        double centromereEnd = -1.0;                 ///< centromere end as a fraction of the length
        double centromereFactor = 0.1;               ///< rate multiplier inside the centromere
    };

    /**
     * @brief Background rate plus Gaussian hotspots, rescaled to a given mean rate.
     *
     * Hotspot centres are uniform outside the centromere; the relative rate is
     * 1 + Σ_h (A_h − 1)·exp(−(x − c_h)²/(2 s_h²)), evaluated on 13 steps across
     * ±3 s.d. of every hotspot, and the whole map is rescaled so that its mean
     * rate equals @p meanRatePerBp (total length = meanRatePerBp·(end − start) M).
     *
     * @param start, end      region (bp)
     * @param meanRatePerBp   mean recombination rate of the result (Morgans per bp per generation)
     * @param hp              hotspot landscape parameters
     * @param seed, stream    identify the random stream (reproducible)
     */
    static GeneticMap withHotspots(double start, double end, double meanRatePerBp, const HotspotParams& hp,
                                   uint64_t seed, uint64_t stream);

    /**
     * @brief Load a genetic map file (plain or .gz), one map per chromosome.
     *
     * Supported formats (auto-detected per line; header lines are skipped):
     *  - 4 columns `chr pos_bp rate_cM/Mb cM`  (HapMap II, Eagle)
     *  - 4 columns `chr id cM pos_bp`          (PLINK .map)
     *  - 3 columns `pos_bp rate_cM/Mb cM`      (1000 Genomes / IMPUTE; chromosome = @p defaultChr)
     * Positions are sorted; chromosome names are normalised (normaliseChr()).
     * @throws std::runtime_error if the file cannot be read or holds no usable entries
     */
    static std::map<std::string, GeneticMap> loadFile(const std::string& path, const std::string& defaultChr = "1");

    /// Chromosome name without a leading "chr"/"Chr"/"CHR" ("chr22" → "22").
    static std::string normaliseChr(std::string c);

    bool empty() const { return pos_.size() < 2; }
    double start() const { return pos_.front(); }                    ///< first breakpoint (bp)
    double end() const { return pos_.back(); }                       ///< last breakpoint (bp)
    double totalMorgans() const { return cum_.back() - cum_.front(); }  ///< map length (M)

    /// Genetic position (M) of physical position @p bp (clamped to the map range).
    double toMorgans(double bp) const;
    /// Physical position (bp) of genetic position @p morgans; on a zero-rate plateau the end of the plateau.
    double toBp(double morgans) const;
    /// Recombination rate (per bp per generation) at @p bp; 0 outside the map.
    double rateAt(double bp) const;

    /// The map restricted to [a, b] (end points interpolated, first point shifted to 0 M).
    GeneticMap clipped(double a, double b) const;

    const std::vector<double>& positions() const { return pos_; }   ///< breakpoints (bp)
    const std::vector<double>& cumulative() const { return cum_; }  ///< cumulative positions (M)

private:
    std::vector<double> pos_;  // bp
    std::vector<double> cum_;  // Morgans
};

}  // namespace hapcoal
