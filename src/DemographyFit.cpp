/**
 * @file DemographyFit.cpp
 * @brief Nelder–Mead fit of a bottleneck + growth demography to observed LD decay.
 */
#include "hapcoal/DemographyFit.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ostream>

namespace hapcoal {

std::vector<std::pair<double, double>> bottleneckGrowthEpochs(const std::vector<double>& lp) {
    const double N0 = std::exp(lp[0]), N1 = std::exp(lp[1]), T1 = std::exp(lp[2]);
    const double NB = std::exp(lp[3]), T2 = T1 + std::exp(lp[4]), NA = std::exp(lp[5]);
    std::vector<std::pair<double, double>> e;
    // Exponential growth from N1 (at T1 generations ago) to N0 (present), as 16
    // piecewise-constant epochs with the size at each epoch's midpoint.
    const int k = 16;
    for (int i = 0; i < k; ++i) {
        const double g0 = T1 * i / k, gm = T1 * (i + 0.5) / k;
        e.emplace_back(g0, N0 * std::pow(N1 / N0, gm / T1));
    }
    e.emplace_back(T1, NB);  // bottleneck
    e.emplace_back(T2, NA);  // ancestral
    return e;
}

DemographyFitResult fitBottleneckGrowth(const LDDecay& obs, SimParams base, int maxEval, std::ostream* log) {
    // Distance classes (bp): every class contributes equally to the objective.
    const std::vector<double> edges = {0, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5};
    base.log = nullptr;
    base.estimateMap = false;
    base.demographyFile.clear();
    base.outPrefix.clear();     // no files: SNP rows stay in memory
    base.keepHaplotypes = false;
    base.ldDecay = true;
    DemographyFitResult res;
    double bestObj = 1e300;

    // Objective f(theta): one simulation, rms of log(mean r2 sim / mean r2 obs) over classes.
    auto f = [&](const std::vector<double>& lp) {
        base.epochs = bottleneckGrowthEpochs(lp);
        CoalescentSimulator sim(base);
        sim.run();
        double obj = 0;
        for (size_t k = 0; k + 1 < edges.size(); ++k) {
            const double rs = sim.decay().meanR2(edges[k], edges[k + 1]);
            const double ro = obs.meanR2(edges[k], edges[k + 1]);
            const double l = (rs > 0 && ro > 0) ? std::log(rs / ro) : 5.0;
            obj += l * l;
        }
        obj = std::sqrt(obj / (edges.size() - 1));
        ++res.evaluations;
        if (obj < bestObj) {
            bestObj = obj;
            res.objective = obj;
            res.epochs = base.epochs;
            res.rmse = decayRmse(sim.decay(), obs, edges.back());
            if (log)
                *log << "  eval " << res.evaluations << ": rms log ratio " << std::setprecision(4) << obj
                     << ", RMSE(r2) " << res.rmse << "  N0=" << std::setprecision(0) << std::fixed << std::exp(lp[0])
                     << " N1=" << std::exp(lp[1]) << " T1=" << std::exp(lp[2]) << " NB=" << std::exp(lp[3])
                     << " T2=" << std::exp(lp[2]) + std::exp(lp[4]) << " NA=" << std::exp(lp[5]) << std::defaultfloat
                     << std::setprecision(6) << "\n";
        }
        return obj;
    };

    // ---- Nelder-Mead simplex on log parameters (reflection 1, expansion 2,
    //      contraction 0.5, shrink 0.5). Start near Tennessen et al. (2012).
    std::vector<double> x0 = {std::log(50000.0), std::log(2000.0), std::log(920.0), std::log(1861.0),
                              std::log(1120.0), std::log(10000.0)};
    const size_t D = x0.size();
    std::vector<std::vector<double>> S(D + 1, x0);
    std::vector<double> F(D + 1);
    for (size_t i = 0; i < D; ++i) S[i + 1][i] += 0.7;
    for (size_t i = 0; i <= D; ++i) F[i] = f(S[i]);
    while (res.evaluations < maxEval) {
        // order vertices by objective
        std::vector<size_t> idx(D + 1);
        for (size_t i = 0; i <= D; ++i) idx[i] = i;
        std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return F[a] < F[b]; });
        std::vector<std::vector<double>> S2;
        std::vector<double> F2;
        for (size_t i : idx) {
            S2.push_back(S[i]);
            F2.push_back(F[i]);
        }
        S.swap(S2);
        F.swap(F2);
        if (F[D] - F[0] < 1e-3) break;  // converged
        // centroid of all but the worst vertex
        std::vector<double> c(D, 0.0);
        for (size_t i = 0; i < D; ++i)
            for (size_t j = 0; j < D; ++j) c[j] += S[i][j] / D;
        auto lin = [&](double a) {  // c + a (worst - c)
            std::vector<double> x(D);
            for (size_t j = 0; j < D; ++j) x[j] = c[j] + a * (S[D][j] - c[j]);
            return x;
        };
        const auto xr = lin(-1.0);  // reflection
        const double fr = f(xr);
        if (fr < F[0]) {
            const auto xe = lin(-2.0);  // expansion
            const double fe = f(xe);
            if (fe < fr) {
                S[D] = xe;
                F[D] = fe;
            } else {
                S[D] = xr;
                F[D] = fr;
            }
        } else if (fr < F[D - 1]) {
            S[D] = xr;
            F[D] = fr;
        } else {
            const auto xc = lin(fr < F[D] ? -0.5 : 0.5);  // outside / inside contraction
            const double fc = f(xc);
            if (fc < std::min(fr, F[D])) {
                S[D] = xc;
                F[D] = fc;
            } else {  // shrink towards the best vertex
                for (size_t i = 1; i <= D; ++i) {
                    for (size_t j = 0; j < D; ++j) S[i][j] = S[0][j] + 0.5 * (S[i][j] - S[0][j]);
                    F[i] = f(S[i]);
                }
            }
        }
    }
    return res;
}

}  // namespace hapcoal
