// test_main.cpp — unit and statistical tests for hapcoal.
//
// Statistical tests compare Monte-Carlo estimates with exact coalescent results
// and use tolerances of ~4-5 standard errors, so they are deterministic (fixed
// seeds) and should practically never fail for a correct implementation.
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <hapcoal/hapcoal.hpp>

#include "HaplotypeStore.hpp"  // internal components are tested directly
#include "LocalTree.hpp"
#include "MacsHistory.hpp"
#include "Rng.hpp"

using namespace hapcoal;
namespace fs = std::filesystem;

static int g_fail = 0, g_pass = 0;
#define CHECK(cond, msg)                                                         \
    do {                                                                         \
        if (cond) ++g_pass;                                                      \
        else { ++g_fail; std::cout << "  FAIL: " << msg << "  [" #cond "]\n"; } \
    } while (0)
#define CHECK_NEAR(val, expect, tol, msg)                                                       \
    do {                                                                                        \
        const double v_ = (val), e_ = (expect), t_ = (tol);                                     \
        if (std::fabs(v_ - e_) <= t_) ++g_pass;                                                 \
        else { ++g_fail; std::cout << "  FAIL: " << msg << ": got " << v_ << ", expected " << e_ \
                                   << " +- " << t_ << "\n"; }                                   \
        std::cout << "    " << msg << ": " << v_ << " (expected " << e_ << " +- " << t_ << ")\n"; \
    } while (0)

static double harmonic(int n) {
    double h = 0;
    for (int i = 1; i <= n; ++i) h += 1.0 / i;
    return h;
}

// Test files go to out/tests (relative to the working directory, i.e. the project
// root for `make test`), so `make rclean` removes them. HAPCOAL_TEST_DIR overrides.
static fs::path tmpDir() {
    const char* env = std::getenv("HAPCOAL_TEST_DIR");
    fs::path d = env && *env ? fs::path(env) : fs::path("out") / "tests";
    fs::create_directories(d);
    return d;
}

// =================================================================== Step 1
static void testRng() {
    std::cout << "[rng]\n";
    Rng r(7);
    const int N = 200000;
    double s = 0;
    for (int i = 0; i < N; ++i) s += r.exp1();
    CHECK_NEAR(s / N, 1.0, 5 * 1.0 / std::sqrt(N), "Exp(1) mean");
    for (double mean : {0.3, 4.0, 37.5, 2500.0}) {
        double m = 0, v = 0;
        const int M = 100000;
        std::vector<double> xs(M);
        for (int i = 0; i < M; ++i) { xs[i] = static_cast<double>(r.poisson(mean)); m += xs[i]; }
        m /= M;
        for (double x : xs) v += (x - m) * (x - m);
        v /= (M - 1);
        CHECK_NEAR(m, mean, 5 * std::sqrt(mean / M), "Poisson(" + std::to_string(mean) + ") mean");
        CHECK_NEAR(v / mean, 1.0, 0.03, "Poisson(" + std::to_string(mean) + ") var/mean");
    }
    std::vector<int> c(7, 0);
    for (int i = 0; i < 70000; ++i) ++c[r.below(7)];
    int mx = 0, mn = 1 << 30;
    for (int x : c) { mx = std::max(mx, x); mn = std::min(mn, x); }
    CHECK(mx - mn < 600, "below(7) roughly uniform");
    Rng a(1, 2), b(1, 2), d(1, 3);
    CHECK(a.next() == b.next(), "same seed/stream reproducible");
    CHECK(a.next() != d.next(), "different streams differ");
}

static void testDemography() {
    std::cout << "[demography]\n";
    Rng r(11);
    Demography constant;
    double s = 0;
    const int N = 200000;
    for (int i = 0; i < N; ++i) s += constant.sampleWaiting(0.0, 3.0, r);
    CHECK_NEAR(s / N, 1.0 / 3.0, 5 * (1.0 / 3.0) / std::sqrt(N), "constant: mean waiting time 1/c");

    // Bottleneck: N = 10000 now, 1000 between 400 and 800 generations, 20000 before.
    const Demography d = Demography::fromGenerations({{0, 10000}, {400, 1000}, {800, 20000}}, 10000, 2);
    CHECK(d.epochs().size() == 3, "three epochs");
    CHECK_NEAR(d.epochs()[1].start, 400.0 / 20000.0, 1e-15, "epoch start in coalescent units");
    CHECK_NEAR(d.lambda(0.03), 0.1, 1e-15, "lambda in bottleneck");
    // Survival function P(T > t) = exp(-H(t)) for pairwise coalescence (c = 1).
    std::vector<double> ts = {0.01, 0.02, 0.03, 0.05, 0.5, 1.5};
    std::vector<int> above(ts.size(), 0);
    for (int i = 0; i < N; ++i) {
        const double T = d.sampleWaiting(0.0, 1.0, r);
        for (size_t k = 0; k < ts.size(); ++k) above[k] += T > ts[k];
    }
    for (size_t k = 0; k < ts.size(); ++k) {
        const double p = std::exp(-d.integratedHazard(0.0, ts[k], 1.0));
        CHECK_NEAR(above[k] / double(N), p, 5 * std::sqrt(p * (1 - p) / N) + 1e-4,
                   "bottleneck P(T2 > " + std::to_string(ts[k]) + ")");
    }
}

static void testGeneticMap() {
    std::cout << "[genetic map]\n";
    const GeneticMap u = GeneticMap::uniform(1000, 2001000, 1e-8);
    CHECK_NEAR(u.totalMorgans(), 0.02, 1e-15, "uniform map length");
    CHECK_NEAR(u.toBp(u.toMorgans(123456.7)), 123456.7, 1e-6, "uniform bp->M->bp");
    const GeneticMap m = GeneticMap::fromPoints({0, 100, 200, 300, 400}, {0, 0, 1, 1, 3});
    CHECK_NEAR(m.toMorgans(150), 0.005, 1e-15, "interpolation");
    CHECK_NEAR(m.toBp(0.015), 325.0, 1e-9, "inverse interpolation");
    CHECK_NEAR(m.toBp(0.01), 300.0, 1e-9, "inverse skips a zero-rate plateau (no recombination inside)");
    CHECK_NEAR(m.rateAt(350), 0.02 / 100, 1e-15, "rate in segment");
    Rng r(3);
    GeneticMap::HotspotParams hp;
    hp.centromereStart = 0.45;
    hp.centromereEnd = 0.55;
    const GeneticMap h = GeneticMap::withHotspots(0, 5e6, 1.2e-8, hp, 3, 0);
    CHECK_NEAR(h.totalMorgans(), 0.06, 1e-9, "hotspot map rescaled to mean rate");
    double mx = 0;
    for (double x = 0; x < 5e6; x += 100) mx = std::max(mx, h.rateAt(x));
    CHECK(mx > 4 * 1.2e-8, "hotspots present (peak rate > 4x mean)");
    CHECK(h.rateAt(2.5e6) < 1.2e-8, "reduced rate in centromere");

    const fs::path d = tmpDir();
    {
        std::ofstream f(d / "hapmap.txt");
        f << "Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)\nchr22\t16051347\t8.1\t0\nchr22\t16052618\t0.1\t0.0103\n"
             "chr22\t17000000\t1.0\t0.1\nchr21\t100\t1\t0\nchr21\t200\t1\t0.5\n";
    }
    auto maps = GeneticMap::loadFile((d / "hapmap.txt").string());
    CHECK(maps.size() == 2 && maps.count("22") && maps.count("21"), "HapMap format: two chromosomes, chr prefix stripped");
    CHECK_NEAR(maps.at("22").totalMorgans(), 0.001, 1e-12, "HapMap total length");
    {
        std::ofstream f(d / "plink.map");
        f << "22 rs1 0.0 1000\n22 rs2 1.5 2000\n22 rs3 2.0 5000\n";
    }
    auto pm = GeneticMap::loadFile((d / "plink.map").string());
    CHECK_NEAR(pm.at("22").toMorgans(1500), 0.0075, 1e-12, "PLINK map format");
    {
        std::ofstream f(d / "impute.txt");
        f << "position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n100 1 0\n1100 1 0.001\n";
    }
    auto im = GeneticMap::loadFile((d / "impute.txt").string(), "7");
    CHECK(im.count("7") == 1, "3-column map assigned to default chromosome");
    const GeneticMap c = maps.at("22").clipped(16100000, 16500000);
    CHECK_NEAR(c.totalMorgans(), maps.at("22").toMorgans(16500000) - maps.at("22").toMorgans(16100000), 1e-15, "clipped map");
}

// =================================================================== Step 2
static void testInitialTree() {
    std::cout << "[local tree, no recombination]\n";
    Rng r(5);
    const int n = 10, R = 20000;
    Demography constant;
    LocalTree t(n);
    double sT = 0, sL = 0;
    std::vector<double> xi(n, 0.0);  // branch length subtending i leaves
    std::vector<int> leaves;
    bool ok = true;
    for (int rep = 0; rep < R; ++rep) {
        t.buildInitial(constant, r);
        if (rep < 50 && !t.validate()) ok = false;
        sT += t.time(t.root());
        sL += t.totalLength();
        for (int v = 0; v < t.numNodes(); ++v) {
            if (v == t.root()) continue;
            leaves.clear();
            t.collectLeaves(v, leaves);
            xi[leaves.size()] += t.time(t.parent(v)) - t.time(v);
        }
    }
    CHECK(ok, "initial trees valid");
    const double eT = 2.0 * (1.0 - 1.0 / n), eL = 2.0 * harmonic(n - 1);
    CHECK_NEAR(sT / R, eT, 5 * 1.0 / std::sqrt(R), "E[TMRCA] = 2(1-1/n)");
    CHECK_NEAR(sL / R, eL, 5 * 2.0 / std::sqrt(R), "E[total length] = 2 H(n-1)");
    for (int i : {1, 2, 5, 9}) CHECK_NEAR(xi[i] / R, 2.0 / i, 0.05 * 2.0 / i + 0.02, "E[xi_" + std::to_string(i) + "] = 2/i (SFS)");

    // Bottleneck: pairwise coalescence time distribution under N(t)
    const Demography d = Demography::fromGenerations({{0, 10000}, {400, 1000}, {800, 20000}}, 10000, 2);
    LocalTree t2(2);
    int above = 0;
    const int R2 = 100000;
    for (int rep = 0; rep < R2; ++rep) {
        t2.buildInitial(d, r);
        above += t2.time(t2.root()) > 0.03;
    }
    const double p = std::exp(-d.integratedHazard(0, 0.03, 1.0));
    CHECK_NEAR(above / double(R2), p, 5 * std::sqrt(p * (1 - p) / R2), "bottleneck: P(T2 > 0.03)");
}

static void testLeafCache() {
    std::cout << "[leaf-set cache]\n";
    Rng r(17);
    const int n = 150;
    Demography d;
    LocalTree t(n);
    CHECK(t.enableLeafCache(1 << 30), "cache enabled");
    t.buildInitial(d, r);
    const size_t W = wordsFor(n);
    std::vector<uint64_t> row(W);
    std::vector<int> leaves;
    bool ok = true;
    for (int ev = 0; ev < 3000 && ok; ++ev) {
        // random SMC'-like regraft
        const int b = t.branchAt(r.uniform() * t.totalLength());
        const int op = t.parent(b), s = t.sibling(b);
        const double tr = r.uniform(t.time(b), t.time(op));
        const double tc = tr + r.exp1();
        int v = t.pickBranchCrossing(tc, b, r);
        if (v == op) v = s;
        t.regraft(b, v, tc);
        if (ev % 10 == 0) {
            for (int q = 0; q < t.numNodes(); ++q) {
                std::fill(row.begin(), row.end(), 0);
                const uint32_t c = t.leafBits(q, row.data());
                leaves.clear();
                t.collectLeaves(q, leaves);
                std::vector<uint64_t> ref(W, 0);
                for (int l : leaves) ref[l >> 6] |= 1ULL << (l & 63);
                ok = ok && c == leaves.size() && ref == row;
            }
        }
    }
    CHECK(ok, "cached leaf sets equal subtree traversal after 3000 regrafts");
}

// =================================================================== Step 3/4
// Collects T(x) and validates trees while a chromosome is simulated.
struct Probe : CoalescentSimulator::TreeSpanCallback {
    std::vector<double> at;      // positions to sample TMRCA (n=2) / total length
    std::vector<double> value;
    size_t next = 0;
    double weightedT = 0, weightedTMRCA = 0, len = 0;
    bool validateTrees = false, ok = true;
    void span(double a, double b, const LocalTree& t) override {
        weightedT += t.totalLength() * (b - a);
        weightedTMRCA += t.time(t.root()) * (b - a);
        len += b - a;
        while (next < at.size() && at[next] < b) {
            if (at[next] >= a) value.push_back(t.time(t.root()));
            ++next;
        }
        if (validateTrees && !t.validate()) ok = false;
    }
};

static ChromStats simulateOnce(SimParams p, double L, Probe& probe, uint64_t stream, size_t* nSnps = nullptr) {
    p.chromosomes = {{"1", 0, L, -1}};
    CoalescentSimulator sim(p);
    ChromosomeSpec spec = p.chromosomes[0];
    const GeneticMap map = sim.buildMap(spec, stream);
    SnpBuffer w(static_cast<size_t>(p.numHaplotypes()));  // in memory
    ChromStats st = sim.simulateChromosome(spec, map, stream, w, nullptr, &probe);
    w.close();
    if (nSnps) *nSnps = w.numSnps();
    return st;
}

static void testStationarity(Model model, double history = 1e5) {
    std::cout << "[" << modelName(model) << (model == Model::MaCS ? " h=" + std::to_string(int(history)) : "")
              << ": stationarity, validity, recombination count]\n";
    SimParams p;
    p.model = model;
    p.nHap = 20;
    p.Ne = 10000;
    p.rho = 1e-8;
    p.mu = 0;
    p.history = history;
    p.seed = 99;
    const double L = 5e6;  // rho = 4 N r L = 2000
    double sumT = 0, sumM = 0, len = 0, rec = 0;
    bool ok = true;
    const int R = 8;
    for (int rep = 0; rep < R; ++rep) {
        Probe pr;
        pr.validateTrees = rep == 0;
        const ChromStats st = simulateOnce(p, L, pr, 1000 + rep);
        ok = ok && pr.ok;
        sumT += pr.weightedT;
        sumM += pr.weightedTMRCA;
        len += pr.len;
        rec += st.recombinations;
    }
    const int n = p.nHap;
    const double eL = 2 * harmonic(n - 1), eM = 2 * (1 - 1.0 / n);
    CHECK(ok, "all local trees valid");
    // Sequence averages over ~R*2000*... nearly independent trees: allow 4% error
    CHECK_NEAR(sumT / len, eL, 0.04 * eL, "sequence-averaged total length = 2 H(n-1) (Kingman marginal)");
    CHECK_NEAR(sumM / len, eM, 0.05 * eM, "sequence-averaged TMRCA = 2(1-1/n)");
    const double eRec = 1e-8 * 2 * p.Ne * L * eL * R;
    CHECK_NEAR(rec, eRec, 0.04 * eRec, "number of recombinations = r*2N*L*E[T_total]");
}

// Covariance of pairwise coalescence times at distance d (n = 2).
// ARG (Griffiths 1981): Cov = (rho+18)/(rho^2+13 rho+18), rho = 4 N r d.
static double argCov(double rho) { return (rho + 18) / (rho * rho + 13 * rho + 18); }

// Reference: the n = 2 SMC / SMC' Markov chain written out directly (independent
// of LocalTree/MacsHistory). Distance in rho units; events at rate T.
static std::vector<double> chainCov(bool smcPrime, const std::vector<double>& rhos, int R, uint64_t seed) {
    Rng g(seed, 7);
    std::vector<double> s0(rhos.size(), 0), s1(rhos.size(), 0), s01(rhos.size(), 0);
    for (int r = 0; r < R; ++r) {
        double T = g.exp1(), x = 0;
        const double T0 = T;
        size_t k = 0;
        for (;;) {
            const double d = g.exp1() / T;
            while (k < rhos.size() && x + d > rhos[k]) { s0[k] += T0; s1[k] += T; s01[k] += T0 * T; ++k; }
            if (k == rhos.size()) break;
            x += d;
            const double tr = g.uniform() * T;
            if (smcPrime) {
                const double w = g.exp1() / 2;  // own branch + sister below T
                if (tr + w < T) { if (g.uniform() >= 0.5) T = tr + w; }
                else T = T + g.exp1();
            } else {
                T = tr + g.exp1();               // sister only
            }
        }
    }
    std::vector<double> c(rhos.size());
    for (size_t k = 0; k < rhos.size(); ++k) c[k] = s01[k] / R - (s0[k] / R) * (s1[k] / R);
    return c;
}

enum class CovRef { SMCPrimeChain, SMCChain, ARG };

static void testTwoLocusCovariance(Model model, double history, CovRef ref) {
    std::cout << "[" << modelName(model) << (model == Model::MaCS ? " h=" + std::to_string(int(history)) : "")
              << ": two-locus covariance of T2]\n";
    SimParams p;
    p.model = model;
    p.nHap = 2;
    p.Ne = 10000;
    p.rho = 1e-8;  // 4 N r = 4e-4 per bp
    p.mu = 0;
    p.history = history;
    p.seed = 1234;
    const std::vector<double> rhos = {0.5, 2.0, 8.0};
    const double perBp = 4 * p.Ne * p.rho;
    const int R = 40000;
    std::vector<double> s0(rhos.size(), 0), s1(rhos.size(), 0), s01(rhos.size(), 0), sab2(rhos.size(), 0);
    for (int rep = 0; rep < R; ++rep) {
        Probe pr;
        pr.at.push_back(0.0);
        for (double rr : rhos) pr.at.push_back(rr / perBp);
        simulateOnce(p, rhos.back() / perBp + 1, pr, 50000 + rep);
        for (size_t k = 0; k < rhos.size(); ++k) {
            const double a = pr.value[0], b = pr.value[k + 1];
            s0[k] += a; s1[k] += b; s01[k] += a * b; sab2[k] += (a - 1) * (a - 1) * (b - 1) * (b - 1);
        }
    }
    const int RC = 2000000;
    std::vector<double> refv;
    std::string refName;
    switch (ref) {
        case CovRef::SMCPrimeChain: refv = chainCov(true, rhos, RC, 31); refName = "exact SMC' chain"; break;
        case CovRef::SMCChain: refv = chainCov(false, rhos, RC, 32); refName = "exact SMC chain"; break;
        case CovRef::ARG: for (double r : rhos) refv.push_back(argCov(r)); refName = "ARG"; break;
    }
    for (size_t k = 0; k < rhos.size(); ++k) {
        const double m0 = s0[k] / R, m1 = s1[k] / R;
        const double cov = s01[k] / R - m0 * m1;
        const double se = std::sqrt(sab2[k] / R / R) + (ref == CovRef::ARG ? 0.0 : std::sqrt(8.0 / RC));
        const double tol = (ref == CovRef::ARG ? 0.035 : 0.0) + 4.5 * se;  // MaCS with finite h is still approximate
        CHECK_NEAR(cov, refv[k], tol,
                   "Cov(T(0),T(d)) at rho=" + std::to_string(rhos[k]).substr(0, 4) + " vs " + refName +
                       " (ARG " + std::to_string(argCov(rhos[k])).substr(0, 6) + ")");
        if (k == 0) CHECK_NEAR(m0, 1.0, 5 * 1.0 / std::sqrt(R), "E[T2] = 1");
    }
}

static void testMacsArg() {
    std::cout << "[MaCS: ARG bookkeeping]\n";
    // Validate ARG and local tree after every event for a small sample.
    Rng rng(77);
    const int n = 8;
    LocalTree tree(n);
    Demography d;
    tree.buildInitial(d, rng);
    MacsHistory arg;
    arg.init(tree);
    bool ok = true;
    std::string err;
    int nonLocal = 0;
    double pos = 0;
    for (int ev = 0; ev < 20000 && ok; ++ev) {
        pos += rng.exp1() * 100;
        arg.prune(pos - 2000);
        const int b = tree.branchAt(rng.uniform() * tree.totalLength());
        const int op = tree.parent(b), s = tree.sibling(b);
        const double tr = rng.uniform(tree.time(b), tree.time(op));
        const auto r = arg.recombine(tree, b, op, s, tr, pos, d, rng);
        nonLocal += r.nonLocal;
        tree.regraft(b, r.u, r.t);
        if (!tree.validate(&err)) { ok = false; std::cout << "    tree: " << err; }
        if (!arg.validate(tree, &err)) { ok = false; std::cout << "    arg: " << err; }
    }
    CHECK(ok, "ARG and local tree consistent after 20000 events");
    CHECK(nonLocal > 0, "non-ancestral coalescences occur with h > 0");
    std::cout << "    non-ancestral coalescences: " << nonLocal << ", retained edges: " << arg.numEdges() << "\n";
    // With h = 0 the ARG never grows beyond the local tree (+ the current event).
    MacsHistory a0;
    tree.buildInitial(d, rng);
    a0.init(tree);
    size_t maxE = 0;
    for (int ev = 0; ev < 5000; ++ev) {
        pos += 1;
        a0.prune(pos);
        maxE = std::max(maxE, a0.numEdges());
        const int b = tree.branchAt(rng.uniform() * tree.totalLength());
        const int op = tree.parent(b), s = tree.sibling(b);
        const auto r = a0.recombine(tree, b, op, s, rng.uniform(tree.time(b), tree.time(op)), pos, d, rng);
        tree.regraft(b, r.u, r.t);
    }
    CHECK(maxE <= static_cast<size_t>(2 * n - 1 + 4), "h=0: retained graph = local tree");
}

// =================================================================== Step 5
static void testTransposeAndStore() {
    std::cout << "[haplotype store]\n";
    // Two chromosomes with random rows; assemble; read back each chromosome.
    Rng r(9);
    const size_t nHap = 130;
    const fs::path d = tmpDir();
    std::vector<ChromData> chroms(2);
    std::vector<std::vector<std::vector<uint8_t>>> truth(2);
    const size_t nS[2] = {1000, 77};
    for (int c = 0; c < 2; ++c) {
        chroms[c].name = std::to_string(c + 1);
        chroms[c].start = 0;
        chroms[c].end = 1e6;
        chroms[c].map = GeneticMap::uniform(0, 1e6, 1e-8);
        // chromosome 1 in a temporary file, chromosome 2 in memory: both storage modes
        chroms[c].snps = std::make_unique<SnpBuffer>(nHap, c == 0 ? (d / "t0.snp").string() : std::string());
        SnpBuffer& w = *chroms[c].snps;
        truth[c].assign(nS[c], std::vector<uint8_t>(nHap));
        for (size_t j = 0; j < nS[c]; ++j) {
            uint64_t* row = w.beginRow();
            uint32_t cnt = 0;
            for (size_t h = 0; h < nHap; ++h)
                if (r.uniform() < 0.3) { row[h / 64] |= 1ULL << (h % 64); truth[c][j][h] = 1; ++cnt; }
            w.commitRow(static_cast<double>(j) * 10, cnt);
        }
        w.close();
        chroms[c].pos = w.positions();
        chroms[c].count = w.counts();
    }
    const std::string prefix = (d / "store").string();
    assembleHaplotypes(prefix, chroms, nHap, 2, true, 2);
    bool ok = true;
    for (int c = 0; c < 2; ++c) {
        const HaplotypeMatrix M = readChromosomeFromBin(prefix + ".hap.bin", chroms[c].rowWordOffset, nS[c]);
        ok = ok && M.nHap == nHap && M.nSnps == nS[c];
        for (size_t h = 0; h < nHap; ++h)
            for (size_t j = 0; j < nS[c]; ++j) ok = ok && M.allele(h, j) == truth[c][j][h];
    }
    CHECK(ok, "binary haplotype-major file reproduces SNP-major input (both chromosomes)");
    bool mok = true;
    for (int c = 0; c < 2; ++c) {
        const HaplotypeMatrix M = toHaplotypeMatrix(*chroms[c].snps, nHap, 3);
        for (size_t h = 0; h < nHap; ++h)
            for (size_t j = 0; j < nS[c]; ++j) mok = mok && M.allele(h, j) == truth[c][j][h];
    }
    CHECK(mok, "in-memory haplotype matrix (file and memory buffers) reproduces the input");
    std::ifstream txt(prefix + ".hap.txt");
    std::string line;
    bool tok = true;
    for (size_t h = 0; h < nHap && std::getline(txt, line); ++h) {
        tok = tok && line.size() == nS[0] + nS[1];
        for (size_t j = 0; j < nS[0]; ++j) tok = tok && (line[j] == '1') == (truth[0][j][h] == 1);
        for (size_t j = 0; j < nS[1]; ++j) tok = tok && (line[nS[0] + j] == '1') == (truth[1][j][h] == 1);
    }
    CHECK(tok, "text haplotypes: chromosomes consecutive per haplotype");
    std::vector<size_t> sel = {0, 5, 999};
    const BitRows br = loadRows(*chroms[0].snps, nHap, sel, chroms[0].pos, chroms[0].count, 0);
    bool lok = br.size() == 3;
    for (size_t k = 0; k < sel.size() && lok; ++k)
        for (size_t h = 0; h < nHap; ++h) lok = lok && (((br.row(k)[h / 64] >> (h % 64)) & 1) == truth[0][sel[k]][h]);
    CHECK(lok, "loadRows selects the right rows");
}

// ============================================================ library API
static void testLibraryApi() {
    std::cout << "[library API: memory / files / streaming / threads give identical results]\n";
    SimParams p;
    p.model = Model::MaCS;
    p.nInd = 60;
    p.chromosomes = {{"1", 0, 2e6, -1}, {"2", 0, 1.5e6, 1.3e-8}};
    p.epochs = {{0, 20000}, {300, 2000}, {1500, 10000}};
    p.hotspots = true;
    p.seed = 4;
    p.threads = 4;
    const SimulationResult mem = CoalescentSimulator(p).run();
    CHECK(mem.chromosomes.size() == 2 && mem.nHap == 120 && !mem.chromosomes[0].haplotypes.empty(),
          "in-memory result holds both chromosomes with haplotypes");
    // files
    SimParams pf = p;
    pf.outPrefix = (tmpDir() / "api").string();
    pf.threads = 1;
    const SimulationResult fil = CoalescentSimulator(pf).run();
    bool same = true;
    for (size_t c = 0; c < 2; ++c) {
        const auto& a = mem.chromosomes[c];
        const HaplotypeMatrix f = readChromosomeFromBin(pf.outPrefix + ".hap.bin", mem.rowWordOffset(c), a.positions.size());
        same = same && f.bits == a.haplotypes.bits && fil.chromosomes[c].haplotypes.bits == a.haplotypes.bits;
        same = same && fil.chromosomes[c].positions == a.positions;
    }
    CHECK(same, "file output == in-memory haplotypes (and 1 thread == 4 threads)");
    // streaming through temporary files
    SimParams ps = p;
    ps.tempDir = (tmpDir() / "api_tmp").string();
    std::vector<std::vector<uint64_t>> got;
    const SimulationResult str = CoalescentSimulator(ps).run([&](ChromosomeResult& c) { got.push_back(c.haplotypes.bits); });
    CHECK(got.size() == 2 && got[0] == mem.chromosomes[0].haplotypes.bits && got[1] == mem.chromosomes[1].haplotypes.bits,
          "callback receives the same haplotypes, in order");
    CHECK(str.chromosomes[0].haplotypes.empty() && str.chromosomes[0].positions == mem.chromosomes[0].positions,
          "after the callback haplotypes are released, other data kept");
    bool noTmp = true;
    for (const auto& e : fs::directory_iterator(ps.tempDir)) { (void)e; noTmp = false; }
    CHECK(noTmp, "temporary SNP files removed");
    // whole-genome row = chromosome rows concatenated
    const auto g = mem.genomeRow(5);
    const auto& c1 = mem.chromosomes[1].haplotypes;
    bool gok = g.size() == mem.rowWordOffset(2);
    for (size_t w = 0; w < c1.words && gok; ++w) gok = g[mem.rowWordOffset(1) + w] == c1.row(5)[w];
    CHECK(gok, "genomeRow concatenates chromosomes with the .hap.bin layout");
    CHECK(std::fabs(mem.ldDecay.meanR2(0, 1e3) - fil.ldDecay.meanR2(0, 1e3)) < 1e-12, "LD decay identical");
}

// =================================================================== Step 6
static void testLD() {
    std::cout << "[LD]\n";
    Rng r(21);
    const size_t nHap = 200, S = 400;
    BitRows rows;
    rows.nHap = nHap;
    rows.words = wordsFor(nHap);
    std::vector<std::vector<int>> g(S, std::vector<int>(nHap));
    // correlated rows: each row copies the previous with 10% flips
    for (size_t j = 0; j < S; ++j) {
        uint32_t c = 0;
        std::vector<uint64_t> row(rows.words, 0);
        for (size_t h = 0; h < nHap; ++h) {
            int v = j == 0 ? (r.uniform() < 0.4) : g[j - 1][h];
            if (j > 0 && r.uniform() < 0.1) v ^= 1;
            g[j][h] = v;
            if (v) { row[h / 64] |= 1ULL << (h % 64); ++c; }
        }
        rows.data.insert(rows.data.end(), row.begin(), row.end());
        rows.pos.push_back(static_cast<double>(j) * 250);
        rows.count.push_back(c);
        rows.globalIdx.push_back(static_cast<int64_t>(j));
    }
    LDOptions o;
    o.maxDist = 20000;
    o.binSize = 1000;
    o.threads = 4;
    LDDecay dec;
    computeLD(rows, "1", o, 0, 1e5, dec, nullptr, nullptr);
    // naive reference
    std::vector<double> sum(dec.sum.size(), 0);
    std::vector<uint64_t> cnt(dec.sum.size(), 0);
    for (size_t i = 0; i < S; ++i)
        for (size_t j = i + 1; j < S; ++j) {
            const double dd = rows.pos[j] - rows.pos[i];
            if (dd > o.maxDist) break;
            double pa = 0, pb = 0, pab = 0;
            for (size_t h = 0; h < nHap; ++h) { pa += g[i][h]; pb += g[j][h]; pab += g[i][h] * g[j][h]; }
            pa /= nHap; pb /= nHap; pab /= nHap;
            const double den = pa * (1 - pa) * pb * (1 - pb);
            if (den <= 0) continue;
            const double D = pab - pa * pb;
            sum[static_cast<size_t>(dd / o.binSize)] += D * D / den;
            ++cnt[static_cast<size_t>(dd / o.binSize)];
        }
    bool ok = true;
    for (size_t b = 0; b < sum.size(); ++b)
        ok = ok && cnt[b] == dec.n[b] && std::fabs(sum[b] - dec.sum[b]) < 1e-9 * (1 + sum[b]);
    CHECK(ok, "binned r2 identical to naive double loop (OpenMP, 4 threads)");
    // genotype (dosage) r2 vs naive Pearson correlation of g = h(2k) + h(2k+1)
    {
        LDOptions og = o;
        og.genotype = true;
        LDDecay dg;
        computeLD(rows, "1", og, 0, 1e5, dg, nullptr, nullptr);
        std::vector<double> gsum(dg.sum.size(), 0);
        const size_t nI = nHap / 2;
        for (size_t i = 0; i < S; ++i)
            for (size_t j = i + 1; j < S; ++j) {
                const double dd = rows.pos[j] - rows.pos[i];
                if (dd > o.maxDist) break;
                double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
                for (size_t k = 0; k < nI; ++k) {
                    const double x = g[i][2 * k] + g[i][2 * k + 1], y = g[j][2 * k] + g[j][2 * k + 1];
                    sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
                }
                const double vx = nI * sxx - sx * sx, vy = nI * syy - sy * sy;
                if (vx <= 0 || vy <= 0) continue;
                const double c = nI * sxy - sx * sy;
                gsum[static_cast<size_t>(dd / o.binSize)] += c * c / (vx * vy);
            }
        bool gok = true;
        for (size_t b = 0; b < gsum.size(); ++b) gok = gok && std::fabs(gsum[b] - dg.sum[b]) < 1e-9 * (1 + gsum[b]);
        CHECK(gok, "genotype r2 identical to naive Pearson correlation of dosages");
    }
    CHECK(dec.sum[0] / dec.n[0] > dec.sum[15] / dec.n[15], "r2 decays with distance in test data");
    LDOptions so;
    so.maf = 0.45;
    const auto sel = selectSnps(rows.count, nHap, so);
    bool sok = true;
    for (size_t i : sel) sok = sok && std::min(rows.count[i], uint32_t(nHap - rows.count[i])) >= 0.45 * nHap;
    CHECK(sok, "MAF filter");
}

// =================================================================== Step 7
static void testMaps() {
    std::cout << "[map estimation]\n";
    CHECK_NEAR(hillWeirR2(0.0, 1e12), 10.0 / 22.0, 1e-9, "Hill-Weir r2 at C=0, n=inf");
    CHECK(hillWeirR2(100, 100) < hillWeirR2(10, 100), "Hill-Weir decreasing in C");
    // Realised map from a full simulation with a 2-level map: rate 1 cM/Mb then 5 cM/Mb.
    SimParams p;
    p.model = Model::SMCPrime;
    p.nHap = 50;
    p.Ne = 10000;
    p.mu = 1.25e-8;
    p.seed = 5;
    p.outPrefix = (tmpDir() / "maptest").string();
    const fs::path mf = tmpDir() / "twolevel.map";
    {
        std::ofstream f(mf);
        f << "chr pos rate cM\n1 0 1 0\n1 5000000 5 5\n1 10000000 0 30\n";
    }
    p.mapFile = mf.string();
    p.estimateMap = true;
    p.ld.windowSize = 2.5e6;
    p.realisedWindow = 2.5e6;
    p.threads = 4;
    CoalescentSimulator sim(p);
    sim.run();
    std::ifstream rm(p.outPrefix + ".map.realised.tsv");
    std::string line;
    std::getline(rm, line);
    std::vector<double> rates;
    while (std::getline(rm, line)) {
        std::istringstream ss(line);
        std::string c;
        double a, b, ev, rate, cm;
        ss >> c >> a >> b >> ev >> rate >> cm;
        rates.push_back(rate);
    }
    CHECK(rates.size() == 4, "realised map windows");
    if (rates.size() == 4) {
        CHECK_NEAR(0.5 * (rates[0] + rates[1]), 1.0, 0.1, "realised rate, 1 cM/Mb region");
        CHECK_NEAR(0.5 * (rates[2] + rates[3]), 5.0, 0.4, "realised rate, 5 cM/Mb region");
    }
    std::ifstream lm(p.outPrefix + ".map.ld.tsv");
    std::getline(lm, line);
    std::vector<double> ldr;
    while (std::getline(lm, line)) {
        std::istringstream ss(line);
        std::string c;
        double a, b, pairs, rho, rate, cm;
        ss >> c >> a >> b >> pairs >> rho >> rate >> cm;
        ldr.push_back(rate);
    }
    if (ldr.size() >= 4) {
        const double lo = 0.5 * (ldr[0] + ldr[1]), hi = 0.5 * (ldr[2] + ldr[3]);
        std::cout << "    LD-based rates: low region " << lo << ", high region " << hi << " cM/Mb\n";
        CHECK(hi > 2.0 * lo, "LD-based map detects the 5x higher rate region");
    } else {
        CHECK(false, "LD-based map windows written");
    }
}

static void testVcf() {
    std::cout << "[VCF reader]\n";
    const fs::path f = tmpDir() / "t.vcf";
    {
        std::ofstream o(f);
        o << "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n"
             "chr22\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t1|1\t0|0\n"
             "chr22\t200\t.\tA\tGT\t.\tPASS\t.\tGT\t0|1\t1|1\t0|0\n"      // indel: skipped
             "chr22\t300\t.\tC\tT\t.\tPASS\tMULTI_ALLELIC\tGT\t0|1\t1|1\t0|0\n"  // skipped
             "chr22\t400\t.\tC\tT\t.\tPASS\t.\tGT:DS\t1|0:1\t0|0:0\t0|1:1\n";
    }
    const VcfData v = readPhasedVcf(f.string(), {"C", "A"}, 0.0);
    CHECK(v.chrom == "22" && v.snpsKept == 2 && v.rows.nHap == 4, "SNP filtering and sample subset");
    // sample order follows the file: A then C -> haplotypes A0 A1 C0 C1
    CHECK(v.rows.count[0] == 1 && (v.rows.row(0)[0] & 0xF) == 0b0010, "first SNP genotypes");
    CHECK(v.rows.count[1] == 2 && (v.rows.row(1)[0] & 0xF) == 0b1001, "second SNP genotypes");
}

int main(int argc, char** argv) {
    const std::string only = argc > 1 ? argv[1] : "";
    const auto run = [&](const std::string& name, const std::function<void()>& f) {
        if (only.empty() || name.find(only) != std::string::npos) f();
    };
    run("rng", testRng);
    run("demography", testDemography);
    run("map", testGeneticMap);
    run("tree", testInitialTree);
    run("tree", testLeafCache);
    run("smc", [] { testStationarity(Model::SMC); });
    run("smcp", [] { testStationarity(Model::SMCPrime); });
    run("macs", [] { testStationarity(Model::MaCS, 0); });
    run("macs", [] { testStationarity(Model::MaCS, 1e5); });
    run("macs", testMacsArg);
    run("cov", [] { testTwoLocusCovariance(Model::SMCPrime, 0, CovRef::SMCPrimeChain); });
    run("cov", [] { testTwoLocusCovariance(Model::SMC, 0, CovRef::SMCChain); });
    run("cov", [] { testTwoLocusCovariance(Model::MaCS, 0, CovRef::SMCPrimeChain); });
    run("cov", [] { testTwoLocusCovariance(Model::MaCS, 1e6, CovRef::ARG); });
    run("store", testTransposeAndStore);
    run("api", testLibraryApi);
    run("ld", testLD);
    run("mapest", testMaps);
    run("vcf", testVcf);
    std::cout << "\n" << g_pass << " checks passed, " << g_fail << " failed\n";
    return g_fail ? 1 : 0;
}
