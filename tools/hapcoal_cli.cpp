/**
 * @file hapcoal_cli.cpp
 * @brief `hapcoal` command-line tool, built on the public library API only.
 *
 *   hapcoal simulate [options]   simulate haplotypes and write files
 *   hapcoal ld [options]         observed LD decay of a phased VCF
 *   hapcoal chr22 [options]      validation: 1000 Genomes CEU chr22 vs SMC, SMC′, MaCS
 */
#include <hapcoal/hapcoal.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace hapcoal;
namespace fs = std::filesystem;

namespace {

// ============================================================== simulate ====
const char* kSimulateUsage = R"(hapcoal simulate [options]

Population genetics
  --model smc|smcp|macs     coalescent approximation (default macs)
  --n-ind N --ploidy P      N individuals of ploidy P (default ploidy 2)
  --n-hap N                 alternatively: N haplotypes
  --Ne N                    reference effective population size N0 (default 10000)
  --mu X                    mutation rate per bp per generation (default 1.25e-8)
  --rho X                   recombination rate per bp per generation if no map (default 1e-8)
  --demography FILE         piecewise-constant size history: lines "generations_ago N"
  --epoch G:N               add an epoch (repeatable), alternative to --demography
  --history BP              MaCS history window in bp (default 100000)

Genome
  --chr NAME:LENGTH         chromosome of given length (repeatable)
  --chr NAME:START-END      chromosome region (1-based, inclusive)
  --chromosomes FILE        lines "name length [rate]"
  --map FILE                genetic map (HapMap/PLINK/1000G formats; regions default to map range)
  --hotspots                without a map: simulate a hotspot map (mean rate = --rho)
  --hotspot-density X       hotspots per Mb (default 20)

Output
  --out PREFIX              output prefix (default out/sim)
  --text                    also write 0/1 text haplotypes (PREFIX.hap.txt)
  --no-ld                   skip LD decay
  --ld-max-dist BP --ld-bin BP --maf X --ld-max-snps N
  --ld-genotype             r2 of unphased diploid dosages instead of haplotypes
  --ld-matrix               sparse LD matrix (PREFIX.ld_matrix.tsv) up to --ld-matrix-dist BP
  --ld-min-r2 X             only write pairs with r2 >= X
  --ld-region CHR:A-B       dense r2 matrix for a region (PREFIX.ld_region.tsv)
  --estimate-map            infer a genetic map from the simulated haplotypes (PREFIX.map.ld.tsv)
  --ld-map-window BP        window for --estimate-map (default 1e6)
  --realised-window BP      window of the realised map (default 1e5)

Run
  --seed N --threads N --keep-temp --quiet
)";

/// Parse the options of `hapcoal simulate` into SimParams.
SimParams parseSimulateArgs(const std::vector<std::string>& a) {
    SimParams p;
    p.outPrefix = "out/sim";
    p.keepHaplotypes = false;  // the tool writes files; no need to hold matrices in memory
    p.log = &std::cerr;
    auto need = [&](size_t& i) -> const std::string& {
        if (i + 1 >= a.size()) throw std::invalid_argument("option " + a[i] + " needs a value");
        return a[++i];
    };
    auto parseChr = [&](const std::string& s) {
        ChromosomeSpec c;
        const auto colon = s.find(':');
        if (colon == std::string::npos) { c.name = GeneticMap::normaliseChr(s); return c; }
        c.name = GeneticMap::normaliseChr(s.substr(0, colon));
        const std::string r = s.substr(colon + 1);
        const auto dash = r.find('-');
        if (dash == std::string::npos) { c.start = 0; c.end = std::stod(r); }
        else { c.start = std::stod(r.substr(0, dash)) - 1; c.end = std::stod(r.substr(dash + 1)); }
        return c;
    };
    for (size_t i = 0; i < a.size(); ++i) {
        const std::string& o = a[i];
        if (o == "--model") p.model = parseModel(need(i));
        else if (o == "--n-ind") p.nInd = std::stoi(need(i));
        else if (o == "--ploidy") p.ploidy = std::stoi(need(i));
        else if (o == "--n-hap") p.nHap = std::stoi(need(i));
        else if (o == "--Ne" || o == "--ne") p.Ne = std::stod(need(i));
        else if (o == "--mu") p.mu = std::stod(need(i));
        else if (o == "--rho" || o == "--r") p.rho = std::stod(need(i));
        else if (o == "--demography") p.demographyFile = need(i);
        else if (o == "--epoch") {
            const std::string& v = need(i);
            const auto c = v.find(':');
            if (c == std::string::npos) throw std::invalid_argument("--epoch expects G:N");
            p.epochs.emplace_back(std::stod(v.substr(0, c)), std::stod(v.substr(c + 1)));
        } else if (o == "--history") p.history = std::stod(need(i));
        else if (o == "--chr") p.chromosomes.push_back(parseChr(need(i)));
        else if (o == "--chromosomes") {
            std::ifstream in(need(i));
            if (!in) throw std::runtime_error("Cannot open " + a[i]);
            std::string line;
            while (std::getline(in, line)) {
                if (line.empty() || line[0] == '#') continue;
                std::istringstream ss(line);
                ChromosomeSpec c;
                double len;
                if (!(ss >> c.name >> len)) continue;
                c.name = GeneticMap::normaliseChr(c.name);
                c.end = len;
                double r;
                if (ss >> r) c.rate = r;
                p.chromosomes.push_back(c);
            }
        } else if (o == "--map") p.mapFile = need(i);
        else if (o == "--hotspots") p.hotspots = true;
        else if (o == "--hotspot-density") { p.hotspots = true; p.hotspotParams.densityPerMb = std::stod(need(i)); }
        else if (o == "--out") p.outPrefix = need(i);
        else if (o == "--text") p.writeText = true;
        else if (o == "--no-ld") p.ldDecay = false;
        else if (o == "--ld-max-dist") p.ld.maxDist = std::stod(need(i));
        else if (o == "--ld-bin") p.ld.binSize = std::stod(need(i));
        else if (o == "--maf") p.ld.maf = std::stod(need(i));
        else if (o == "--ld-max-snps") p.ld.maxSnps = static_cast<size_t>(std::stoll(need(i)));
        else if (o == "--ld-genotype") p.ld.genotype = true;
        else if (o == "--ld-matrix") p.ldMatrix = true;
        else if (o == "--ld-matrix-dist") { p.ldMatrix = true; p.ld.matrixMaxDist = std::stod(need(i)); }
        else if (o == "--ld-min-r2") p.ld.matrixMinR2 = std::stod(need(i));
        else if (o == "--ld-region") p.ldRegion = need(i);
        else if (o == "--estimate-map") p.estimateMap = true;
        else if (o == "--ld-map-window") p.ld.windowSize = std::stod(need(i));
        else if (o == "--realised-window") p.realisedWindow = std::stod(need(i));
        else if (o == "--seed") p.seed = static_cast<uint64_t>(std::stoull(need(i)));
        else if (o == "--threads") p.threads = std::stoi(need(i));
        else if (o == "--keep-temp") p.keepTemp = true;
        else if (o == "--quiet") p.log = nullptr;
        else throw std::invalid_argument("unknown option " + o + "\n\n" + kSimulateUsage);
    }
    return p;
}

// ==================================================================== ld ====
struct LdArgs {
    std::string vcf, panel, pop, out = "out/obs";
    double regionStart = 0, regionEnd = 1e18;
    LDOptions ld;
};

LdArgs parseLdArgs(const std::vector<std::string>& a) {
    LdArgs r;
    auto need = [&](size_t& i) -> const std::string& {
        if (i + 1 >= a.size()) throw std::invalid_argument("option " + a[i] + " needs a value");
        return a[++i];
    };
    for (size_t i = 0; i < a.size(); ++i) {
        const std::string& o = a[i];
        if (o == "--vcf") r.vcf = need(i);
        else if (o == "--panel") r.panel = need(i);
        else if (o == "--pop") r.pop = need(i);
        else if (o == "--out") r.out = need(i);
        else if (o == "--maf") r.ld.maf = std::stod(need(i));
        else if (o == "--ld-max-dist") r.ld.maxDist = std::stod(need(i));
        else if (o == "--ld-bin") r.ld.binSize = std::stod(need(i));
        else if (o == "--ld-max-snps") r.ld.maxSnps = static_cast<size_t>(std::stoll(need(i)));
        else if (o == "--ld-genotype") r.ld.genotype = true;
        else if (o == "--region") {
            const std::string& v = need(i);
            const auto d = v.find('-');
            r.regionStart = std::stod(v.substr(0, d));
            r.regionEnd = std::stod(v.substr(d + 1));
        } else if (o == "--threads") r.ld.threads = std::stoi(need(i));
        else throw std::invalid_argument("unknown option " + o);
    }
    if (r.vcf.empty()) throw std::invalid_argument("ld: --vcf is required");
    if (!r.panel.empty() && r.pop.empty()) throw std::invalid_argument("ld: --panel needs --pop");
    return r;
}

/// Observed LD decay of a VCF; writes OUT.ld_decay.tsv and reports to stderr.
ObservedLD runObservedLD(const LdArgs& a) {
    std::vector<std::string> samples;
    if (!a.panel.empty()) samples = readPanelSamples(a.panel, a.pop);
    std::cerr << "reading " << a.vcf
              << (samples.empty() ? "" : " (" + std::to_string(samples.size()) + " samples of " + a.pop + ")") << " ...\n";
    ObservedLD r = observedLD(a.vcf, samples, a.regionStart, a.regionEnd, a.ld);
    std::cerr << "  " << r.stats.linesRead << " records, " << r.stats.snpsKept << " biallelic SNPs with MAF >= "
              << a.ld.maf << ", " << r.nHap << " haplotypes";
    if (r.stats.unphased) std::cerr << ", WARNING: " << r.stats.unphased << " unphased genotypes";
    std::cerr << "\n";
    const fs::path pre(a.out);
    if (pre.has_parent_path()) fs::create_directories(pre.parent_path());
    std::ofstream out(a.out + ".ld_decay.tsv");
    out << "chrom\tbin_start_bp\tbin_mid_bp\tmean_r2\tn_pairs\n";
    writeDecay(out, r.chrom, r.decay);
    return r;
}

// ================================================================= chr22 ====
// Validation: observed CEU LD decay vs simulations under the Tennessen et al.
// (2012) demography and under a bottleneck + growth model fitted to the data.
int runChr22(const std::vector<std::string>& args) {
    std::string data = "data", out = "out/chr22", models = "smc,smcp,macs", pop = "CEU", demoFile;
    int threads = 0, nInd = 0;
    double history = 1e5, Ne = 10000, mu = 2.36e-8;
    uint64_t seed = 2026;
    double maxDist = 500000;
    bool genotype = false, constant = false;
    int fitEvals = 120;  // max simulations for the demography fit
    for (size_t i = 0; i < args.size(); ++i) {
        auto need = [&]() -> const std::string& {
            if (i + 1 >= args.size()) throw std::invalid_argument("option needs a value");
            return args[++i];
        };
        if (args[i] == "--data") data = need();
        else if (args[i] == "--out") out = need();
        else if (args[i] == "--models") models = need();
        else if (args[i] == "--pop") pop = need();
        else if (args[i] == "--threads") threads = std::stoi(need());
        else if (args[i] == "--n-ind") nInd = std::stoi(need());
        else if (args[i] == "--history") history = std::stod(need());
        else if (args[i] == "--seed") seed = std::stoull(need());
        else if (args[i] == "--ld-max-dist") maxDist = std::stod(need());
        else if (args[i] == "--ld-genotype") genotype = true;
        else if (args[i] == "--demography") demoFile = need();
        else if (args[i] == "--constant") constant = true;
        else if (args[i] == "--Ne") Ne = std::stod(need());
        else if (args[i] == "--mu") mu = std::stod(need());
        else if (args[i] == "--fit-evals") fitEvals = std::stoi(need());
        else throw std::invalid_argument("unknown option " + args[i]);
    }
    const std::string mapFile = data + "/genetic_map_GRCh37_chr22.txt";
    const std::string vcf = data + "/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz";
    const std::string panel = data + "/integrated_call_samples_v3.20130502.ALL.panel";
    if (!fs::exists(mapFile)) throw std::runtime_error("missing " + mapFile + " (run: make data)");
    fs::create_directories(out);

    const auto maps = GeneticMap::loadFile(mapFile, "22");
    const GeneticMap& m22 = maps.at("22");
    std::cerr << "chr22 map: " << std::fixed << std::setprecision(0) << m22.start() << "-" << m22.end() << " bp, "
              << std::setprecision(2) << m22.totalMorgans() * 100 << " cM\n" << std::defaultfloat << std::setprecision(6);

    // ---- observed LD
    LDDecay obs;
    bool haveObs = false;
    size_t nObsHap = 0;
    if (fs::exists(vcf) && fs::exists(panel)) {
        LdArgs la;
        la.vcf = vcf;
        la.panel = panel;
        la.pop = pop;
        la.out = out + "/observed_" + pop;
        la.regionStart = m22.start() + 1;
        la.regionEnd = m22.end();
        la.ld.maxDist = maxDist;
        la.ld.genotype = genotype;
        la.ld.threads = threads;
        const ObservedLD o = runObservedLD(la);
        obs = o.decay;
        nObsHap = o.nHap;
        haveObs = true;
    } else {
        std::cerr << "warning: 1000 Genomes VCF/panel not found in " << data << " - observed LD skipped\n";
    }
    if (nInd <= 0) nInd = haveObs ? static_cast<int>(nObsHap / 2) : 99;

    // ---- common parameters
    SimParams base;
    base.nInd = nInd;
    base.ploidy = 2;
    base.Ne = Ne;
    base.mu = mu;
    base.ld.genotype = genotype;
    base.ld.maxDist = maxDist;
    base.mapFile = mapFile;
    base.chromosomes = {{"22", 0, 0, -1}};
    base.history = history;
    base.seed = seed;
    base.threads = threads;
    base.keepHaplotypes = false;

    // ---- demographies to compare
    struct Demo {
        std::string label;
        std::vector<std::pair<double, double>> epochs;
        std::string file;
    };
    std::vector<Demo> demos;
    if (!demoFile.empty()) demos.push_back({"user", {}, demoFile});
    else if (constant) demos.push_back({"constant", {}, ""});
    else {
        demos.push_back({"tennessen2012", Demography::tennessen2012CEU(), ""});
        if (haveObs && fitEvals > 0) {
            std::cerr << "fitting a bottleneck + growth demography to the observed LD decay (SMC', Nelder-Mead) ...\n";
            SimParams fp = base;
            fp.model = Model::SMCPrime;
            const DemographyFitResult fit = fitBottleneckGrowth(obs, fp, fitEvals, &std::cerr);
            std::ofstream df(out + "/fitted_demography.tsv");
            df << "# LD-fitted piecewise-constant demography for " << pop << " chr22 (generations_ago  N)\n";
            for (const auto& e : fit.epochs)
                df << std::fixed << std::setprecision(1) << e.first << '\t' << std::setprecision(0) << e.second << '\n';
            std::cerr << "  best fit RMSE " << fit.rmse << " -> " << out << "/fitted_demography.tsv\n";
            demos.push_back({"ld_fit", fit.epochs, ""});
        }
    }

    // ---- simulations and comparison with the observed decay
    std::ofstream cmp(out + "/comparison.tsv");
    cmp << "model\tdemography\thaplotypes\tsnps\trecombinations\tinvisible\tnon_ancestral\tsim_seconds\trmse_r2\tmean_diff_r2\n";
    for (const auto& dm : demos) {
        std::stringstream ms(models);
        for (std::string mname; std::getline(ms, mname, ',');) {
            SimParams p = base;
            p.model = parseModel(mname);
            p.epochs = dm.epochs;
            p.demographyFile = dm.file;
            p.outPrefix = out + "/sim_" + mname + "_" + dm.label;
            p.estimateMap = true;
            p.log = &std::cerr;
            CoalescentSimulator sim(p);
            const SimulationResult res = sim.run();
            const ChromStats& st = res.chromosomes.front().stats;
            double md = 0;
            const double e = haveObs ? decayRmse(res.ldDecay, obs, maxDist, &md) : NAN;
            cmp << modelName(p.model) << '\t' << dm.label << '\t' << p.numHaplotypes() << '\t' << st.snps << '\t'
                << st.recombinations << '\t' << st.invisible << '\t' << st.nonLocal << '\t' << st.seconds << '\t' << e
                << '\t' << md << '\n';
            std::cerr << modelName(p.model) << " [" << dm.label << "]: RMSE(r2) vs observed = " << e
                      << " (mean sim-obs " << md << ")\n";
        }
    }
    std::cerr << "comparison written to " << out << "/comparison.tsv\n";
    return 0;
}

void help() {
    std::cout << "hapcoal " HAPCOAL_VERSION_STRING " — SMC / SMC' / MaCS coalescent haplotype simulator\n\n"
                 "  hapcoal simulate [options]   simulate haplotypes\n"
                 "  hapcoal ld --vcf F [--panel P --pop CEU] [--region A-B] [--maf X] [--out PREFIX]\n"
                 "                              observed LD decay of a phased VCF\n"
                 "  hapcoal chr22 [--data DIR] [--out DIR] [--models smc,smcp,macs] [--history BP]\n"
                 "               [--fit-evals N] [--ld-genotype] [--demography FILE | --constant]\n"
                 "                              human chr22 / 1000 Genomes CEU validation test case\n\n"
              << kSimulateUsage;
}

}  // namespace

int main(int argc, char** argv) {
    std::vector<std::string> args(argv + 1, argv + argc);
    if (args.empty() || args[0] == "-h" || args[0] == "--help") {
        help();
        return args.empty() ? 1 : 0;
    }
    try {
        const std::string cmd = args[0];
        args.erase(args.begin());
        for (const auto& a : args)
            if (a == "-h" || a == "--help") {
                help();
                return 0;
            }
        if (cmd == "simulate") {
            CoalescentSimulator sim(parseSimulateArgs(args));
            sim.run();
        } else if (cmd == "ld") {
            runObservedLD(parseLdArgs(args));
        } else if (cmd == "chr22") {
            return runChr22(args);
        } else {
            help();
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 2;
    }
    return 0;
}
