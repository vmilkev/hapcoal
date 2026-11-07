/**
 * @file CoalescentSimulator.cpp
 * @brief SMC / SMC′ / MaCS simulation of one chromosome and the run() pipeline.
 *
 * Pipeline of run():
 *  1. build the genetic map of every chromosome;
 *  2. simulate all chromosomes in parallel (one independent RNG stream each),
 *     appending SNP rows to a SnpBuffer (memory or temporary file);
 *  3. if an output prefix is set: write the haplotype file, index files and maps;
 *  4. per chromosome, in input order: LD decay / matrices / LD-based map,
 *     haplotype-major matrix, then hand the ChromosomeResult to the caller.
 */
#include "hapcoal/CoalescentSimulator.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "HaplotypeStore.hpp"
#include "LocalTree.hpp"
#include "MacsHistory.hpp"
#include "Rng.hpp"

namespace hapcoal {

namespace fs = std::filesystem;

std::string modelName(Model m) {
    switch (m) {
        case Model::SMC: return "SMC";
        case Model::SMCPrime: return "SMC'";
        case Model::MaCS: return "MaCS";
    }
    return "?";
}

Model parseModel(const std::string& s0) {
    std::string s;
    for (char c : s0) s.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    if (s == "smc") return Model::SMC;
    if (s == "smcp" || s == "smc'" || s == "smcprime" || s == "smc_prime") return Model::SMCPrime;
    if (s == "macs") return Model::MaCS;
    throw std::invalid_argument("Unknown model '" + s0 + "' (use smc, smcp or macs)");
}

// ======================================================== SimulationResult ====
size_t SimulationResult::numSnps() const {
    size_t n = 0;
    for (const auto& c : chromosomes) n += c.positions.size();
    return n;
}

size_t SimulationResult::rowWordOffset(size_t c) const {
    size_t off = 0;
    for (size_t i = 0; i < c && i < chromosomes.size(); ++i) off += wordsFor(chromosomes[i].positions.size());
    return off;
}

std::vector<uint64_t> SimulationResult::genomeRow(size_t h) const {
    std::vector<uint64_t> row(rowWordOffset(chromosomes.size()), 0);
    size_t off = 0;
    for (const auto& c : chromosomes) {
        if (c.haplotypes.empty() && !c.positions.empty())
            throw std::logic_error("genomeRow: haplotypes were not kept (keepHaplotypes=false or callback used)");
        std::copy(c.haplotypes.row(h), c.haplotypes.row(h) + c.haplotypes.words, row.begin() + static_cast<long>(off));
        off += wordsFor(c.positions.size());
    }
    return row;
}

// ==================================================================== setup ====
CoalescentSimulator::CoalescentSimulator(SimParams p) : p_(std::move(p)) {
    if (p_.ploidy < 1) throw std::invalid_argument("ploidy must be >= 1");
    if (p_.numHaplotypes() < 2) throw std::invalid_argument("need at least 2 haplotypes (nInd or nHap)");
    if (p_.Ne <= 0 || p_.mu < 0 || p_.rho < 0) throw std::invalid_argument("Ne must be > 0, mu and rho >= 0");
    if (p_.history < 0) throw std::invalid_argument("history must be >= 0");
    if (!p_.demographyFile.empty())
        demo_ = Demography::fromFile(p_.demographyFile, p_.Ne, p_.ploidy);
    else if (!p_.epochs.empty())
        demo_ = Demography::fromGenerations(p_.epochs, p_.Ne, p_.ploidy);
    if (!p_.mapFile.empty()) {
        // a 3-column map file has no chromosome column: assign it to the only chromosome given
        const std::string def = p_.chromosomes.size() == 1 ? p_.chromosomes[0].name : "1";
        userMaps_ = GeneticMap::loadFile(p_.mapFile, def);
    }
    if (p_.chromosomes.empty()) {
        if (userMaps_.empty()) throw std::invalid_argument("no chromosomes given (chromosomes or mapFile)");
        for (const auto& kv : userMaps_) p_.chromosomes.push_back({kv.first, 0.0, 0.0, -1.0});
    }
}

GeneticMap CoalescentSimulator::buildMap(ChromosomeSpec& spec, uint64_t stream) const {
    auto it = userMaps_.find(GeneticMap::normaliseChr(spec.name));
    if (it != userMaps_.end()) {
        const GeneticMap& m = it->second;
        if (spec.end <= 0) {  // whole map range
            spec.start = m.start();
            spec.end = m.end();
            return m;
        }
        if ((spec.start < m.start() || spec.end > m.end()) && p_.log)
            *p_.log << "warning: chromosome " << spec.name << " region exceeds the genetic map; clipped to ["
                    << m.start() << ", " << m.end() << "]\n";
        spec.start = std::max(spec.start, m.start());
        spec.end = std::min(spec.end, m.end());
        return m.clipped(spec.start, spec.end);
    }
    if (!(spec.end > spec.start))
        throw std::invalid_argument("chromosome " + spec.name + ": no length given and no genetic map");
    const double rate = spec.rate >= 0 ? spec.rate : p_.rho;
    if (p_.hotspots)
        return GeneticMap::withHotspots(spec.start, spec.end, rate, p_.hotspotParams, p_.seed ^ 0xA5A5A5A5ULL, stream);
    return GeneticMap::uniform(spec.start, spec.end, rate);
}

// ============================================================ one chromosome ====
//
// The local tree at position x is a Kingman-distributed genealogy. Moving
// along the chromosome, recombination events occur at rate r(x)·p·N0·T(x) per bp
// (T = total branch length in coalescent units), i.e. at rate p·N0·T per Morgan;
// the distance to the next event is therefore exponential on the genetic scale.
// Between events, mutations fall on the current tree as a Poisson process with
// rate mu·p·N0·T per bp, each on a branch chosen proportionally to its length.
// At an event a point is chosen uniformly on the tree (branch b, time tr); the
// lineage above it "floats" and coalesces back according to the model.
ChromStats CoalescentSimulator::simulateChromosome(const ChromosomeSpec& spec, const GeneticMap& map,
                                                   uint64_t stream, SnpBuffer& snps, RealisedMap* realised,
                                                   TreeSpanCallback* cb) const {
    const auto t0 = std::chrono::steady_clock::now();
    Rng rng(p_.seed, stream);
    const int n = p_.numHaplotypes();
    const double scale = static_cast<double>(p_.ploidy) * p_.Ne;  // generations per coalescent unit
    const Model model = p_.model;

    LocalTree tree(n);
    tree.enableLeafCache(512ull << 20);  // per-node leaf bitsets if they fit in 512 MB
    tree.buildInitial(demo_, rng);
    std::unique_ptr<MacsHistory> arg;
    if (model == Model::MaCS) {
        arg = std::make_unique<MacsHistory>();
        arg->init(tree);
    }

    ChromStats st;
    st.name = spec.name;
    st.length = spec.end - spec.start;
    st.lengthCM = (map.toMorgans(spec.end) - map.toMorgans(spec.start)) * 100.0;

    double x = spec.start;             // current position (bp)
    double g = map.toMorgans(x);       // current genetic position (M)
    const double gEnd = map.toMorgans(spec.end);
    std::vector<double> mpos;

    for (;;) {
        const double T = tree.totalLength();
        // --- distance to the next recombination, drawn on the genetic scale
        const double gNext = g + rng.exp1() / (T * scale);
        const bool last = gNext >= gEnd;
        const double xNext = last ? spec.end : std::max(x, map.toBp(gNext));

        // --- mutations on the current tree over [x, xNext) (infinite sites)
        const double span = xNext - x;
        if (span > 0) {
            const int64_t m = rng.poisson(p_.mu * scale * T * span);
            if (m > 0) {
                mpos.resize(static_cast<size_t>(m));
                for (auto& q : mpos) q = x + rng.uniform() * span;
                std::sort(mpos.begin(), mpos.end());
                for (double q : mpos) {
                    const int v = tree.branchAt(rng.uniform() * T);  // branch ∝ length
                    uint64_t* row = snps.beginRow();
                    snps.commitRow(q, tree.leafBits(v, row));       // carriers = leaves below v
                }
            }
            if (realised) realised->addSpan(x, xNext, T);
            if (cb) cb->span(x, xNext, tree);
        }
        if (last) break;
        x = xNext;
        g = gNext;
        if (realised) realised->addEvent(x);
        ++st.recombinations;

        // --- recombination point on the local tree: branch proportional to length, time uniform
        const int b = tree.branchAt(rng.uniform() * T);
        const int op = tree.parent(b);
        const int s = tree.sibling(b);
        const double tr = rng.uniform(tree.time(b), tree.time(op));

        if (model == Model::MaCS) {
            // MaCS: drop ARG edges unused for more than `history` bp, then let the
            // floating lineage coalesce with any retained edge (MacsHistory).
            arg->prune(x - p_.history);
            const MacsHistory::Result r = arg->recombine(tree, b, op, s, tr, x, demo_, rng);
            if (r.nonLocal) ++st.nonLocal;
            st.maxArgEdges = std::max(st.maxArgEdges, arg->numEdges());
            if (r.u == s && r.t == tree.time(op)) {
                ++st.invisible;  // lineage re-joined its old position: local tree unchanged
            } else {
                tree.regraft(b, r.u, r.t);
            }
        } else {
            // SMC / SMC': coalescence of the floating lineage with the lineages of
            // the tree. k(t) = K(t) for SMC' (K counts the original tree, including
            // the lineage's own old branch above tr); for SMC the own branch is
            // excluded, i.e. k(t) = K(t) - 1 below t_op and K(t) above it, which is
            // obtained by starting at K-1 and skipping op's coalescence event.
            // The hazard k(t)/lambda(t) is piecewise constant between node times and
            // epoch boundaries; the coalescence time is found by inverting its integral.
            const bool smc = model == Model::SMC;
            const TimeIndex& TI = tree.times();
            double t = tr;
            int k = tree.lineagesAt(tr) - (smc ? 1 : 0);
            TimeIndex::Cursor c = TI.upperBound(tr);
            double E = rng.exp1();
            double tc;
            for (;;) {
                const double tn = TI.valid(c) ? TI.get(c).t : Demography::INF;  // next node time
                const double te = demo_.nextBoundary(t);                         // next epoch boundary
                const double tb = std::min(tn, te);
                const double rate = k / demo_.lambda(t);
                if (E < rate * (tb - t)) {
                    tc = t + E / rate;
                    break;
                }
                E -= rate * (tb - t);
                t = tb;
                if (tb == tn) {
                    const auto& item = TI.get(c);
                    if (!(smc && item.id == op)) k += item.delta;
                    TI.advance(c);
                }
            }
            // target lineage: uniform among the k lineages crossing tc
            const int exclude = (smc && tc < tree.time(op)) ? b : -1;
            int v = tree.pickBranchCrossing(tc, exclude, rng);
            if (v == b) {          // SMC': back onto its own branch, tree unchanged
                ++st.invisible;
                continue;
            }
            if (v == op) v = s;    // above op the remaining lineage is s's
            tree.regraft(b, v, tc);
        }
#ifdef HAPCOAL_DEBUG
        std::string err;
        if (!tree.validate(&err)) throw std::logic_error("invalid local tree: " + err);
        if (arg && !arg->validate(tree, &err)) throw std::logic_error("invalid ARG: " + err);
#endif
    }
    st.snps = snps.numSnps();
    st.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return st;
}

// ====================================================================== run ====
SimulationResult CoalescentSimulator::run(const ChromosomeCallback& onChromosome) {
    const auto tStart = std::chrono::steady_clock::now();
    const size_t nHap = static_cast<size_t>(p_.numHaplotypes());
    int threads = p_.threads;
#ifdef _OPENMP
    if (threads <= 0) threads = omp_get_max_threads();
#else
    threads = 1;
#endif
    p_.ld.threads = threads;

    // ---- where SNP rows live during the run
    const bool files = !p_.outPrefix.empty();
    fs::path tmpDir;
    bool ownTmpDir = false;
    if (!p_.tempDir.empty()) {
        tmpDir = p_.tempDir;
    } else if (files) {
        tmpDir = fs::path(p_.outPrefix + "_tmp");
        ownTmpDir = true;
    }
    const bool onDisk = !tmpDir.empty();
    if (files) {
        const fs::path prefix(p_.outPrefix);
        if (prefix.has_parent_path()) fs::create_directories(prefix.parent_path());
    }
    if (onDisk) fs::create_directories(tmpDir);

    std::ostream* log = p_.log;
    if (log) {
        *log << std::setprecision(6);
        *log << "hapcoal: model " << modelName(p_.model) << ", " << nHap << " haplotypes";
        if (p_.nInd > 0) *log << " (" << p_.nInd << " individuals x ploidy " << p_.ploidy << ")";
        *log << ", N0=" << p_.Ne << ", mu=" << p_.mu << ", threads=" << threads;
        if (p_.model == Model::MaCS) *log << ", history=" << p_.history << " bp";
        *log << "\n";
        if (!demo_.isConstant()) *log << "demography (piecewise constant):\n" << demo_.describe(p_.Ne, p_.ploidy);
    }

    // ---- 1. maps and per-chromosome storage
    std::vector<ChromData> chroms(p_.chromosomes.size());
    std::vector<RealisedMap> realised(p_.chromosomes.size());
    for (size_t i = 0; i < p_.chromosomes.size(); ++i) {
        auto& spec = p_.chromosomes[i];
        chroms[i].map = buildMap(spec, i + 1);
        chroms[i].name = spec.name;
        chroms[i].start = spec.start;
        chroms[i].end = spec.end;
        chroms[i].snps = std::make_unique<SnpBuffer>(
            nHap, onDisk ? (tmpDir / ("chr" + spec.name + ".snp")).string() : std::string());
        if (p_.keepTemp && onDisk) chroms[i].snps->keepFile();
        realised[i].init(spec.start, spec.end, p_.realisedWindow);
    }
    stats_.assign(chroms.size(), ChromStats{});

    // ---- 2. simulate (chromosomes in parallel, longest first)
    std::vector<size_t> order(chroms.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(),
              [&](size_t a, size_t b) { return chroms[a].end - chroms[a].start > chroms[b].end - chroms[b].start; });
    std::string error;
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (long long oi = 0; oi < static_cast<long long>(order.size()); ++oi) {
        const size_t i = order[static_cast<size_t>(oi)];
        try {
            SnpBuffer& w = *chroms[i].snps;
            stats_[i] = simulateChromosome(p_.chromosomes[i], chroms[i].map, i + 1, w, &realised[i]);
            w.close();
            chroms[i].pos = std::move(w.positions());
            chroms[i].count = std::move(w.counts());
            if (log) {
#pragma omp critical(log)
                *log << "  chr " << stats_[i].name << ": " << std::fixed << std::setprecision(1)
                     << stats_[i].length / 1e6 << " Mb, " << std::setprecision(2) << stats_[i].lengthCM << " cM, "
                     << stats_[i].recombinations << " recombinations (" << stats_[i].invisible << " invisible, "
                     << stats_[i].nonLocal << " non-ancestral), " << stats_[i].snps << " SNPs, "
                     << std::setprecision(2) << stats_[i].seconds << " s\n"
                     << std::defaultfloat << std::setprecision(6);
            }
        } catch (const std::exception& e) {
#pragma omp critical(err)
            error = "chromosome " + chroms[i].name + ": " + e.what();
        }
    }
    if (!error.empty()) throw std::runtime_error(error);

    auto phase = [&](const char* name, std::chrono::steady_clock::time_point& t0) {
        const auto t1 = std::chrono::steady_clock::now();
        if (log)
            *log << "  [" << name << ": " << std::fixed << std::setprecision(1)
                 << std::chrono::duration<double>(t1 - t0).count() << " s]\n" << std::defaultfloat << std::setprecision(6);
        t0 = t1;
    };
    auto tPhase = tStart;
    phase("simulation", tPhase);

    // Global SNP indices (as in the files; also used by the sparse LD matrix).
    {
        size_t total = 0;
        for (auto& c : chroms) {
            c.firstSnp = total;
            total += c.pos.size();
        }
    }

    // ---- 3. files: haplotypes, index files, true and realised maps
    if (files) {
        assembleHaplotypes(p_.outPrefix, chroms, nHap, p_.ploidy, p_.writeText, threads);
        phase("haplotype output", tPhase);
        std::ofstream tm(p_.outPrefix + ".map.true.tsv");
        tm << "chrom\tpos_bp\trate_cM_per_Mb\tpos_cM\n";
        for (const auto& c : chroms) writeTrueMap(tm, c.name, c.map);
        std::ofstream rm(p_.outPrefix + ".map.realised.tsv");
        rm << "chrom\tstart_bp\tend_bp\tevents\trate_cM_per_Mb\tpos_cM\n";
        for (size_t i = 0; i < chroms.size(); ++i)
            writeRealisedMap(rm, chroms[i].name, realised[i], static_cast<double>(p_.ploidy) * p_.Ne);
    }

    // ---- 4. per chromosome: LD, LD-based map, haplotype matrix, delivery
    const bool doLD = p_.ldDecay || p_.ldMatrix || p_.estimateMap || !p_.ldRegion.empty();
    std::unique_ptr<std::ofstream> mat, ldm, dec;
    std::string regChr;
    double regA = 0, regB = 0;
    if (doLD && files) {
        if (p_.ldMatrix) {
            mat = std::make_unique<std::ofstream>(p_.outPrefix + ".ld_matrix.tsv");
            *mat << "chrom\tsnp_i\tsnp_j\tpos_i\tpos_j\tr2\n";
        }
        if (p_.estimateMap) {
            ldm = std::make_unique<std::ofstream>(p_.outPrefix + ".map.ld.tsv");
            *ldm << "chrom\tstart_bp\tend_bp\tpairs\trho_per_kb\trate_cM_per_Mb\tpos_cM\n";
        }
        dec = std::make_unique<std::ofstream>(p_.outPrefix + ".ld_decay.tsv");
        *dec << "chrom\tbin_start_bp\tbin_mid_bp\tmean_r2\tn_pairs\n";
        if (!p_.ldRegion.empty()) {
            const auto colon = p_.ldRegion.find(':'), dash = p_.ldRegion.find('-');
            if (colon == std::string::npos || dash == std::string::npos)
                throw std::invalid_argument("ldRegion must be chr:start-end");
            regChr = GeneticMap::normaliseChr(p_.ldRegion.substr(0, colon));
            regA = std::stod(p_.ldRegion.substr(colon + 1, dash - colon - 1)) - 1;
            regB = std::stod(p_.ldRegion.substr(dash + 1)) - 1;
        }
    }
    decay_ = LDDecay{};
    SimulationResult result;
    result.nHap = nHap;
    result.ploidy = p_.ploidy;
    for (size_t i = 0; i < chroms.size(); ++i) {
        ChromData& c = chroms[i];
        ChromosomeResult cr;
        if (doLD) {
            LDOptions o = p_.ld;
            o.matrix = p_.ldMatrix && files;
            o.windows = p_.estimateMap;
            o.seed = p_.seed + std::hash<std::string>{}(c.name);  // thinning seed per chromosome
            const auto sel = selectSnps(c.count, nHap, o);
            const BitRows rows = loadRows(*c.snps, nHap, sel, c.pos, c.count, c.firstSnp);
            LDWindows win;
            computeLD(rows, c.name, o, c.start, c.end, cr.ldDecay, p_.estimateMap ? &win : nullptr, mat.get());
            decay_.add(cr.ldDecay);
            if (dec) writeDecay(*dec, c.name, cr.ldDecay);
            if (p_.estimateMap) {
                cr.ldMap = estimateLDMap(win, c.end, nHap, p_.Ne, p_.ploidy);
                if (ldm) writeLDMap(*ldm, c.name, cr.ldMap);
            }
            if (!regChr.empty() && GeneticMap::normaliseChr(c.name) == regChr)
                writeDenseMatrix(rows, regA, regB, p_.outPrefix + ".ld_region.tsv");
            if (log) *log << "  LD chr " << c.name << ": " << rows.size() << " SNPs with MAF >= " << o.maf << "\n";
        }
        cr.name = c.name;
        cr.start = c.start;
        cr.end = c.end;
        cr.positions = std::move(c.pos);
        cr.altCounts = std::move(c.count);
        cr.map = c.map;
        cr.realisedMap = realised[i];
        cr.stats = stats_[i];
        if (p_.keepHaplotypes) cr.haplotypes = toHaplotypeMatrix(*c.snps, nHap, threads);
        if (!p_.keepTemp) c.snps->release();  // free memory / delete the temporary file now
        if (onChromosome) {
            onChromosome(cr);
            cr.haplotypes = HaplotypeMatrix{};  // delivered: do not keep (bounded memory)
        }
        result.chromosomes.push_back(std::move(cr));
    }
    if (dec && chroms.size() > 1) writeDecay(*dec, "all", decay_);
    result.ldDecay = decay_;
    phase("maps + LD", tPhase);

    // ---- run log and cleanup
    if (files) {
        std::ofstream lg(p_.outPrefix + ".log");
        lg << "model\t" << modelName(p_.model) << "\nhaplotypes\t" << nHap << "\nploidy\t" << p_.ploidy
           << "\nN0\t" << p_.Ne << "\nmu\t" << p_.mu << "\nseed\t" << p_.seed << "\n";
        if (p_.model == Model::MaCS) lg << "history_bp\t" << p_.history << "\n";
        lg << "demography\n" << demo_.describe(p_.Ne, p_.ploidy);
        lg << "chrom\tlength_bp\tlength_cM\trecombinations\tinvisible\tnon_ancestral\tsnps\tmax_arg_edges\tseconds\n";
        for (const auto& s : stats_)
            lg << s.name << '\t' << s.length << '\t' << s.lengthCM << '\t' << s.recombinations << '\t' << s.invisible
               << '\t' << s.nonLocal << '\t' << s.snps << '\t' << s.maxArgEdges << '\t' << s.seconds << '\n';
        lg << "total_seconds\t"
           << std::chrono::duration<double>(std::chrono::steady_clock::now() - tStart).count() << "\n";
    }
    if (!p_.keepTemp) {
        chroms.clear();  // SnpBuffer destructors remove any remaining temporary files
        if (ownTmpDir) {
            std::error_code ec;
            fs::remove_all(tmpDir, ec);
        }
    }
    result.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - tStart).count();
    if (log) {
        *log << "done in " << std::fixed << std::setprecision(1) << result.seconds << " s";
        if (files)
            *log << "; outputs: " << p_.outPrefix << ".{hap.bin,chrom.tsv,snps.tsv,samples.tsv,map.*.tsv"
                 << (p_.ldDecay ? ",ld_decay.tsv" : "") << "}";
        *log << "\n" << std::defaultfloat << std::setprecision(6);
    }
    return result;
}

}  // namespace hapcoal
