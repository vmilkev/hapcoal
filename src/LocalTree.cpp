/**
 * @file LocalTree.cpp
 * @brief (internal) Local tree: Kingman initial tree, Fenwick sampling, prune-and-regraft, leaf-set cache.
 */
#include "LocalTree.hpp"

#include <cassert>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

namespace hapcoal {

LocalTree::LocalTree(int nSamples) : n_(nSamples) {
    if (n_ < 2) throw std::invalid_argument("LocalTree: need at least 2 samples");
    const int N = numNodes();
    parent_.assign(N, -1);
    left_.assign(N, -1);
    right_.assign(N, -1);
    time_.assign(N, 0.0);
    fw_.assign(N + 1, 0.0);
    w_.assign(N, 0.0);
    fwLog_ = 1;
    while ((1 << fwLog_) <= N) ++fwLog_;
    stack_.reserve(N);
    crossBuf_.reserve(N);
}

// Kingman coalescent: with k lineages the next coalescence happens at hazard
// C(k,2)/lambda(t); a uniformly chosen pair merges into a new node.
void LocalTree::buildInitial(const Demography& demo, Rng& rng) {
    const int N = numNodes();
    std::fill(parent_.begin(), parent_.end(), -1);
    std::fill(left_.begin(), left_.end(), -1);
    std::fill(right_.begin(), right_.end(), -1);
    std::fill(time_.begin(), time_.end(), 0.0);
    times_.clear();

    std::vector<int> lineages(n_);
    for (int i = 0; i < n_; ++i) lineages[i] = i;
    double t = 0.0;
    int next = n_;
    while (lineages.size() > 1) {
        const double k = static_cast<double>(lineages.size());
        t = demo.sampleWaiting(t, 0.5 * k * (k - 1.0), rng);
        const size_t i = rng.below(lineages.size());
        size_t j = rng.below(lineages.size() - 1);
        if (j >= i) ++j;
        const int a = lineages[i], b = lineages[j];
        const int p = next++;
        time_[p] = t;
        left_[p] = a;
        right_[p] = b;
        parent_[a] = p;
        parent_[b] = p;
        times_.insert(t, -1, p);
        // remove i and j, add p
        const size_t hi = std::max(i, j), lo = std::min(i, j);
        lineages[hi] = lineages.back();
        lineages.pop_back();
        lineages[lo] = p;
    }
    root_ = lineages[0];
    assert(next == N);
    (void)N;
    rebuildFenwick();
    if (cacheOn_)
        for (int v = n_; v < N; ++v) valid_[v] = 0;
}

// Rebuild the Fenwick (binary indexed) tree over branch lengths from scratch in
// O(N); done initially and periodically to remove floating-point drift.
void LocalTree::rebuildFenwick() {
    const int N = numNodes();
    std::fill(fw_.begin(), fw_.end(), 0.0);
    total_ = 0.0;
    for (int v = 0; v < N; ++v) {
        w_[v] = branchLength(v);
        total_ += w_[v];
        fw_[v + 1] += w_[v];
    }
    for (int i = 1; i <= N; ++i) {
        const int j = i + (i & -i);
        if (j <= N) fw_[j] += fw_[i];
    }
    updatesSinceRebuild_ = 0;
}

// Update the Fenwick weight of node v after its branch length changed (O(log N)).
void LocalTree::refresh(int v) {
    const double nw = branchLength(v);
    const double d = nw - w_[v];
    if (d == 0.0) return;
    w_[v] = nw;
    total_ += d;
    const int N = numNodes();
    for (int i = v + 1; i <= N; i += i & -i) fw_[i] += d;
    if (++updatesSinceRebuild_ > (1u << 20)) rebuildFenwick();
}

// Fenwick descent: the smallest index whose prefix sum of branch lengths exceeds u,
// i.e. the branch containing the point u of the tree laid out branch by branch.
int LocalTree::branchAt(double u) const {
    const int N = numNodes();
    int pos = 0;
    for (int step = 1 << (fwLog_ - 1); step > 0; step >>= 1) {
        const int nxt = pos + step;
        if (nxt <= N && fw_[nxt] <= u) {
            pos = nxt;
            u -= fw_[nxt];
        }
    }
    int v = pos < N ? pos : N - 1;
    // Guard against rounding landing on a zero-length branch (the root).
    if (w_[v] <= 0.0) {
        int a = v;
        while (a >= 0 && w_[a] <= 0.0) --a;
        if (a >= 0) return a;
        a = v;
        while (a < N && w_[a] <= 0.0) ++a;
        v = a;
    }
    return v;
}

// Uniform choice among the K branches crossing t. If they are a large fraction of
// all branches (K*8 >= N) rejection sampling over all nodes is O(1) expected;
// otherwise the K branches are enumerated by descending from the root (O(K)).
int LocalTree::pickBranchCrossing(double t, int exclude, Rng& rng) {
    if (t >= time_[root_]) return root_;
    const int K = lineagesAt(t);
    const int N = numNodes();
    if (static_cast<int64_t>(K) * 8 >= N) {
        for (;;) {
            const int v = static_cast<int>(rng.below(static_cast<uint64_t>(N)));
            if (v != exclude && crosses(v, t)) return v;
        }
    }
    crossBuf_.clear();
    stack_.clear();
    stack_.push_back(root_);
    while (!stack_.empty()) {
        const int w = stack_.back();
        stack_.pop_back();
        for (int c : {left_[w], right_[w]}) {
            if (time_[c] <= t) {
                if (c != exclude) crossBuf_.push_back(c);
            } else {
                stack_.push_back(c);
            }
        }
    }
    if (crossBuf_.empty()) throw std::logic_error("pickBranchCrossing: no branch crosses t");
    return crossBuf_[rng.below(crossBuf_.size())];
}

int LocalTree::regraft(int b, int u, double t) {
    const int op = parent_[b];
    const int s = sibling(b);
    const int gp = parent_[op];
    assert(u != op && u != b);
    if (cacheOn_) invalidateUp(op);  // op and its old ancestors lose b's leaves

    // --- prune: remove op, connect s to gp
    if (gp >= 0) {
        if (left_[gp] == op) left_[gp] = s;
        else right_[gp] = s;
        parent_[s] = gp;
    } else {
        parent_[s] = -1;
        root_ = s;
    }
    times_.erase(time_[op], op);

    // --- regraft: new node (reusing id op) on the branch above u at time t
    const int p = parent_[u];
    assert(t >= time_[u] && t >= time_[b]);
    assert(p < 0 || t <= time_[p]);
    time_[op] = t;
    left_[op] = b;
    right_[op] = u;
    parent_[op] = p;
    parent_[b] = op;
    parent_[u] = op;
    if (p >= 0) {
        if (left_[p] == u) left_[p] = op;
        else right_[p] = op;
    } else {
        root_ = op;
    }
    times_.insert(t, -1, op);
    if (cacheOn_) {
        valid_[op] = 0;
        invalidateUp(p);  // new ancestors gain b's leaves
    }

    refresh(b);
    refresh(s);
    refresh(u);
    refresh(op);
    return op;
}

bool LocalTree::enableLeafCache(size_t maxBytes) {
    words_ = (static_cast<size_t>(n_) + 63) / 64;
    const size_t need = static_cast<size_t>(numNodes()) * words_ * 8;
    if (need > maxBytes) { cacheOn_ = false; return false; }
    cacheOn_ = true;
    bits_.assign(static_cast<size_t>(numNodes()) * words_, 0);
    nLeaves_.assign(numNodes(), 0);
    valid_.assign(numNodes(), 0);
    for (int i = 0; i < n_; ++i) {
        bits_[static_cast<size_t>(i) * words_ + (i >> 6)] = 1ULL << (i & 63);
        nLeaves_[i] = 1;
        valid_[i] = 1;
    }
    return true;
}

void LocalTree::ensureValid(int v) {
    if (valid_[v]) return;
    // post-order over the invalid part of the subtree
    cstack_.clear();
    cstack_.push_back(v);
    while (!cstack_.empty()) {
        const int w = cstack_.back();
        const int a = left_[w], b = right_[w];
        if (!valid_[a]) { cstack_.push_back(a); continue; }
        if (!valid_[b]) { cstack_.push_back(b); continue; }
        uint64_t* __restrict d = bits_.data() + static_cast<size_t>(w) * words_;
        const uint64_t* __restrict x = bits_.data() + static_cast<size_t>(a) * words_;
        const uint64_t* __restrict y = bits_.data() + static_cast<size_t>(b) * words_;
        for (size_t k = 0; k < words_; ++k) d[k] = x[k] | y[k];
        nLeaves_[w] = nLeaves_[a] + nLeaves_[b];
        valid_[w] = 1;
        cstack_.pop_back();
    }
}

uint32_t LocalTree::leafBits(int v, uint64_t* row) {
    if (cacheOn_) {
        ensureValid(v);
        const uint64_t* s = bits_.data() + static_cast<size_t>(v) * words_;
        for (size_t k = 0; k < words_; ++k) row[k] = s[k];
        return nLeaves_[v];
    }
    uint32_t c = 0;
    stack_.clear();
    stack_.push_back(v);
    while (!stack_.empty()) {
        const int w = stack_.back();
        stack_.pop_back();
        if (w < n_) {
            row[w >> 6] |= 1ULL << (w & 63);
            ++c;
        } else {
            stack_.push_back(right_[w]);
            stack_.push_back(left_[w]);
        }
    }
    return c;
}

void LocalTree::collectLeaves(int v, std::vector<int>& out) const {
    stack_.clear();
    stack_.push_back(v);
    while (!stack_.empty()) {
        const int w = stack_.back();
        stack_.pop_back();
        if (w < n_) {
            out.push_back(w);
        } else {
            stack_.push_back(right_[w]);
            stack_.push_back(left_[w]);
        }
    }
}

bool LocalTree::validate(std::string* err) const {
    std::ostringstream os;
    const int N = numNodes();
    int roots = 0;
    for (int v = 0; v < N; ++v) {
        if (parent_[v] < 0) {
            ++roots;
            if (v != root_) os << "node " << v << " has no parent but is not root\n";
        } else {
            const int p = parent_[v];
            if (left_[p] != v && right_[p] != v) os << "parent " << p << " does not list child " << v << "\n";
            if (!(time_[p] >= time_[v])) os << "time inversion at " << v << "\n";
        }
        if (v >= n_) {
            if (left_[v] < 0 || right_[v] < 0) os << "internal node " << v << " lacks children\n";
            else if (parent_[left_[v]] != v || parent_[right_[v]] != v) os << "child link broken at " << v << "\n";
        }
        if (std::fabs(w_[v] - branchLength(v)) > 1e-9 * (1 + std::fabs(w_[v]))) os << "fenwick weight stale at " << v << "\n";
    }
    if (roots != 1) os << "number of roots = " << roots << "\n";
    std::vector<int> leaves;
    collectLeaves(root_, leaves);
    if (static_cast<int>(leaves.size()) != n_) os << "root reaches " << leaves.size() << " leaves\n";
    if (static_cast<int>(times_.size()) != n_ - 1) os << "time index has " << times_.size() << " items\n";
    double tot = 0;
    for (int v = 0; v < N; ++v) tot += branchLength(v);
    if (std::fabs(tot - total_) > 1e-7 * (1 + tot)) os << "total length drift " << tot << " vs " << total_ << "\n";
    if (err) *err = os.str();
    return os.str().empty();
}

}  // namespace hapcoal
