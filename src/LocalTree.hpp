/**
 * @file LocalTree.hpp
 * @brief (internal) the marginal genealogy at the current position.
 *
 * Nodes 0..n-1 are samples (time 0); nodes n..2n-2 are coalescence nodes.
 * Node ids are stable: a prune-and-regraft reuses the removed parent's id
 * for the new coalescence node, so the tree always has exactly 2n-1 nodes.
 *
 * Data structures:
 * * SoA arrays parent/left/right/time
 * * Fenwick tree over branch lengths -> O(log n) sampling of a point on the
 * tree proportional to length (recombination and mutation placement)
 * * TimeIndex of internal-node times -> lineage count K(t) and ordered
 * traversal of coalescence times above t
 */
#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "hapcoal/Demography.hpp"
#include "Rng.hpp"
#include "TimeIndex.hpp"

namespace hapcoal {

/**
 * @brief Marginal genealogy (local tree) of n samples at the current position.
 *
 * Supports the operations of the sequential process: sampling a point on the
 * tree proportionally to branch length, counting and enumerating the lineages
 * crossing a time, and moving a subtree (prune and regraft).
 */
class LocalTree {
public:
    /// @param nSamples number of sampled haplotypes n (>= 2)
    explicit LocalTree(int nSamples);

    /// Kingman coalescent under the given demography.
    void buildInitial(const Demography& demo, Rng& rng);

    int numSamples() const { return n_; }
    int numNodes() const { return 2 * n_ - 1; }
    int root() const { return root_; }
    int parent(int v) const { return parent_[v]; }
    int left(int v) const { return left_[v]; }
    int right(int v) const { return right_[v]; }
    double time(int v) const { return time_[v]; }
    int sibling(int v) const {
        const int p = parent_[v];
        return left_[p] == v ? right_[p] : left_[p];
    }
    bool isLeaf(int v) const { return v < n_; }

    /// Total branch length (coalescent units).
    double totalLength() const { return total_; }
    /// Node whose branch (to its parent) contains the point `u` in [0, totalLength()).
    int branchAt(double u) const;

    /// Number of lineages of the tree crossing time t (K = 1 above the root).
    int lineagesAt(double t) const { return n_ + times_.prefix(t); }
    const TimeIndex& times() const { return times_; }

    /// Branch crossing time t (node v with time(v) <= t < time(parent(v)), or the
    /// root if t >= time(root)), chosen uniformly; `exclude` is never returned.
    int pickBranchCrossing(double t, int exclude, Rng& rng);

    /// Prune the branch above b (removing its parent op) and regraft it onto the
    /// branch above u at time t. `u` refers to the tree after pruning; passing
    /// the old parent itself is not allowed (map it to the sibling first).
    /// Returns the id of the new coalescence node (== old parent id).
    int regraft(int b, int u, double t);

    /// Append all leaves below v to `out`.
    void collectLeaves(int v, std::vector<int>& out) const;

    /// Enable the lazy per-node leaf-set cache (bitsets of nSamples bits) if it
    /// needs at most maxBytes. Returns whether it is enabled.
    bool enableLeafCache(size_t maxBytes);
    /// Write the leaf set of v as a bit row (words = ceil(n/64)); returns the number of leaves.
    uint32_t leafBits(int v, uint64_t* row);

    /// Full consistency check (debug / tests). Returns false on error.
    bool validate(std::string* err = nullptr) const;

private:
    int n_;
    int root_ = -1;
    std::vector<int> parent_, left_, right_;
    std::vector<double> time_;
    TimeIndex times_;

    // Fenwick tree over branch lengths
    std::vector<double> fw_, w_;
    double total_ = 0.0;
    int fwLog_ = 0;
    uint64_t updatesSinceRebuild_ = 0;
    mutable std::vector<int> stack_;
    std::vector<int> crossBuf_;
    // Leaf-set cache. Invariant: a valid node has only valid descendants, so an
    // invalid node always has invalid ancestors (invalidation can stop early).
    bool cacheOn_ = false;
    size_t words_ = 0;
    std::vector<uint64_t> bits_;
    std::vector<uint32_t> nLeaves_;
    std::vector<uint8_t> valid_;
    std::vector<int> cstack_;
    void invalidateUp(int v) {
        while (v >= 0 && valid_[v]) { valid_[v] = 0; v = parent_[v]; }
    }
    void ensureValid(int v);

    double branchLength(int v) const { return parent_[v] < 0 ? 0.0 : time_[parent_[v]] - time_[v]; }
    void refresh(int v);
    void rebuildFenwick();
    bool crosses(int v, double t) const {
        return time_[v] <= t && (parent_[v] < 0 || time_[parent_[v]] > t);
    }
};

}  // namespace hapcoal
