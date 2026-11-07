/**
 * @file MacsHistory.hpp
 * @brief (internal) the retained ancestral recombination graph (ARG) of MaCS.
 *
 * MaCS (Chen, Marjoram & Wall 2009, Genome Res 19:136) stores part of the ARG:
 * all edges used by local trees within the last h base pairs. A recombination
 * happens on the *local* tree, but the new lineage may coalesce with *any*
 * retained edge, including non-ancestral ones (edges no longer in the local
 * tree). Material at the current position always follows the most recent
 * edge above a crossover node. Edges unused for more than h bp are pruned.
 * With h = 0 the retained graph is just the local tree plus the branch being
 * replaced, i.e. SMC'.
 *
 * This class keeps the ARG next to a LocalTree. Each ARG edge records its
 * `owner`: the local-tree node whose branch it is part of (-1 if it carries no
 * material at the current position). recombine() returns the local branch and
 * time where the detached lineage really joins the remaining tree, which may be
 * higher than the coalescence point if the lineage hit a non-ancestral edge.
 */
#pragma once
#include <cstdint>
#include <limits>
#include <queue>
#include <string>
#include <vector>

#include "hapcoal/Demography.hpp"
#include "LocalTree.hpp"
#include "Rng.hpp"
#include "TimeIndex.hpp"

namespace hapcoal {

/**
 * @brief Retained ARG of MaCS, kept consistent with a LocalTree.
 *
 * Nodes are samples, coalescences (2 children, 1 parent) and crossovers
 * (1 child, 2 parents: older and newer). Edges carry their owner (local branch)
 * and the last position where they were ancestral; non-ancestral edges older
 * than the history window are removed and redundant nodes fused.
 */
class MacsHistory {
public:
    struct Result {
        int u;           // local node whose (pruned-tree) branch receives the lineage
        double t;        // attachment time in the local tree
        double tCoal;    // time of the coalescence in the ARG (<= t)
        bool nonLocal;   // coalesced with a non-ancestral edge
    };

    /// Copy the local tree into an ARG (no history yet).
    void init(const LocalTree& tree);

    /// Recombination on local branch b at time tr at genomic position pos.
    /// op / s are b's parent and sibling in the local tree *before* the event.
    /// `tree` is the local tree before the event (used for fast edge sampling).
    Result recombine(LocalTree& tree, int b, int op, int s, double tr, double pos, const Demography& demo, Rng& rng);

    /// Remove non-local edges last used before `cutoff` (bp).
    void prune(double cutoff);

    size_t numEdges() const { return edges_.size() - freeEdges_.size(); }
    size_t numNodes() const { return nodes_.size() - freeNodes_.size(); }

    /// Consistency check against the local tree (tests / debug).
    bool validate(const LocalTree& tree, std::string* err = nullptr) const;

    static constexpr double INF = std::numeric_limits<double>::infinity();

private:
    enum : uint8_t { SAMPLE = 0, COAL = 1, XOVER = 2, DEAD = 3 };
    struct Node {
        double t = 0.0;
        uint8_t type = DEAD;
        int bot[2] = {-1, -1};
        int top[2] = {-1, -1};  // XOVER: top[0] older, top[1] newer
    };
    struct Edge {
        int bot = -1, top = -1;  // top == -1: open-ended edge above the grand MRCA
        int owner = -1;
        double last = INF;       // last position (bp) where the edge was ancestral
        uint32_t ver = 0;
        bool alive = false;
    };
    struct HeapItem {
        double last;
        int e;
        uint32_t ver;
        bool operator>(const HeapItem& o) const { return last > o.last; }
    };

    std::vector<Node> nodes_;
    std::vector<Edge> edges_;
    std::vector<int> freeNodes_, freeEdges_;
    std::vector<int> argOf_;  // local node -> ARG node
    int extEdge_ = -1;        // the edge with open top
    TimeIndex times_;
    std::priority_queue<HeapItem, std::vector<HeapItem>, std::greater<HeapItem>> heap_;
    std::vector<int> stack_, buf_, work_;
    // Non-ancestral edges (owner == -1) in a dense list for O(1) uniform sampling.
    std::vector<int> nonLocal_;
    std::vector<int> nlPos_;          // edge -> index in nonLocal_ (-1 if absent)
    std::vector<uint32_t> oldChain_;  // edge -> event stamp if on b's replaced path
    uint32_t event_ = 0;
    void setOwner(int e, int owner);

    double topTime(int e) const { return edges_[e].top < 0 ? INF : nodes_[edges_[e].top].t; }
    double botTime(int e) const { return nodes_[edges_[e].bot].t; }
    static int deltaOf(uint8_t type) { return type == COAL ? -1 : +1; }

    // ---- allocation (free lists; ids are reused)
    int newNode(double t, uint8_t type);
    void freeNode(int v);
    int newEdge(int bot, int top, int owner, double last);
    void freeEdge(int e);
    /// Mark edge @p e non-ancestral since position @p pos (enters the pruning heap).
    void setNonLocal(int e, double pos);
    int upEdge(int e) const;  // next edge on the material path above e
    /// Replace a child (bottom) edge slot of @p node.
    void replaceBot(int node, int oldE, int newE);
    /// Replace a parent (top) edge slot of @p node.
    void replaceTop(int node, int oldE, int newE);
    int pickEdgeCrossing(double t, int K, LocalTree& tree, Rng& rng);
    int chainEdgeAt(int v, double t) const;  // ARG edge of local branch v crossing t

    // ---- pruning
    /// Detach and free edge @p e; its end nodes are queued for clean-up.
    void removeEdge(int e);
    /// Clean up queued nodes: fuse pass-through nodes, cascade removal of dead branches.
    void fixNodes();
    /// Remove a node with one child and one parent edge by joining the two edges.
    void fuse(int v);
};

}  // namespace hapcoal
