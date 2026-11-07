/**
 * @file MacsHistory.cpp
 * @brief (internal) MaCS retained ARG: recombination, coalescence with retained edges, pruning.
 */
#include "MacsHistory.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace hapcoal {

// ------------------------------------------------------------ allocation ----
int MacsHistory::newNode(double t, uint8_t type) {
    int v;
    if (!freeNodes_.empty()) {
        v = freeNodes_.back();
        freeNodes_.pop_back();
    } else {
        v = static_cast<int>(nodes_.size());
        nodes_.emplace_back();
    }
    Node& n = nodes_[v];
    n.t = t;
    n.type = type;
    n.bot[0] = n.bot[1] = n.top[0] = n.top[1] = -1;
    return v;
}

void MacsHistory::freeNode(int v) {
    Node& n = nodes_[v];
    if (n.type == DEAD) return;
    times_.erase(n.t, v);
    n.type = DEAD;
    freeNodes_.push_back(v);
}

int MacsHistory::newEdge(int bot, int top, int owner, double last) {
    int e;
    if (!freeEdges_.empty()) {
        e = freeEdges_.back();
        freeEdges_.pop_back();
    } else {
        e = static_cast<int>(edges_.size());
        edges_.emplace_back();
        nlPos_.push_back(-1);
        oldChain_.push_back(0);
    }
    Edge& E = edges_[e];
    E.bot = bot;
    E.top = top;
    E.owner = -2;  // forces list bookkeeping in setOwner
    E.last = last;
    ++E.ver;
    E.alive = true;
    oldChain_[e] = 0;
    setOwner(e, owner);
    if (top < 0) extEdge_ = e;
    return e;
}

void MacsHistory::setOwner(int e, int owner) {
    Edge& E = edges_[e];
    if (E.owner == owner) return;
    if (owner == -1 && nlPos_[e] < 0) {
        nlPos_[e] = static_cast<int>(nonLocal_.size());
        nonLocal_.push_back(e);
    } else if (owner != -1 && nlPos_[e] >= 0) {
        const int i = nlPos_[e], last = nonLocal_.back();
        nonLocal_[i] = last;
        nlPos_[last] = i;
        nonLocal_.pop_back();
        nlPos_[e] = -1;
    }
    E.owner = owner;
}

void MacsHistory::freeEdge(int e) {
    Edge& E = edges_[e];
    setOwner(e, -3);  // leave the non-local list
    E.alive = false;
    ++E.ver;
    freeEdges_.push_back(e);
}

void MacsHistory::setNonLocal(int e, double pos) {
    Edge& E = edges_[e];
    setOwner(e, -1);
    E.last = pos;
    ++E.ver;
    heap_.push({pos, e, E.ver});
}

int MacsHistory::upEdge(int e) const {
    const Node& T = nodes_[edges_[e].top];
    if (T.type == XOVER) return T.top[1] >= 0 ? T.top[1] : T.top[0];
    return T.top[0];
}

void MacsHistory::replaceBot(int node, int oldE, int newE) {
    Node& n = nodes_[node];
    if (n.bot[0] == oldE) n.bot[0] = newE;
    else if (n.bot[1] == oldE) n.bot[1] = newE;
    else throw std::logic_error("MacsHistory::replaceBot: edge not found");
}

void MacsHistory::replaceTop(int node, int oldE, int newE) {
    Node& n = nodes_[node];
    if (n.top[0] == oldE) n.top[0] = newE;
    else if (n.top[1] == oldE) n.top[1] = newE;
    else throw std::logic_error("MacsHistory::replaceTop: edge not found");
}

// ------------------------------------------------------------------ init ----
void MacsHistory::init(const LocalTree& tree) {
    nodes_.clear();
    edges_.clear();
    freeNodes_.clear();
    freeEdges_.clear();
    nonLocal_.clear();
    nlPos_.clear();
    oldChain_.clear();
    times_.clear();
    heap_ = decltype(heap_)();
    const int N = tree.numNodes();
    argOf_.assign(N, -1);
    for (int v = 0; v < N; ++v) {
        const int a = newNode(tree.time(v), tree.isLeaf(v) ? SAMPLE : COAL);
        assert(a == v);
        argOf_[v] = a;
        times_.insert(tree.time(v), deltaOf(nodes_[a].type), a);
    }
    for (int v = 0; v < N; ++v) {
        const int p = tree.parent(v);
        const int e = newEdge(v, p, v, INF);  // p == -1 creates the open-ended edge
        nodes_[v].top[0] = e;
        if (p >= 0) {
            Node& P = nodes_[p];
            if (P.bot[0] < 0) P.bot[0] = e;
            else P.bot[1] = e;
        }
    }
}

// ------------------------------------------------------- edge selection ----
int MacsHistory::chainEdgeAt(int v, double t) const {
    int e = nodes_[argOf_[v]].top[0];
    while (topTime(e) <= t) e = upEdge(e);
    return e;
}

// Uniform choice among the K retained edges crossing t. Every such edge lies
// either on the path of a branch of the (pre-event) local tree - K_tree(t) of
// them, sampled through the LocalTree in O(log n) - or is a purely non-ancestral
// edge (sampled from the dense non-local list, excluding b's replaced path,
// which is already represented by b's local branch).
int MacsHistory::pickEdgeCrossing(double t, int K, LocalTree& tree, Rng& rng) {
    const int Kt = tree.lineagesAt(t);
    const int Ko = K - Kt;
    if (Ko < 0) throw std::logic_error("MacsHistory: lineage count below local tree count");
    if (Ko == 0 || static_cast<int>(rng.below(static_cast<uint64_t>(K))) < Kt)
        return chainEdgeAt(tree.pickBranchCrossing(t, -1, rng), t);
    const auto other = [&](int e) {
        const Edge& E = edges_[e];
        return oldChain_[e] != event_ && nodes_[E.bot].t <= t && t < topTime(e);
    };
    const size_t N = nonLocal_.size();
    if (N <= static_cast<size_t>(Ko) * 32) {
        for (;;) {
            const int e = nonLocal_[rng.below(N)];
            if (other(e)) return e;
        }
    }
    buf_.clear();
    for (int e : nonLocal_)
        if (other(e)) buf_.push_back(e);
    if (static_cast<int>(buf_.size()) != Ko)
        throw std::logic_error("MacsHistory: non-local lineage count mismatch (" + std::to_string(buf_.size()) +
                               " vs " + std::to_string(Ko) + ")");
    return buf_[rng.below(buf_.size())];
}

// ---------------------------------------------------------- recombination ----
MacsHistory::Result MacsHistory::recombine(LocalTree& tree, int b, int op, int s, double tr, double pos,
                                           const Demography& demo, Rng& rng) {
    ++event_;
    // 1. Edge of b's material path that contains tr.
    int e = nodes_[argOf_[b]].top[0];
    while (topTime(e) <= tr) e = upEdge(e);
    assert(edges_[e].owner == b);

    // 2. Insert the crossover node x on e.
    const int T0 = edges_[e].top;
    const int x = newNode(tr, XOVER);
    const int eup = newEdge(x, T0, b, INF);
    if (T0 >= 0) replaceBot(T0, e, eup);
    edges_[e].top = x;
    nodes_[x].bot[0] = e;
    nodes_[x].top[0] = eup;

    // 3. The old path above x no longer carries material: it becomes non-local
    //    up to A(op), where it met the sibling's material.
    int cur = eup;
    int Aop = -1;
    for (;;) {
        setNonLocal(cur, pos);
        oldChain_[cur] = event_;
        const int Tn = edges_[cur].top;
        if (Tn < 0) throw std::logic_error("MacsHistory: b path reached the open edge");
        const Node& N = nodes_[Tn];
        if (N.type == COAL) {
            const int other = N.bot[0] == cur ? N.bot[1] : N.bot[0];
            if (edges_[other].owner != -1) {
                Aop = Tn;
                break;
            }
            cur = N.top[0];
        } else {
            cur = N.top[1] >= 0 ? N.top[1] : N.top[0];
        }
    }
    // op is no longer a coalescence node of the local tree: its path now belongs to s.
    for (cur = nodes_[Aop].top[0]; cur >= 0 && edges_[cur].owner == op;) {
        setOwner(cur, s);
        if (edges_[cur].top < 0) break;
        cur = upEdge(cur);
    }

    // 4. Coalescence time: hazard k(t)/lambda(t), k = number of retained edges at t
    //    (x is not yet in the time index, so the new lineage is not counted).
    double t = tr;
    int k = times_.prefix(tr);
    TimeIndex::Cursor c = times_.upperBound(tr);
    double E = rng.exp1();
    double tc;
    for (;;) {
        const double tn = times_.valid(c) ? times_.get(c).t : INF;
        const double te = demo.nextBoundary(t);
        const double tb = std::min(tn, te);
        const double rate = k / demo.lambda(t);
        if (E < rate * (tb - t)) {
            tc = t + E / rate;
            break;
        }
        E -= rate * (tb - t);
        t = tb;
        if (tb == tn) {
            k += times_.get(c).delta;
            times_.advance(c);
        }
    }

    // 5. Target edge, chosen uniformly among the retained edges crossing tc.
    const int tgt = pickEdgeCrossing(tc, k, tree, rng);

    // 6. Insert the coalescence node cn on tgt and connect the new lineage F: x -> cn.
    const int Tt = edges_[tgt].top;
    const int cn = newNode(tc, COAL);
    const int ehigh = newEdge(cn, Tt, edges_[tgt].owner, edges_[tgt].last);
    if (Tt >= 0) replaceBot(Tt, tgt, ehigh);
    edges_[tgt].top = cn;
    const int F = newEdge(x, cn, b, INF);
    nodes_[x].top[1] = F;
    nodes_[cn].bot[0] = tgt;
    nodes_[cn].bot[1] = F;
    nodes_[cn].top[0] = ehigh;
    times_.insert(tr, +1, x);
    times_.insert(tc, -1, cn);
    if (edges_[ehigh].owner == -1) heap_.push({edges_[ehigh].last, ehigh, edges_[ehigh].ver});

    Result res{-1, tc, tc, false};
    int meet, aboveMeet;
    if (edges_[tgt].owner != -1) {
        // Coalesced directly with an ancestral edge: ordinary prune-and-regraft.
        res.u = edges_[tgt].owner;
        meet = cn;
        aboveMeet = ehigh;
    } else {
        // Coalesced with a non-ancestral edge: the material follows that edge's
        // path upwards (newest edge at each crossover) until it meets material
        // of the remaining sample lineages.
        res.nonLocal = true;
        cur = ehigh;
        setOwner(cur, b);
        edges_[cur].last = INF;
        ++edges_[cur].ver;
        for (;;) {
            const int Tn = edges_[cur].top;
            if (Tn < 0) throw std::logic_error("MacsHistory: non-local path reached the open edge");
            const Node& N = nodes_[Tn];
            if (N.type == COAL) {
                const int other = N.bot[0] == cur ? N.bot[1] : N.bot[0];
                if (edges_[other].owner != -1) {
                    res.u = edges_[other].owner;
                    res.t = N.t;
                    meet = Tn;
                    aboveMeet = N.top[0];
                    break;
                }
                cur = N.top[0];
            } else {
                cur = N.top[1] >= 0 ? N.top[1] : N.top[0];
            }
            assert(edges_[cur].owner == -1);
            setOwner(cur, b);
            edges_[cur].last = INF;
            ++edges_[cur].ver;
        }
    }
    // 7. Above the meeting point the path of u now belongs to the new local node (id op).
    for (cur = aboveMeet; cur >= 0 && edges_[cur].owner == res.u;) {
        setOwner(cur, op);
        if (edges_[cur].top < 0) break;
        cur = upEdge(cur);
    }
    argOf_[op] = meet;
    return res;
}

// --------------------------------------------------------------- pruning ----
void MacsHistory::prune(double cutoff) {
    while (!heap_.empty() && heap_.top().last < cutoff) {
        const HeapItem h = heap_.top();
        heap_.pop();
        const Edge& E = edges_[h.e];
        if (!E.alive || E.ver != h.ver || E.owner != -1) continue;
        removeEdge(h.e);
        fixNodes();
    }
}

void MacsHistory::removeEdge(int e) {
    Edge& E = edges_[e];
    assert(E.alive);
    assert(E.owner == -1 && "MacsHistory: attempted to prune an ancestral edge");
    const int B = E.bot, T = E.top;
    {
        Node& nb = nodes_[B];
        if (nb.top[0] == e) { nb.top[0] = nb.top[1]; nb.top[1] = -1; }
        else if (nb.top[1] == e) nb.top[1] = -1;
    }
    if (T >= 0) {
        Node& nt = nodes_[T];
        if (nt.bot[0] == e) { nt.bot[0] = nt.bot[1]; nt.bot[1] = -1; }
        else if (nt.bot[1] == e) nt.bot[1] = -1;
        work_.push_back(T);
    }
    freeEdge(e);
    work_.push_back(B);
}

void MacsHistory::fuse(int v) {
    Node& N = nodes_[v];
    const int a = N.bot[0] >= 0 ? N.bot[0] : N.bot[1];
    const int c = N.top[0] >= 0 ? N.top[0] : N.top[1];
    const int T = edges_[c].top;
    Edge& A = edges_[a];
    const Edge& C = edges_[c];
    A.top = T;
    if (T >= 0) replaceBot(T, c, a);
    else extEdge_ = a;
    if (A.owner == -1) setOwner(a, C.owner);
    A.last = std::max(A.last, C.last);
    ++A.ver;
    if (A.owner == -1) heap_.push({A.last, a, A.ver});
    freeEdge(c);
    freeNode(v);
}

void MacsHistory::fixNodes() {
    while (!work_.empty()) {
        const int v = work_.back();
        work_.pop_back();
        Node& N = nodes_[v];
        if (N.type == DEAD || N.type == SAMPLE) continue;
        const int nb = (N.bot[0] >= 0) + (N.bot[1] >= 0);
        const int nt = (N.top[0] >= 0) + (N.top[1] >= 0);
        const int wantTop = N.type == XOVER ? 2 : 1;
        const int wantBot = N.type == XOVER ? 1 : 2;
        if (nb == wantBot && nt == wantTop) continue;
        if (nb == 1 && nt == 1) {
            fuse(v);
            continue;
        }
        if (nb == 0) {
            const int t0 = N.top[0], t1 = N.top[1];
            freeNode(v);
            if (t0 >= 0) removeEdge(t0);
            if (t1 >= 0) removeEdge(t1);
            continue;
        }
        if (nt == 0) {
            const int b0 = N.bot[0], b1 = N.bot[1];
            freeNode(v);
            if (b0 >= 0) removeEdge(b0);
            if (b1 >= 0) removeEdge(b1);
            continue;
        }
    }
}

// ------------------------------------------------------------ validation ----
bool MacsHistory::validate(const LocalTree& tree, std::string* err) const {
    std::ostringstream os;
    // Each local branch must be a chain of edges owned by that node.
    for (int v = 0; v < tree.numNodes(); ++v) {
        const int a = argOf_[v];
        if (a < 0 || nodes_[a].type == DEAD) { os << "local node " << v << " has no ARG node\n"; continue; }
        if (std::fabs(nodes_[a].t - tree.time(v)) > 1e-12) os << "time mismatch at local node " << v << "\n";
        int e = nodes_[a].top[0];
        int steps = 0;
        const int p = tree.parent(v);
        while (true) {
            if (e < 0 || !edges_[e].alive) { os << "broken chain above local node " << v << "\n"; break; }
            if (edges_[e].owner != v) { os << "edge owner " << edges_[e].owner << " != " << v << "\n"; break; }
            if (p < 0) { if (edges_[e].top < 0) break; }
            else if (edges_[e].top == argOf_[p]) break;
            if (edges_[e].top < 0) { os << "chain of " << v << " ends at open edge\n"; break; }
            e = upEdge(e);
            if (++steps > 1000000) { os << "cycle above " << v << "\n"; break; }
        }
    }
    // Lineage count from the time index must equal the number of edges crossing t.
    for (double t : {0.01, 0.1, 0.5, 1.0, 2.0}) {
        int cnt = 0;
        for (size_t e = 0; e < edges_.size(); ++e)
            if (edges_[e].alive && nodes_[edges_[e].bot].t <= t && t < topTime(static_cast<int>(e))) ++cnt;
        if (cnt != times_.prefix(t)) os << "count mismatch at t=" << t << ": " << cnt << " vs " << times_.prefix(t) << "\n";
    }
    if (err) *err = os.str();
    return os.str().empty();
}

}  // namespace hapcoal
