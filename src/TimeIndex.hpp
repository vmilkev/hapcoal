/**
 * @file TimeIndex.hpp
 * @brief (internal) ordered multiset of node events (time, delta, id).
 *
 * The number of lineages crossing time t in a genealogy (local tree or ARG)
 * is a step function: every node changes it by a fixed delta (sample +1,
 * coalescence -1, crossover +1). TimeIndex stores the node times in sorted
 * blocks with per-block delta sums so that
 * * insert / erase              O(B + the number of blocks)
 * * prefix(t) = sum of deltas   O(the number of blocks + B)
 * * forward iteration from t    O(1) per step
 * with block size B ~ 256 this is fast for tens of thousands of nodes.
 */
#pragma once
#include <algorithm>
#include <cassert>
#include <limits>
#include <vector>

namespace hapcoal {

/// Ordered multiset of node events with per-block delta sums (see file description).
class TimeIndex {
public:
    /// One node event.
    struct Item {
        double t;   ///< node time (coalescent units)
        int delta;  ///< change of the lineage count at t (+1 sample / crossover, -1 coalescence)
        int id;     ///< node id (to find the item again on erase)
    };
    /// Position of an item: block b, index i inside the block.
    struct Cursor {
        size_t b = 0, i = 0;
    };

    void clear() { blocks_.clear(); size_ = 0; }
    size_t size() const { return size_; }

    /// Insert an event (O(B + log the number of blocks)); splits a block that grows beyond 2B.
    void insert(double t, int delta, int id) {
        if (blocks_.empty()) blocks_.emplace_back();
        size_t b = blockFor(t);
        auto& it = blocks_[b].items;
        auto pos = std::upper_bound(it.begin(), it.end(), t, [](double v, const Item& x) { return v < x.t; });
        it.insert(pos, Item{t, delta, id});
        blocks_[b].sum += delta;
        ++size_;
        if (it.size() > 2 * kBlock) split(b);
    }

    /// Remove the item with the given (time, id). Returns false if not found.
    bool erase(double t, int id) {
        if (blocks_.empty()) return false;
        size_t b = blockFor(t);
        // Equal times may straddle blocks: search this and neighbouring blocks.
        for (size_t bb = (b > 0 ? b - 1 : 0); bb < blocks_.size() && bb <= b + 1; ++bb) {
            auto& it = blocks_[bb].items;
            auto lo = std::lower_bound(it.begin(), it.end(), t, [](const Item& x, double v) { return x.t < v; });
            for (auto p = lo; p != it.end() && p->t == t; ++p)
                if (p->id == id) {
                    blocks_[bb].sum -= p->delta;
                    it.erase(p);
                    --size_;
                    if (it.empty() && blocks_.size() > 1) blocks_.erase(blocks_.begin() + static_cast<long>(bb));
                    return true;
                }
        }
        return false;
    }

    /// Sum of deltas of all items with time <= t (the lineage count K(t) up to a constant).
    int prefix(double t) const {
        int s = 0;
        for (const auto& bl : blocks_) {
            if (bl.items.empty()) continue;
            if (bl.items.back().t <= t) { s += bl.sum; continue; }
            for (const auto& x : bl.items) {
                if (x.t > t) break;
                s += x.delta;
            }
            break;
        }
        return s;
    }

    /// Cursor at the first item with time > t.
    Cursor upperBound(double t) const {
        Cursor c;
        for (c.b = 0; c.b < blocks_.size(); ++c.b) {
            const auto& it = blocks_[c.b].items;
            if (!it.empty() && it.back().t > t) {
                c.i = static_cast<size_t>(std::upper_bound(it.begin(), it.end(), t,
                                                           [](double v, const Item& x) { return v < x.t; }) -
                                          it.begin());
                return c;
            }
        }
        c.i = 0;
        return c;
    }
    bool valid(const Cursor& c) const { return c.b < blocks_.size() && c.i < blocks_[c.b].items.size(); }
    const Item& get(const Cursor& c) const { return blocks_[c.b].items[c.i]; }
    void advance(Cursor& c) const {
        if (++c.i >= blocks_[c.b].items.size()) {
            ++c.b;
            c.i = 0;
            while (c.b < blocks_.size() && blocks_[c.b].items.empty()) ++c.b;
        }
    }

    template <class F>
    void forEach(F f) const {
        for (const auto& bl : blocks_)
            for (const auto& x : bl.items) f(x);
    }

private:
    static constexpr size_t kBlock = 256;
    struct Block {
        std::vector<Item> items;
        int sum = 0;
    };
    std::vector<Block> blocks_;
    size_t size_ = 0;

    // Block whose range should contain t: first block whose last time >= t, else the last block.
    size_t blockFor(double t) const {
        size_t lo = 0, hi = blocks_.size();
        while (lo < hi) {
            size_t mid = (lo + hi) / 2;
            const auto& it = blocks_[mid].items;
            if (!it.empty() && it.back().t < t) lo = mid + 1;
            else hi = mid;
        }
        return lo < blocks_.size() ? lo : blocks_.size() - 1;
    }

    void split(size_t b) {
        Block nb;
        auto& it = blocks_[b].items;
        const size_t half = it.size() / 2;
        nb.items.assign(it.begin() + static_cast<long>(half), it.end());
        it.resize(half);
        int s = 0;
        for (const auto& x : nb.items) s += x.delta;
        nb.sum = s;
        blocks_[b].sum -= s;
        blocks_.insert(blocks_.begin() + static_cast<long>(b) + 1, std::move(nb));
    }
};

}  // namespace hapcoal
