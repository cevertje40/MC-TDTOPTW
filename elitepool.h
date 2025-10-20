#pragma once
#include <cstdint>
#include <unordered_set>
#include <vector>
#include <random>
#include "solution.h"

class ElitePool 
{
    public:
    struct Config 
    {
        int   max_size = 8;      // max # of elites kept
        double min_dist = 0.30;   // min Jaccard distance (on directed arcs) for diversity
        double dist_w = 0.5;   // weight of distance in selection
        double qual_w = 0.5;   // weight of quality (normalized score) in selection
        int   cooldown = 20;     // iterations to avoid immediate reselection
    };

    explicit ElitePool(const Config& cfg = Config());

    // Add 's' to the pool if it’s diverse/competitive.
    void consider(const Sol& s, int maxV);

    // Pick index of an elite to restart from (biased to far & good). -1 if none.
    int pick_idx(const Sol& current, int maxV, std::mt19937& eng) const;

    int pick_farthest_idx(const Sol& current, int maxV) const;

    // Mark chosen elite as used (starts cooldown).
    void mark_used(int idx);

    // Decrement cooldowns by 1 (call each iteration).
    void tick();

    static double arc_distance_frac(const std::unordered_set<std::uint64_t>& A,const std::unordered_set<std::uint64_t>& B);

    static void fill_arcs(const Sol& s, int maxV, std::unordered_set<std::uint64_t>& out);

    // Accessors
    int size() const { return static_cast<int>(pool_.size()); }
    const Sol& get(int idx) const;   // throws std::out_of_range if bad idx
    int best_score() const;          // best score or INT_MIN if empty
    void clear();                    // drop all elites

    private:
    struct Elite 
    {
        Sol sol;
        int score = 0;
        std::unordered_set<std::uint64_t> arcs; // directed arc ids (u->v)
        int cooldown = 0;
    };

    static inline std::uint64_t arc_id(int u, int v, int maxV) {
        return (std::uint64_t)u * (std::uint64_t)maxV + (std::uint64_t)v;
    }
    
    Config cfg_;
    std::vector<Elite> pool_;
};
