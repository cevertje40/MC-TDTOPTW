#include "ElitePool.h"

ElitePool::ElitePool(const Config& cfg) : cfg_(cfg) {}

// Jaccard distance on directed arc sets
double ElitePool::arc_distance_frac(const std::unordered_set<std::uint64_t>& A,const std::unordered_set<std::uint64_t>& B)
{
    if (A.empty() && B.empty()) return 0.0;
    std::size_t inter = 0;
    if (A.size() < B.size()) 
    {
        for (auto x : A) if (B.count(x)) ++inter;
    }
    else 
    {
        for (auto x : B) if (A.count(x)) ++inter;
    }
    const double uni = double(A.size() + B.size() - inter);
    return uni > 0.0 ? (1.0 - double(inter) / uni) : 0.0;
}

void ElitePool::fill_arcs(const Sol& s, int maxV, std::unordered_set<std::uint64_t>& out)
{
    out.clear();
    for (const auto& t : s.tours) {
        const int n = static_cast<int>(t.seq.size());
        for (int k = 0; k + 1 < n; ++k) {
            const int u = t.seq[k]->index;
            const int v = t.seq[k + 1]->index;
            if (u != v) out.insert(arc_id(u, v, maxV));
        }
    }
}

void ElitePool::consider(const Sol& s, int maxV)
{
    Elite cand;
    cand.sol = s;               // deep copy intentionally
    cand.score = s.score;
    fill_arcs(s, maxV, cand.arcs);

    // Reject if too similar to an equal/better elite
    for (const auto& q : pool_) 
    {
        const double d = arc_distance_frac(cand.arcs, q.arcs);
        if (d < cfg_.min_dist&& cand.score <= q.score) return;
    }

    // Insert & trim by quality (and tie-break by #arcs)
    pool_.push_back(std::move(cand));
    std::sort(pool_.begin(), pool_.end(), [](const Elite& a, const Elite& b) {
        if (a.score != b.score) return a.score > b.score;
        return a.arcs.size() > b.arcs.size();
        });
    if (static_cast<int>(pool_.size()) > cfg_.max_size) pool_.pop_back();
}

int ElitePool::pick_idx(const Sol& current, int maxV, std::mt19937& eng) const
{
    if (pool_.empty()) return -1;

    std::unordered_set<std::uint64_t> A;
    fill_arcs(current, maxV, A);

    struct Item { int idx; double w; };
    std::vector<Item> cand;
    cand.reserve(pool_.size());

    int best = pool_.front().score;
    int worst = pool_.back().score;
    double range = std::max(1, best - worst);

    for (int i = 0; i < static_cast<int>(pool_.size()); ++i) {
        if (pool_[i].cooldown > 0) continue;
        double dist = arc_distance_frac(A, pool_[i].arcs);          // [0,1]
        double qual = (pool_[i].score - worst) / range;             // [0,1]
        double w = cfg_.dist_w * dist + cfg_.qual_w * qual;      // linear mix
        if (!(w > 0.0) || !std::isfinite(w)) w = 0.0;
        cand.push_back({ i, w });
    }
    if (cand.empty()) return -1;

    // Build weights for MSVC-friendly discrete_distribution
    std::vector<double> weights; weights.reserve(cand.size());
    for (auto& it : cand) weights.push_back(std::max(1e-12, it.w));

    double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
    if (!(sum > 0.0)) {
        std::uniform_int_distribution<int> ud(0, static_cast<int>(cand.size()) - 1);
        return cand[ud(eng)].idx;
    }

    std::discrete_distribution<int> dd(weights.begin(), weights.end());
    int pos = dd(eng);
    return cand[pos].idx;
}

int ElitePool::pick_farthest_idx(const Sol& current, int maxV) const 
{
    if (pool_.empty()) return -1;
    std::unordered_set<uint64_t> A; fill_arcs(current, maxV, A);
    int best_i = -1; double best_d = -1.0;

    for (int i = 0; i < (int)pool_.size(); ++i) 
    {
        if (pool_[i].cooldown > 0) continue;
        double d = arc_distance_frac(A, pool_[i].arcs);
        if (d > best_d) { best_d = d; best_i = i; }
    }
    return best_i;
}

void ElitePool::mark_used(int idx)
{
    if (idx < 0 || idx >= static_cast<int>(pool_.size())) return;
    pool_[idx].cooldown = cfg_.cooldown;
}

void ElitePool::tick()
{
    for (auto& e : pool_) if (e.cooldown > 0) --e.cooldown;
}

const Sol& ElitePool::get(int idx) const
{
    if (idx < 0 || idx >= static_cast<int>(pool_.size()))
        throw std::out_of_range("ElitePool::get invalid index");
    return pool_[idx].sol;
}

int ElitePool::best_score() const
{
    return pool_.empty() ? std::numeric_limits<int>::min() : pool_.front().score;
}

void ElitePool::clear()
{
    pool_.clear();
}
