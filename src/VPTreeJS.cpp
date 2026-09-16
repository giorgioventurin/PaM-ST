#include "VPTreeJS.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace {

double square(const double x) {
    return x * x;
}

}  // namespace

VPTreeJS::VPTreeJS(const std::vector<double>& points,
                   const std::vector<double>& entropies,
                   const int dim,
                   const std::vector<int>& weights)
    : points_(points), entropies_(entropies), weights_(weights),
      n_(static_cast<int>(weights.size())), dim_(dim)
{
    if (points_.size() != static_cast<std::size_t>(n_) * dim_ ||
        entropies_.size() != weights_.size()) {
        throw std::runtime_error("VPTreeJS needs one entropy and one weight per point");
    }
    indices_.resize(n_);
    std::iota(indices_.begin(), indices_.end(), 0);
    nodes_.reserve(n_);
    build(0, n_);

    node_points_.resize(nodes_.size() * dim_);
    node_entropies_.resize(nodes_.size());
    for (std::size_t node = 0; node < nodes_.size(); ++node) {
        const int point = nodes_[node].point;
        std::copy(row(point), row(point) + dim_, node_points_.data() + node * dim_);
        node_entropies_[node] = entropies_[point];
    }
}

std::vector<int> VPTreeJS::weighted_neighbor_counts(const double radius,
                                                    const int threads) const {
    std::vector<int> counts(n_, 0);
    parallel_for(n_, threads, [&](const int lo, const int hi) {
        for (int i = lo; i < hi; ++i) {
            counts[i] = count(0, i, row(i), entropies_[i], radius * radius, radius);
        }
    });
    return counts;
}

int VPTreeJS::weighted_count_within(const double* query,
                                    const double query_entropy,
                                    const double radius) const {
    return count(0, -1, query, query_entropy, radius * radius, radius);
}

const double* VPTreeJS::row(const int point) const {
    return points_.data() + static_cast<std::size_t>(point) * dim_;
}

double VPTreeJS::distance2(const int a, const int b) const {
    if (a == b) return 0.0;
    return jensen_shannon_divergence_from_entropy(
        row(a), entropies_[a], row(b), entropies_[b], dim_);
}

// Among up to 5 evenly spaced candidates, pick the one whose distances to 24
// evenly spaced samples vary most (ties: larger mean).
int VPTreeJS::select_vantage(const int lo, const int hi) const {
    const int count = hi - lo;
    if (count < 64) return lo;

    const int candidate_count = std::min(5, count);
    const int sample_count = std::min(24, count);
    auto sampled_pos = [lo, count](const int offset, const int total) {
        if (total <= 1) return lo;
        return lo + static_cast<int>(
            (static_cast<long long>(offset) * (count - 1)) / (total - 1));
    };

    int best_pos = lo;
    double best_score = -1.0;
    double best_mean = -1.0;
    for (int c = 0; c < candidate_count; ++c) {
        const int candidate_pos = sampled_pos(c, candidate_count);
        const int candidate = indices_[candidate_pos];

        double sum = 0.0;
        double sum2 = 0.0;
        int used = 0;
        for (int s = 0; s < sample_count; ++s) {
            const int sample_pos = sampled_pos(s, sample_count);
            if (sample_pos == candidate_pos) continue;

            const double d2 = distance2(candidate, indices_[sample_pos]);
            sum += d2;
            sum2 += d2 * d2;
            ++used;
        }
        if (used == 0) continue;

        const double mean = sum / used;
        const double variance = std::max(0.0, sum2 / used - mean * mean);
        if (variance > best_score || (variance == best_score && mean > best_mean)) {
            best_score = variance;
            best_mean = mean;
            best_pos = candidate_pos;
        }
    }
    return best_pos;
}

int VPTreeJS::build(const int lo, const int hi) {
    if (lo >= hi) return -1;

    std::swap(indices_[lo], indices_[select_vantage(lo, hi)]);
    const int node_idx = static_cast<int>(nodes_.size());
    nodes_.push_back({indices_[lo], lo, hi});
    int weight_sum = 0;
    for (int p = lo; p < hi; ++p) weight_sum += weights_[indices_[p]];
    node_weight_sums_.push_back(weight_sum);
    if (hi - lo == 1) return node_idx;

    // Split the remaining points at the median distance to the vantage point.
    const int first = lo + 1;
    const int mid = first + (hi - first) / 2;
    const int vantage = indices_[lo];
    std::vector<std::pair<double, int>> distances;
    distances.reserve(hi - first);
    double max_distance2 = 0.0;
    for (int p = first; p < hi; ++p) {
        const double d2 = distance2(vantage, indices_[p]);
        distances.emplace_back(d2, indices_[p]);
        max_distance2 = std::max(max_distance2, d2);
    }
    const int middle_offset = mid - first;
    std::nth_element(distances.begin(), distances.begin() + middle_offset, distances.end(),
                     [](const auto& a, const auto& b) { return a.first < b.first; });
    for (int offset = 0; offset < static_cast<int>(distances.size()); ++offset) {
        indices_[first + offset] = distances[offset].second;
    }

    nodes_[node_idx].threshold = std::sqrt(distances[middle_offset].first);
    nodes_[node_idx].max_radius = std::sqrt(max_distance2);
    const int left = build(first, mid);
    const int right = build(mid, hi);
    nodes_[node_idx].left = left;
    nodes_[node_idx].right = right;
    return node_idx;
}

// Triangle-inequality pruning: skip a subtree that cannot reach the query, take
// it whole when it is entirely inside, and otherwise test the vantage point and
// descend. eps absorbs rounding in the squared comparisons.
int VPTreeJS::count(const int node_idx,
                    const int point,
                    const double* query,
                    const double query_entropy,
                    const double radius2,
                    const double radius) const {
    if (node_idx < 0) return 0;

    constexpr double eps = 1e-12;
    const Node& node = nodes_[node_idx];
    const double d2 = point == node.point
        ? 0.0
        : jensen_shannon_divergence_from_entropy(
              query, query_entropy,
              node_points_.data() + static_cast<std::size_t>(node_idx) * dim_,
              node_entropies_[node_idx], dim_);

    if (d2 > square(radius + node.max_radius) + eps) return 0;
    if (radius >= node.max_radius && d2 <= square(radius - node.max_radius) + eps) {
        return node_weight_sums_[node_idx];
    }

    int total = d2 <= radius2 + eps ? weights_[node.point] : 0;
    if (node.left >= 0 && d2 <= square(node.threshold + radius) + eps) {
        total += count(node.left, point, query, query_entropy, radius2, radius);
    }
    if (node.right >= 0 &&
        (node.threshold <= radius + eps || d2 + eps >= square(node.threshold - radius))) {
        total += count(node.right, point, query, query_entropy, radius2, radius);
    }
    return total;
}
