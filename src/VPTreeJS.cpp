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
    half_terms_.resize(points_.size());
    self_terms_.resize(points_.size());
    roots_.resize(points_.size());
    for (int i = 0; i < n_; ++i) {
        const std::size_t offset = static_cast<std::size_t>(i) * dim_;
        js_row_terms(row(i), dim_, half_terms_.data() + offset, self_terms_.data() + offset);
    }
    for (std::size_t j = 0; j < points_.size(); ++j) roots_[j] = std::sqrt(points_[j]);
    indices_.resize(n_);
    std::iota(indices_.begin(), indices_.end(), 0);
    nodes_.reserve(n_);
    build(0, n_);

    node_points_.resize(nodes_.size() * dim_);
    node_half_terms_.resize(nodes_.size() * dim_);
    node_roots_.resize(nodes_.size() * dim_);
    node_entropies_.resize(nodes_.size());
    for (std::size_t node = 0; node < nodes_.size(); ++node) {
        const int point = nodes_[node].point;
        std::copy(row(point), row(point) + dim_, node_points_.data() + node * dim_);
        std::copy(half_row(point), half_row(point) + dim_, node_half_terms_.data() + node * dim_);
        std::copy(root_row(point), root_row(point) + dim_, node_roots_.data() + node * dim_);
        node_entropies_[node] = entropies_[point];
    }

    // Boxes bottom-up: nodes are stored in preorder, so children follow parents.
    box_lo_ = node_roots_;
    box_hi_ = node_roots_;
    for (std::size_t node = nodes_.size(); node-- > 0;) {
        double* lo = box_lo_.data() + node * dim_;
        double* hi = box_hi_.data() + node * dim_;
        for (const int child : {nodes_[node].left, nodes_[node].right}) {
            if (child < 0) continue;
            const double* child_lo = box_lo_.data() + static_cast<std::size_t>(child) * dim_;
            const double* child_hi = box_hi_.data() + static_cast<std::size_t>(child) * dim_;
            for (int d = 0; d < dim_; ++d) {
                lo[d] = std::min(lo[d], child_lo[d]);
                hi[d] = std::max(hi[d], child_hi[d]);
            }
        }
    }
}

std::vector<int> VPTreeJS::weighted_neighbor_counts(const double radius,
                                                    const int threads) const {
    std::vector<int> counts(n_, 0);
    parallel_for(n_, threads, [&](const int lo, const int hi) {
        for (int i = lo; i < hi; ++i) {
            const Query query{row(i), half_row(i), self_row(i), root_row(i), entropies_[i]};
            counts[i] = count(0, i, query, radius * radius, radius);
        }
    });
    return counts;
}

int VPTreeJS::weighted_count_within(const double* query,
                                    const double query_entropy,
                                    const double radius) const {
    std::vector<double> terms(3 * static_cast<std::size_t>(dim_));
    double* half = terms.data();
    double* self = half + dim_;
    double* roots = self + dim_;
    js_row_terms(query, dim_, half, self);
    for (int d = 0; d < dim_; ++d) roots[d] = std::sqrt(query[d]);
    return count(0, -1, Query{query, half, self, roots, query_entropy}, radius * radius, radius);
}

const double* VPTreeJS::row(const int point) const {
    return points_.data() + static_cast<std::size_t>(point) * dim_;
}

const double* VPTreeJS::half_row(const int point) const {
    return half_terms_.data() + static_cast<std::size_t>(point) * dim_;
}

const double* VPTreeJS::self_row(const int point) const {
    return self_terms_.data() + static_cast<std::size_t>(point) * dim_;
}

const double* VPTreeJS::root_row(const int point) const {
    return roots_.data() + static_cast<std::size_t>(point) * dim_;
}

double VPTreeJS::distance2(const int a, const int b) const {
    if (a == b) return 0.0;
    return jensen_shannon_divergence_from_terms(
        row(a), half_row(a), self_row(a), entropies_[a],
        row(b), half_row(b), entropies_[b], dim_);
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

namespace {

// What a node does with a query: skip its subtree, take the subtree whole, or
// count its vantage point and descend into either half.
enum VisitFlags { kPrune = 1, kTakeAll = 2, kOwn = 4, kLeft = 8, kRight = 16 };

constexpr double kEps = 1e-12;  // absorbs rounding in the squared comparisons

// The limits one node compares the squared distance against, formed once per
// visit exactly as the comparisons below would form them.
struct NodeTests {
    double prune_above;
    double take_all_upto;
    double own_upto;
    double left_upto;
    double right_from;
    bool can_take_all;
    bool has_left;
    bool has_right;
    bool always_right;
};

// Triangle-inequality pruning: skip a subtree that cannot reach the query, take
// it whole when it is entirely inside, and otherwise test the vantage point and
// descend. Every test is monotone in d2.
int visit_flags(const NodeTests& tests, const double d2) {
    if (d2 > tests.prune_above) return kPrune;
    if (tests.can_take_all && d2 <= tests.take_all_upto) return kTakeAll;
    return (d2 <= tests.own_upto ? kOwn : 0) |
           (tests.has_left && d2 <= tests.left_upto ? kLeft : 0) |
           (tests.has_right && (tests.always_right || d2 + kEps >= tests.right_from) ? kRight : 0);
}

}  // namespace

// Since every test is monotone in d2, both ends of an interval holding d2
// leading to the same visit means d2 itself leads to it. The divergence is
// first bounded cheaply through square roots, then more tightly, and computed
// exactly only when the bounds still straddle a test; the visits are those of
// the exact divergence throughout.
int VPTreeJS::count(const int node_idx,
                    const int point,
                    const Query& query,
                    const double radius2,
                    const double radius) const {
    if (node_idx < 0) return 0;

    // A subtree adds to the count only through a vantage point within the
    // radius: its own test needs one, and taking a subtree whole needs
    // d2 <= (radius - max_radius)^2 + eps, which is no larger. When the box of
    // square roots puts every point of the subtree beyond the radius, the
    // subtree adds nothing and is skipped without changing the count.
    const std::size_t offset = static_cast<std::size_t>(node_idx) * dim_;
    const double own_upto = radius2 + kEps;
    {
        const double* lo = box_lo_.data() + offset;
        const double* hi = box_hi_.data() + offset;
        double gap2 = 0.0;
        for (int d = 0; d < dim_; ++d) {
            // At most one of the two differences is positive, since lo <= hi.
            const double below = positive_part(lo[d] - query.roots[d]);
            const double above = positive_part(query.roots[d] - hi[d]);
            gap2 += below * below + above * above;
        }
        if (widen_divergence_bounds(0.34657359027997264 * gap2, 0.0, dim_).lo > own_upto) {
            return 0;
        }
    }

    const Node& node = nodes_[node_idx];
    NodeTests tests;
    tests.prune_above = square(radius + node.max_radius) + kEps;
    tests.can_take_all = radius >= node.max_radius;
    tests.take_all_upto = square(radius - node.max_radius) + kEps;
    tests.own_upto = own_upto;
    tests.has_left = node.left >= 0;
    tests.left_upto = square(node.threshold + radius) + kEps;
    tests.has_right = node.right >= 0;
    tests.always_right = node.threshold <= radius + kEps;
    tests.right_from = square(node.threshold - radius);

    int visit = 0;
    if (point == node.point) {
        visit = visit_flags(tests, 0.0);
    } else {
        const double* vantage = node_points_.data() + offset;
        auto settled = [&](const DivergenceBounds& bounds) {
            visit = visit_flags(tests, bounds.lo);
            return visit == visit_flags(tests, bounds.hi);
        };
        if (!settled(jensen_shannon_bounds_from_roots(
                query.roots, node_roots_.data() + offset, dim_)) &&
            !settled(jensen_shannon_bounds(query.row, vantage, dim_))) {
            visit = visit_flags(tests, jensen_shannon_divergence_from_terms(
                query.row, query.half, query.self, query.entropy,
                vantage, node_half_terms_.data() + offset, node_entropies_[node_idx], dim_));
        }
    }

    if (visit & kPrune) return 0;
    if (visit & kTakeAll) return node_weight_sums_[node_idx];
    int total = (visit & kOwn) ? weights_[node.point] : 0;
    if (visit & kLeft) total += count(node.left, point, query, radius2, radius);
    if (visit & kRight) total += count(node.right, point, query, radius2, radius);
    return total;
}
