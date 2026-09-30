#include "KDTree.h"
#include "utils.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>

KDTreeND::KDTreeND(const std::vector<double>& points,
                   const int dim,
                   const std::vector<int>& weights)
    : points_(points), weights_(weights), n_(static_cast<int>(weights.size())), dim_(dim)
{
    if (points_.size() != static_cast<std::size_t>(n_) * dim_) {
        throw std::runtime_error("KDTreeND needs one weight per point");
    }
    indices_.resize(n_);
    std::iota(indices_.begin(), indices_.end(), 0);
    build(0, n_);

    ordered_points_.resize(points_.size());
    ordered_weights_.resize(n_);
    weight_prefix_.assign(n_ + 1, 0);
    for (int p = 0; p < n_; ++p) {
        std::copy(row(indices_[p]), row(indices_[p]) + dim_,
                  ordered_points_.data() + static_cast<std::size_t>(p) * dim_);
        ordered_weights_[p] = weights_[indices_[p]];
        weight_prefix_[p + 1] = weight_prefix_[p] + ordered_weights_[p];
    }

    // Every point is at least half the box diagonal from its farthest corner,
    // so a box with a longer half-diagonal than the radius is never inside it.
    half_diagonal2_.resize(nodes_.size());
    for (std::size_t node = 0; node < nodes_.size(); ++node) {
        const double* minv = bbox_min_.data() + node * dim_;
        const double* maxv = bbox_max_.data() + node * dim_;
        double sum = 0.0;
        for (int d = 0; d < dim_; ++d) {
            const double half = 0.5 * (maxv[d] - minv[d]);
            sum += half * half;
        }
        half_diagonal2_[node] = sum;
    }
}

// The margin covers the rounding in either sum many times over, so a box
// failing this check would fail bbox_inside_radius2 too.
bool KDTreeND::may_be_inside(const int node_idx, const double radius2) const {
    return !(half_diagonal2_[node_idx] > radius2 * (1.0 + 1e-9));
}

std::vector<int> KDTreeND::weighted_neighbor_counts(const double radius2,
                                                    const int threads) const {
    const int workers = normalize_thread_count(threads, std::max(n_, 1));
    std::vector<PairCounts> per_worker(
        workers, PairCounts{std::vector<int>(n_, 0), std::vector<int>(n_ + 1, 0)});

    // One chunk per worker; each worker strides through the positions so that
    // the shrinking amount of work per position stays evenly spread.
    parallel_for(workers, workers, [&](const int lo, const int hi) {
        for (int worker = lo; worker < hi; ++worker) {
            PairCounts& counts = per_worker[worker];
            for (int p = worker; p < n_; p += workers) {
                count_pairs_above(0, p, ordered_row(p), radius2, counts);
            }
        }
    });

    // Every point also counts itself.
    std::vector<int> ordered_counts(ordered_weights_);
    std::vector<int> ranges(n_ + 1, 0);
    for (const PairCounts& counts : per_worker) {
        for (int p = 0; p < n_; ++p) ordered_counts[p] += counts.own[p];
        for (int p = 0; p <= n_; ++p) ranges[p] += counts.ranges[p];
    }
    std::vector<int> counts(n_, 0);
    int running = 0;
    for (int p = 0; p < n_; ++p) {
        running += ranges[p];
        counts[indices_[p]] = ordered_counts[p] + running;
    }
    return counts;
}

int KDTreeND::weighted_count_within(const double* query, const double radius2) const {
    return count_within(0, query, radius2);
}

const double* KDTreeND::row(const int point) const {
    return points_.data() + static_cast<std::size_t>(point) * dim_;
}

const double* KDTreeND::ordered_row(const int ordered_pos) const {
    return ordered_points_.data() + static_cast<std::size_t>(ordered_pos) * dim_;
}

int KDTreeND::build(const int lo, const int hi) {
    if (lo >= hi) return -1;

    const int node_idx = static_cast<int>(nodes_.size());
    nodes_.push_back({lo, hi});
    bbox_min_.resize(nodes_.size() * dim_);
    bbox_max_.resize(nodes_.size() * dim_);
    double* minv = bbox_min_.data() + static_cast<std::size_t>(node_idx) * dim_;
    double* maxv = bbox_max_.data() + static_cast<std::size_t>(node_idx) * dim_;
    std::copy(row(indices_[lo]), row(indices_[lo]) + dim_, minv);
    std::copy(row(indices_[lo]), row(indices_[lo]) + dim_, maxv);
    for (int p = lo; p < hi; ++p) {
        const double* point = row(indices_[p]);
        for (int d = 0; d < dim_; ++d) {
            minv[d] = std::min(minv[d], point[d]);
            maxv[d] = std::max(maxv[d], point[d]);
        }
    }
    if (hi - lo <= kLeafSize) return node_idx;

    // Split the widest dimension at its median.
    int split_dim = 0;
    for (int d = 1; d < dim_; ++d) {
        if (maxv[d] - minv[d] > maxv[split_dim] - minv[split_dim]) split_dim = d;
    }
    const int mid = (lo + hi) / 2;
    std::nth_element(indices_.begin() + lo, indices_.begin() + mid, indices_.begin() + hi,
        [&](const int a, const int b) { return row(a)[split_dim] < row(b)[split_dim]; });

    const int left = build(lo, mid);
    const int right = build(mid, hi);
    nodes_[node_idx].left = left;
    nodes_[node_idx].right = right;
    return node_idx;
}

// The pair test is the running sum of squared differences in dimension order.
// Since the partial sums only grow, testing the whole sum once is the same
// test as stopping at the first partial sum above radius2. The tests here are
// written as "not above", as the stepwise ones were, so a NaN radius still
// counts every pair.
bool KDTreeND::within_radius2(const int ordered_pos,
                              const double* query,
                              const double radius2) const {
    const double* candidate = ordered_row(ordered_pos);
    double distance2 = 0.0;
    for (int d = 0; d < dim_; ++d) {
        const double diff = candidate[d] - query[d];
        distance2 += diff * diff;
    }
    return !(distance2 > radius2);
}

// The box bounds below compare, dimension by dimension, a difference at least
// (or at most) as large as that of any pair they stand for. Rounding is
// monotone, so every rounded term, and so every rounded partial sum, bounds
// the pair's own the same way: a bound beyond radius2 proves every pair
// beyond it, and one within radius2 proves every pair within.

// Lower bound on the squared distance from the query to any point in the box.
bool KDTreeND::bbox_outside_radius2(const int node_idx,
                                    const double* query,
                                    const double radius2) const {
    double out = 0.0;
    const double* minv = bbox_min_.data() + static_cast<std::size_t>(node_idx) * dim_;
    const double* maxv = bbox_max_.data() + static_cast<std::size_t>(node_idx) * dim_;
    for (int d = 0; d < dim_; ++d) {
        // At most one side is positive, and both are exact maxima with zero,
        // so their sum is exactly the gap on that side.
        const double gap = positive_part(minv[d] - query[d]) + positive_part(query[d] - maxv[d]);
        out += gap * gap;
    }
    return out > radius2;
}

// Upper bound on the squared distance from the query to any point in the box.
bool KDTreeND::bbox_inside_radius2(const int node_idx,
                                   const double* query,
                                   const double radius2) const {
    double out = 0.0;
    const double* minv = bbox_min_.data() + static_cast<std::size_t>(node_idx) * dim_;
    const double* maxv = bbox_max_.data() + static_cast<std::size_t>(node_idx) * dim_;
    for (int d = 0; d < dim_; ++d) {
        const double min_diff = query[d] - minv[d];
        const double max_diff = query[d] - maxv[d];
        const double min_term = min_diff * min_diff;
        const double max_term = max_diff * max_diff;
        out += min_term > max_term ? min_term : max_term;
    }
    return !(out > radius2);
}

// Visits only positions above `ordered_pos`, crediting every pair found to both
// its points: directly for a leaf hit, through the difference array when a whole
// subtree lies inside the radius.
void KDTreeND::count_pairs_above(const int node_idx,
                                 const int ordered_pos,
                                 const double* query,
                                 const double radius2,
                                 PairCounts& counts) const {
    if (node_idx < 0) return;
    const Node& node = nodes_[node_idx];
    if (node.hi <= ordered_pos + 1) return;
    if (bbox_outside_radius2(node_idx, query, radius2)) return;

    const int first = std::max(node.lo, ordered_pos + 1);
    if (node.left < 0 && node.right < 0) {
        for (int p = first; p < node.hi; ++p) {
            if (!within_radius2(p, query, radius2)) continue;
            counts.own[ordered_pos] += ordered_weights_[p];
            counts.own[p] += ordered_weights_[ordered_pos];
        }
        return;
    }
    if (may_be_inside(node_idx, radius2) && bbox_inside_radius2(node_idx, query, radius2)) {
        counts.own[ordered_pos] += weight_prefix_[node.hi] - weight_prefix_[first];
        counts.ranges[first] += ordered_weights_[ordered_pos];
        counts.ranges[node.hi] -= ordered_weights_[ordered_pos];
        return;
    }
    count_pairs_above(node.left, ordered_pos, query, radius2, counts);
    count_pairs_above(node.right, ordered_pos, query, radius2, counts);
}

// Whole-subtree shortcut uses the prefix sums, so no per-node totals are needed.
int KDTreeND::count_within(const int node_idx, const double* query, const double radius2) const {
    if (node_idx < 0 || nodes_.empty()) return 0;
    if (bbox_outside_radius2(node_idx, query, radius2)) return 0;
    const Node& node = nodes_[node_idx];
    if (node.left < 0 && node.right < 0) {
        int total = 0;
        for (int p = node.lo; p < node.hi; ++p) {
            if (within_radius2(p, query, radius2)) total += ordered_weights_[p];
        }
        return total;
    }
    if (may_be_inside(node_idx, radius2) && bbox_inside_radius2(node_idx, query, radius2)) {
        return weight_prefix_[node.hi] - weight_prefix_[node.lo];
    }
    return count_within(node.left, query, radius2) +
           count_within(node.right, query, radius2);
}
