#include "utils.h"
#include "KDTree.h"
#include "VPTreeJS.h"
#include "null_model.h"
#include "significance.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <queue>
#include <random>
#include <unordered_map>

namespace {

// A neighbourhood-search bucket. The sample is part of the key, so cells of
// different tissues, which share coordinate ranges, are never neighbours.
struct GridKey {
    int sample = 0;
    std::int64_t x = 0;
    std::int64_t y = 0;
    bool operator==(const GridKey& other) const = default;
};

struct GridKeyHash {
    std::size_t operator()(const GridKey& key) const {
        const auto a = static_cast<std::uint64_t>(key.x) * 0x9E3779B185EBCA87ULL;
        const auto b = static_cast<std::uint64_t>(key.y) * 0xC2B2AE3D27D4EB4FULL;
        const auto s = static_cast<std::uint64_t>(key.sample) * 0x165667B19E3779F9ULL;
        return static_cast<std::size_t>(a ^ (b + (a << 6U) + (a >> 2U)) ^ s);
    }
};

using Neighbors = std::vector<std::vector<int>>;

// Distinct compositions (count vectors equal up to a common factor).
struct CompositionClasses {
    std::vector<int> rows;            // first count vector seen for each class, flattened
    std::vector<int> multiplicities;  // neighbourhoods in each class
    std::vector<int> first_index;     // neighbourhood where each class first appears
};

struct MotifCandidate {
    int count = 0;
    int first_index = 0;
    std::vector<int> pattern;
};

// The candidate neighbourhoods of one analysis: centers[i] owns neighbors[rows[i]].
struct CandidateNeighborhoods {
    const std::vector<int>& centers;
    const Neighbors& neighbors;
    const std::vector<int>& rows;
    const std::vector<int>& observed_histograms;
};

const double* row_of(const std::vector<double>& flat, const int index, const int k) {
    return flat.data() + static_cast<std::size_t>(index) * k;
}

const int* row_of(const std::vector<int>& flat, const int index, const int k) {
    return flat.data() + static_cast<std::size_t>(index) * k;
}

// ---------------------------------------------------------------------------
// Cells, centers and neighbourhoods

std::vector<unsigned char> frozen_label_mask(const Dataset& data,
                                             const std::vector<std::string>& names) {
    std::vector<unsigned char> mask(data.label_names.size(), 0);
    for (const std::string& name : names) {
        const auto it = std::find(data.label_names.begin(), data.label_names.end(), name);
        if (it == data.label_names.end()) {
            throw std::runtime_error("Unknown frozen cell type: " + name);
        }
        unsigned char& frozen = mask[it - data.label_names.begin()];
        if (frozen) throw std::runtime_error("Duplicate frozen cell type: " + name);
        frozen = 1;
    }
    return mask;
}

std::vector<int> non_frozen_cells(const std::vector<int>& labels,
                                  const std::vector<unsigned char>& frozen_labels) {
    std::vector<int> cells;
    for (int i = 0; i < static_cast<int>(labels.size()); ++i) {
        if (!frozen_labels[labels[i]]) cells.push_back(i);
    }
    return cells;
}

// Greedy set cover: repeatedly take the candidate covering the most uncovered
// cells (ties: less overlap, then lower id), re-scoring stale queue entries.
std::vector<int> covering_centers(const Neighbors& neighbors,
                                  const std::vector<int>& candidates,
                                  const int n_cells) {
    struct Entry {
        int uncovered = 0;
        int overlap = 0;
        int center = 0;
        bool operator<(const Entry& rhs) const {
            if (uncovered != rhs.uncovered) return uncovered < rhs.uncovered;
            if (overlap != rhs.overlap) return overlap > rhs.overlap;
            return center > rhs.center;
        }
    };
    std::priority_queue<Entry> queue;
    for (const int center : candidates) {
        queue.push({static_cast<int>(neighbors[center].size()), 0, center});
    }

    std::vector<unsigned char> covered(n_cells, 0);
    std::vector<int> selected;
    while (!queue.empty()) {
        const Entry cached = queue.top();
        queue.pop();
        const int size = static_cast<int>(neighbors[cached.center].size());
        const int uncovered = static_cast<int>(std::count_if(
            neighbors[cached.center].begin(), neighbors[cached.center].end(),
            [&](const int cell) { return !covered[cell]; }));
        if (uncovered != cached.uncovered || size - uncovered != cached.overlap) {
            queue.push({uncovered, size - uncovered, cached.center});
            continue;
        }
        if (uncovered == 0) break;
        selected.push_back(cached.center);
        for (const int cell : neighbors[cached.center]) covered[cell] = 1;
    }
    std::sort(selected.begin(), selected.end());
    return selected;
}

std::vector<int> abundance_order(const Dataset& data, const std::vector<int>& counts) {
    std::vector<int> order(data.label_names.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](const int a, const int b) {
        if (counts[a] != counts[b]) return counts[a] > counts[b];
        return data.label_names[a] < data.label_names[b];
    });
    return order;
}

// neighbors[i] lists the cells of the same sample within `radius` of
// centers[i], using a hash grid with bucket size `radius`.
Neighbors build_neighbors(const Dataset& data,
                          const std::vector<int>& centers,
                          const double radius,
                          const int threads) {
    const std::vector<std::array<double, 2>>& coords = data.coords;
    if (radius < 0.0) throw std::runtime_error("radius must be non-negative");
    const int n_centers = static_cast<int>(centers.size());
    Neighbors neighbors(n_centers);
    if (radius == 0.0) {
        for (int i = 0; i < n_centers; ++i) neighbors[i].push_back(centers[i]);
        return neighbors;
    }

    auto bucket_of = [&](const int cell) {
        const std::array<double, 2>& p = coords[cell];
        return GridKey{data.sample_of(cell),
                       static_cast<std::int64_t>(std::floor(p[0] / radius)),
                       static_cast<std::int64_t>(std::floor(p[1] / radius))};
    };
    std::unordered_map<GridKey, std::vector<int>, GridKeyHash> grid;
    grid.reserve(static_cast<std::size_t>(coords.size() * 1.3));
    for (int i = 0; i < static_cast<int>(coords.size()); ++i) {
        grid[bucket_of(i)].push_back(i);
    }

    const double radius2 = radius * radius;
    parallel_for(n_centers, threads, [&](const int lo, const int hi) {
        for (int i = lo; i < hi; ++i) {
            const auto& center = coords[centers[i]];
            const GridKey c = bucket_of(centers[i]);
            for (std::int64_t dx = -1; dx <= 1; ++dx) {
                for (std::int64_t dy = -1; dy <= 1; ++dy) {
                    const auto it = grid.find(GridKey{c.sample, c.x + dx, c.y + dy});
                    if (it == grid.end()) continue;
                    for (const int j : it->second) {
                        const double xdiff = center[0] - coords[j][0];
                        const double ydiff = center[1] - coords[j][1];
                        if (xdiff * xdiff + ydiff * ydiff <= radius2) neighbors[i].push_back(j);
                    }
                }
            }
        }
    });
    return neighbors;
}

// Row i of the flattened result counts the labels of neighbors[rows[i]].
std::vector<int> label_histograms(const Neighbors& neighbors,
                                  const std::vector<int>& rows,
                                  const std::vector<int>& labels,
                                  const int n_labels,
                                  const int threads) {
    std::vector<int> histograms(rows.size() * n_labels, 0);
    parallel_for(static_cast<int>(rows.size()), threads, [&](const int lo, const int hi) {
        for (int i = lo; i < hi; ++i) {
            int* row = histograms.data() + static_cast<std::size_t>(i) * n_labels;
            for (const int cell : neighbors[rows[i]]) ++row[labels[cell]];
        }
    });
    return histograms;
}

std::vector<int> select_histogram_rows(const std::vector<int>& histograms,
                                       const std::vector<int>& rows,
                                       const int n_labels) {
    std::vector<int> selected;
    selected.reserve(rows.size() * n_labels);
    for (const int row : rows) {
        selected.insert(selected.end(), row_of(histograms, row, n_labels),
                        row_of(histograms, row, n_labels) + n_labels);
    }
    return selected;
}

// ---------------------------------------------------------------------------
// Compositions and motif search

int composition_divisor(const int* row, const int k) {
    int divisor = 0;
    for (int d = 0; d < k && divisor != 1; ++d) divisor = std::gcd(divisor, row[d]);
    return divisor;
}

// Open-addressing index of count vectors up to a common factor: each vector is
// keyed by itself divided by the gcd of its entries. Classes are numbered in
// order of first insertion, so the numbering does not depend on the hash.
class CompositionTable {
public:
    CompositionTable(const int k, const int capacity) : k_(k), reduced_(k) {
        std::size_t slots = 16;
        while (slots < 2 * static_cast<std::size_t>(std::max(capacity, 1))) slots *= 2;
        slots_.assign(slots, -1);
        mask_ = slots - 1;
        keys_.reserve(static_cast<std::size_t>(std::max(capacity, 1)) * k);
    }

    // Class of the row, and whether the row started it.
    std::pair<int, bool> insert(const int* row) {
        std::size_t slot = locate(row);
        if (slots_[slot] >= 0) return {slots_[slot], false};
        if (2 * (static_cast<std::size_t>(size_) + 1) > slots_.size()) {
            grow();
            slot = locate(row);
        }
        slots_[slot] = size_;
        keys_.insert(keys_.end(), reduced_.begin(), reduced_.end());
        return {size_++, true};
    }

    // Class of the row, or -1 when no inserted row has its composition.
    int find(const int* row) { return slots_[locate(row)]; }

private:
    int k_;
    int size_ = 0;
    std::size_t mask_ = 0;
    std::vector<int> slots_;  // class index, or -1 for an empty slot
    std::vector<int> keys_;   // reduced vector of each class, flattened
    std::vector<int> reduced_;

    std::size_t hash(const int* key) const {
        std::uint64_t h = 0x9E3779B97F4A7C15ULL;
        for (int d = 0; d < k_; ++d) {
            h = (h ^ static_cast<std::uint32_t>(key[d])) * 0xFF51AFD7ED558CCDULL;
        }
        return static_cast<std::size_t>(h ^ (h >> 29U));
    }

    // Reduces the row into reduced_ and returns its slot: the one holding its
    // class, or the empty one where the class would go.
    std::size_t locate(const int* row) {
        const int divisor = composition_divisor(row, k_);
        for (int d = 0; d < k_; ++d) reduced_[d] = divisor > 1 ? row[d] / divisor : row[d];
        std::size_t slot = hash(reduced_.data()) & mask_;
        while (slots_[slot] >= 0) {
            const int* key = keys_.data() + static_cast<std::size_t>(slots_[slot]) * k_;
            if (std::equal(key, key + k_, reduced_.begin())) break;
            slot = (slot + 1) & mask_;
        }
        return slot;
    }

    void grow() {
        slots_.assign(slots_.size() * 2, -1);
        mask_ = slots_.size() - 1;
        for (int c = 0; c < size_; ++c) {
            std::size_t slot = hash(keys_.data() + static_cast<std::size_t>(c) * k_) & mask_;
            while (slots_[slot] >= 0) slot = (slot + 1) & mask_;
            slots_[slot] = c;
        }
    }
};

bool same_composition(const int* lhs, const int* rhs, const int k) {
    const int lhs_divisor = composition_divisor(lhs, k);
    const int rhs_divisor = composition_divisor(rhs, k);
    if (lhs_divisor == 0 || rhs_divisor == 0) return lhs_divisor == rhs_divisor;
    for (int d = 0; d < k; ++d) {
        if (lhs[d] / lhs_divisor != rhs[d] / rhs_divisor) return false;
    }
    return true;
}

CompositionClasses composition_classes(const std::vector<int>& histograms,
                                       const int n,
                                       const int k) {
    CompositionTable table(k, n);
    CompositionClasses out;
    for (int i = 0; i < n; ++i) {
        const int* row = row_of(histograms, i, k);
        const auto [index, inserted] = table.insert(row);
        if (inserted) {
            out.rows.insert(out.rows.end(), row, row + k);
            out.multiplicities.push_back(1);
            out.first_index.push_back(i);
        } else {
            ++out.multiplicities[index];
        }
    }
    return out;
}

std::vector<double> normalized_rows(const std::vector<int>& rows, const int n, const int k) {
    std::vector<double> out(static_cast<std::size_t>(n) * k, 0.0);
    for (int i = 0; i < n; ++i) {
        const int* row = row_of(rows, i, k);
        double norm2 = 0.0;
        for (int d = 0; d < k; ++d) norm2 += static_cast<double>(row[d]) * row[d];
        if (norm2 == 0.0) continue;
        const double inv_norm = 1.0 / std::sqrt(norm2);
        double* dst = out.data() + static_cast<std::size_t>(i) * k;
        for (int d = 0; d < k; ++d) dst[d] = row[d] * inv_norm;
    }
    return out;
}

std::vector<double> probability_rows(const std::vector<int>& rows, const int n, const int k) {
    std::vector<double> out(static_cast<std::size_t>(n) * k, 0.0);
    for (int i = 0; i < n; ++i) {
        const int* row = row_of(rows, i, k);
        const int total = std::accumulate(row, row + k, 0);
        if (total == 0) continue;
        const double inv_total = 1.0 / static_cast<double>(total);
        double* dst = out.data() + static_cast<std::size_t>(i) * k;
        for (int d = 0; d < k; ++d) dst[d] = row[d] * inv_total;
    }
    return out;
}

// Ranked by count, then by first appearance.
bool motif_rank_better(const int a, const int b, const std::vector<int>& counts,
                       const CompositionClasses& classes) {
    if (counts[a] != counts[b]) return counts[a] > counts[b];
    return classes.first_index[a] < classes.first_index[b];
}

// Highest-ranked classes whose pairwise distances exceed min_separation.
std::vector<int> select_top_motif_indices(const std::vector<int>& counts,
                                          const CompositionClasses& classes,
                                          const std::vector<double>& points,
                                          const std::vector<double>& entropies,
                                          const int dim,
                                          const DistanceMetric metric,
                                          const int limit,
                                          const double min_separation) {
    std::vector<int> order(counts.size());
    std::iota(order.begin(), order.end(), 0);
    auto better = [&](const int a, const int b) { return motif_rank_better(a, b, counts, classes); };
    if (min_separation <= 0.0) {
        const int top = std::min(limit, static_cast<int>(order.size()));
        std::partial_sort(order.begin(), order.begin() + top, order.end(), better);
        order.resize(top);
        return order;
    }
    std::sort(order.begin(), order.end(), better);

    const double min_separation2 = min_separation * min_separation;
    std::vector<int> selected;
    for (const int candidate : order) {
        const double* candidate_row = row_of(points, candidate, dim);
        const bool separated = std::all_of(selected.begin(), selected.end(), [&](const int accepted) {
            const double* accepted_row = row_of(points, accepted, dim);
            double distance2 = 0.0;
            if (metric == DistanceMetric::Euclidean) {
                for (int d = 0; d < dim; ++d) {
                    const double diff = candidate_row[d] - accepted_row[d];
                    distance2 += diff * diff;
                    if (distance2 > min_separation2) break;
                }
            } else {
                distance2 = jensen_shannon_divergence_from_entropy(
                    candidate_row, entropies[candidate], accepted_row, entropies[accepted], dim);
            }
            return distance2 > min_separation2 + 1e-12;
        });
        if (!separated) continue;
        selected.push_back(candidate);
        if (static_cast<int>(selected.size()) == limit) break;
    }
    return selected;
}

// Every distinct composition of one histogram matrix, with the number of
// neighbourhoods within rho of it. Shared by the frequency search and the
// significance search.
struct ScoredClasses {
    CompositionClasses classes;
    std::vector<int> counts;
    std::vector<double> points;     // query rows; empty when rho == 0
    std::vector<double> entropies;  // JS only
    int size() const { return static_cast<int>(counts.size()); }
};

ScoredClasses score_classes(const std::vector<int>& histograms,
                            const int n,
                            const int n_labels,
                            const double rho,
                            const DistanceMetric metric,
                            const int threads) {
    ScoredClasses out;
    out.classes = composition_classes(histograms, n, n_labels);
    out.counts = out.classes.multiplicities;
    if (rho == 0.0) return out;  // NaN rho continues, as in the original search

    const int m = out.size();
    const bool js = metric == DistanceMetric::JensenShannon;
    out.points = js ? probability_rows(out.classes.rows, m, n_labels)
                    : normalized_rows(out.classes.rows, m, n_labels);
    if (js) {
        out.entropies.resize(m);
        for (int i = 0; i < m; ++i) {
            out.entropies[i] = shannon_entropy(row_of(out.points, i, n_labels), n_labels);
        }
    }
    const double max_distance = js ? std::sqrt(std::log(2.0)) : std::sqrt(2.0);
    if (rho >= max_distance) {
        out.counts.assign(m, n);
    } else if (js) {
        out.counts = VPTreeJS(out.points, out.entropies, n_labels, out.classes.multiplicities)
                         .weighted_neighbor_counts(rho, threads);
    } else {
        out.counts = KDTreeND(out.points, n_labels, out.classes.multiplicities)
                         .weighted_neighbor_counts(rho * rho, threads);
    }
    return out;
}

// The top `max_motifs` compositions by frequency (the original search).
std::vector<MotifCandidate> top_motifs(const std::vector<int>& histograms,
                                       const int n,
                                       const int n_labels,
                                       const double rho,
                                       const DistanceMetric metric,
                                       const int threads,
                                       const int max_motifs,
                                       const double min_separation) {
    if (n == 0 || max_motifs <= 0) return {};
    const ScoredClasses scored = score_classes(histograms, n, n_labels, rho, metric, threads);
    const int limit = std::min(max_motifs, scored.size());
    const std::vector<int> top = select_top_motif_indices(
        scored.counts, scored.classes, scored.points, scored.entropies, n_labels, metric,
        limit, rho != 0.0 ? min_separation : 0.0);

    std::vector<MotifCandidate> motifs;
    for (const int idx : top) {
        const int* row = row_of(scored.classes.rows, idx, n_labels);
        motifs.push_back({scored.counts[idx], scored.classes.first_index[idx],
                          std::vector<int>(row, row + n_labels)});
    }
    return motifs;
}

// Histogram rows prepared once for matching against many patterns. Each row is
// normalised (l2) or turned into proportions with its entropy (js) with the
// same arithmetic as a one-off comparison, so every match is unchanged; only
// the per-row work is no longer repeated for every pattern.
class RowMatcher {
public:
    RowMatcher(const std::vector<int>& histograms, const int n, const int k,
               const double rho, const DistanceMetric metric)
        : histograms_(histograms), n_(n), k_(k), rho_(rho),
          js_(metric == DistanceMetric::JensenShannon) {
        if (rho == 0.0) return;
        if (!js_) {
            points_ = normalized_rows(histograms, n, k);
            return;
        }
        points_ = probability_rows(histograms, n, k);
        entropies_.resize(n);
        half_.resize(points_.size());
        self_.resize(points_.size());
        for (int i = 0; i < n; ++i) {
            const std::size_t offset = static_cast<std::size_t>(i) * k;
            entropies_[i] = shannon_entropy(points_.data() + offset, k);
            js_row_terms(points_.data() + offset, k, half_.data() + offset, self_.data() + offset);
        }
    }

    // Marks the rows lying within rho of the pattern.
    std::vector<unsigned char> match(const int* pattern) const {
        const int k = k_;
        std::vector<unsigned char> out(n_, 0);
        if (rho_ == 0.0) {
            for (int i = 0; i < n_; ++i) {
                out[i] = same_composition(row_of(histograms_, i, k), pattern, k);
            }
            return out;
        }

        const double rho2 = rho_ * rho_;
        std::vector<double> pn(k, 0.0);
        if (js_) {
            const int pattern_total = std::accumulate(pattern, pattern + k, 0);
            if (pattern_total > 0) {
                const double inv = 1.0 / static_cast<double>(pattern_total);
                for (int d = 0; d < k; ++d) pn[d] = pattern[d] * inv;
            }
            const double pattern_entropy = shannon_entropy(pn.data(), k);
            std::vector<double> pattern_terms(2 * static_cast<std::size_t>(k));
            js_row_terms(pn.data(), k, pattern_terms.data(), pattern_terms.data() + k);
            const double limit = rho2 + 1e-12;
            for (int i = 0; i < n_; ++i) {
                // The divergence is computed only when its bounds straddle the limit.
                const std::size_t offset = static_cast<std::size_t>(i) * k;
                const double* row = points_.data() + offset;
                const DivergenceBounds bounds = jensen_shannon_bounds(row, pn.data(), k);
                if (bounds.lo > limit) continue;
                out[i] = bounds.hi <= limit ||
                    jensen_shannon_divergence_from_terms(
                        row, half_.data() + offset, self_.data() + offset, entropies_[i],
                        pn.data(), pattern_terms.data(), pattern_entropy, k) <= limit;
            }
            return out;
        }

        double pnorm2 = 0.0;
        for (int d = 0; d < k; ++d) pnorm2 += static_cast<double>(pattern[d]) * pattern[d];
        if (pnorm2 > 0.0) {
            const double inv = 1.0 / std::sqrt(pnorm2);
            for (int d = 0; d < k; ++d) pn[d] = pattern[d] * inv;
        }
        for (int i = 0; i < n_; ++i) {
            const double* row = points_.data() + static_cast<std::size_t>(i) * k;
            // Partial sums only grow, so stopping once one exceeds rho2 is exact.
            double dist2 = 0.0;
            for (int d = 0; d < k && dist2 <= rho2; ++d) {
                const double diff = row[d] - pn[d];
                dist2 += diff * diff;
            }
            out[i] = dist2 <= rho2;
        }
        return out;
    }

private:
    const std::vector<int>& histograms_;
    int n_;
    int k_;
    double rho_;
    bool js_;
    std::vector<double> points_;     // normalised or probability rows
    std::vector<double> entropies_;  // js only
    std::vector<double> half_;       // js only, see js_row_terms
    std::vector<double> self_;
};

// Scatters a row mask onto cell ids, for the match files.
std::vector<unsigned char> motif_matches_for_centers(const std::vector<int>& histograms,
                                                     const int k,
                                                     const std::vector<int>& pattern,
                                                     const double rho,
                                                     const DistanceMetric metric,
                                                     const std::vector<int>& centers,
                                                     const int n_cells) {
    const std::vector<unsigned char> rows =
        RowMatcher(histograms, static_cast<int>(centers.size()), k, rho, metric)
            .match(pattern.data());
    std::vector<unsigned char> out(n_cells, 0);
    for (int i = 0; i < static_cast<int>(centers.size()); ++i) out[centers[i]] = rows[i];
    return out;
}

// ---------------------------------------------------------------------------
// Significance search: a fixed candidate family, each member with its own null

struct CandidateFamily {
    std::vector<int> patterns;      // size x n_labels raw count vectors
    std::vector<double> points;     // size x n_labels query rows
    std::vector<double> entropies;  // JS only
    std::vector<int> observed;      // neighbourhoods within rho of each candidate
    std::vector<int> first_index;
    int size = 0;
};

// Counts every candidate in one (permuted) tissue.
void count_family(const std::vector<int>& histograms, const int n, const int n_labels,
                  const CandidateFamily& family, const AnalysisConfig& config,
                  const int threads, int* out) {
    const double rho = config.rho;
    const bool js = config.metric == DistanceMetric::JensenShannon;

    if (rho == 0.0) {
        // Exact proportional match: look each candidate up among the classes.
        CompositionTable table(n_labels, n);
        std::vector<int> multiplicities;
        for (int i = 0; i < n; ++i) {
            const auto [index, inserted] = table.insert(row_of(histograms, i, n_labels));
            if (inserted) multiplicities.push_back(1);
            else ++multiplicities[index];
        }
        for (int c = 0; c < family.size; ++c) {
            const int index = table.find(row_of(family.patterns, c, n_labels));
            out[c] = index < 0 ? 0 : multiplicities[index];
        }
        return;
    }

    const double max_distance = js ? std::sqrt(std::log(2.0)) : std::sqrt(2.0);
    if (rho >= max_distance) {
        std::fill(out, out + family.size, n);
        return;
    }
    const CompositionClasses classes = composition_classes(histograms, n, n_labels);
    const int m = static_cast<int>(classes.multiplicities.size());
    const std::vector<double> points = js ? probability_rows(classes.rows, m, n_labels)
                                          : normalized_rows(classes.rows, m, n_labels);
    if (js) {
        std::vector<double> entropies(m);
        for (int i = 0; i < m; ++i) {
            entropies[i] = shannon_entropy(row_of(points, i, n_labels), n_labels);
        }
        const VPTreeJS tree(points, entropies, n_labels, classes.multiplicities);
        parallel_for(family.size, threads, [&](const int lo, const int hi) {
            for (int c = lo; c < hi; ++c) {
                out[c] = tree.weighted_count_within(
                    row_of(family.points, c, n_labels), family.entropies[c], rho);
            }
        });
        return;
    }
    const KDTreeND tree(points, n_labels, classes.multiplicities);
    parallel_for(family.size, threads, [&](const int lo, const int hi) {
        for (int c = lo; c < hi; ++c) {
            out[c] = tree.weighted_count_within(row_of(family.points, c, n_labels), rho * rho);
        }
    });
}

// One sample's testing centres: positions into the candidate arrays, the
// neighbourhood rows they own, and their observed histograms.
// The shape of the circles that match a motif: their mean size and their mean
// composition as proportions. Every circle counts once, whatever its size, so
// one dense neighbourhood cannot outvote a dozen sparse ones.
struct MatchedShape {
    double mean_cells = 0.0;
    std::vector<double> mean_composition;
};

MatchedShape matched_shape(const std::vector<int>& histograms,
                           const std::vector<unsigned char>& matched,
                           const int n_rows,
                           const int n_labels) {
    MatchedShape shape;
    shape.mean_composition.assign(n_labels, 0.0);
    int hits = 0;
    for (int i = 0; i < n_rows; ++i) {
        if (!matched[i]) continue;
        const int* row = row_of(histograms, i, n_labels);
        const int total = std::accumulate(row, row + n_labels, 0);
        if (total <= 0) continue;
        ++hits;
        shape.mean_cells += total;
        for (int d = 0; d < n_labels; ++d) {
            shape.mean_composition[d] += static_cast<double>(row[d]) / total;
        }
    }
    if (hits > 0) {
        shape.mean_cells /= hits;
        for (double& value : shape.mean_composition) value /= hits;
    }
    return shape;
}

struct SampleTestSet {
    std::vector<int> positions;
    std::vector<int> rows;
    std::vector<int> histograms;
};

// The null draws. `pooled` is candidate-major, pooled[c * B + b], summed over
// samples; it is what the family-wise test uses. With several samples, running
// per-sample totals are kept too, so each sample's own null mean, spread and
// p-value can be reported without storing a matrix per sample.
struct NullDraws {
    std::vector<int> pooled;
    std::vector<std::vector<long long>> sum;     // [sample][candidate]
    std::vector<std::vector<long long>> sum_sq;  // [sample][candidate]
    std::vector<std::vector<int>> at_least;      // draws >= that sample's observed count
};

// ---------------------------------------------------------------------------
// Permutations

// How the permutations share the threads: `workers` permutations are counted at
// once, each with `threads` threads of its own. Whole permutations parallelise
// best, but each one in flight holds its own histograms, classes and search
// tree, so with large inputs fewer run at once and each gets more threads.
struct PermutationWorkers {
    int workers = 1;
    int threads = 1;
};

// Rough working memory of counting one permutation over `rows` neighbourhoods:
// histograms, composition classes, query rows and search tree, plus its labels.
std::size_t permutation_bytes(const std::size_t rows, const int n_labels, const std::size_t cells) {
    return rows * (64 * static_cast<std::size_t>(n_labels) + 96) + cells * sizeof(int);
}

PermutationWorkers permutation_workers(const AnalysisConfig& config,
                                       const std::size_t bytes_per_permutation) {
    constexpr std::size_t kMemoryBudget = std::size_t{2} << 30;  // for permutations in flight
    const int total = normalize_thread_count(config.threads, std::numeric_limits<int>::max());
    const std::size_t by_memory =
        std::max<std::size_t>(1, kMemoryBudget / std::max<std::size_t>(1, bytes_per_permutation));
    PermutationWorkers pool;
    pool.workers = static_cast<int>(std::min<std::size_t>(
        {static_cast<std::size_t>(total), by_memory,
         static_cast<std::size_t>(std::max(1, config.permutations))}));
    pool.threads = std::max(1, total / pool.workers);
    return pool;
}

// Draws config.permutations label shuffles from one generator seeded with
// config.seed and calls fn(b, worker, labels, threads) for each. The shuffles
// are drawn under a lock in permutation order, each continuing from the last,
// so permutation b sees exactly the labels of a sequential loop; only the
// counting runs concurrently. fn may write only to what belongs to b or to
// worker (in [0, pool.workers)).
template <typename Fn>
void for_each_permutation(const std::vector<int>& labels,
                          const LabelShufflePlan& plan,
                          const AnalysisConfig& config,
                          const PermutationWorkers& pool,
                          Fn fn) {
    const int permutations = config.permutations;
    std::vector<int> shuffled = labels;
    std::mt19937 rng(config.seed);
    std::mutex mutex;
    int next = 0;
    int done = 0;
    std::exception_ptr error;
    parallel_for_each(pool.workers, pool.workers, [&](const int worker, int) {
        std::vector<int> own;
        while (true) {
            int b = 0;
            {
                const std::lock_guard<std::mutex> lock(mutex);
                if (next == permutations || error) return;
                b = next++;
                shuffle_labels_within_groups(plan, shuffled, rng);
                own = shuffled;
            }
            try {
                fn(b, worker, own, pool.threads);
            } catch (...) {
                const std::lock_guard<std::mutex> lock(mutex);
                if (!error) error = std::current_exception();
                return;
            }
            const std::lock_guard<std::mutex> lock(mutex);
            ++done;
            if (done % 50 == 0 || done == permutations) {
                std::cerr << "  permutation " << done << '/' << permutations << '\n';
            }
        }
    });
    if (error) std::rethrow_exception(error);
}

NullDraws null_family_counts(const Dataset& data,
                             const AnalysisConfig& config,
                             const LabelShufflePlan& shuffle_plan,
                             const Neighbors& neighbors,
                             const std::vector<SampleTestSet>& sets,
                             const CandidateFamily& family,
                             const std::vector<std::vector<int>>& observed_by_sample) {
    const int n_labels = static_cast<int>(data.label_names.size());
    const int n_samples = static_cast<int>(sets.size());
    const int permutations = config.permutations;
    const bool per_sample = n_samples > 1;
    NullDraws draws;
    draws.pooled.assign(static_cast<std::size_t>(family.size) * permutations, 0);

    std::size_t largest_set = 0;
    for (const SampleTestSet& set : sets) largest_set = std::max(largest_set, set.rows.size());
    const PermutationWorkers pool = permutation_workers(
        config, permutation_bytes(largest_set, n_labels, data.labels.size()) +
                    static_cast<std::size_t>(family.size) * 3 * sizeof(double));

    // Per-sample totals are kept per worker and added up at the end. They are
    // integers, so the order in which the draws arrive does not matter.
    std::vector<NullDraws> partial(per_sample ? pool.workers : 0);
    for (NullDraws& part : partial) {
        part.sum.assign(n_samples, std::vector<long long>(family.size, 0));
        part.sum_sq.assign(n_samples, std::vector<long long>(family.size, 0));
        part.at_least.assign(n_samples, std::vector<int>(family.size, 0));
    }
    for_each_permutation(data.labels, shuffle_plan, config, pool,
                         [&](const int b, const int worker, const std::vector<int>& labels,
                             const int threads) {
        // One shuffle moves every sample at once; each sample is then counted
        // on its own, and the pooled count is the sum.
        std::vector<int> draw(family.size, 0);
        for (int sample = 0; sample < n_samples; ++sample) {
            const SampleTestSet& set = sets[sample];
            if (set.rows.empty()) continue;
            count_family(label_histograms(neighbors, set.rows, labels, n_labels, threads),
                         static_cast<int>(set.rows.size()), n_labels, family, config, threads,
                         draw.data());
            for (int c = 0; c < family.size; ++c) {
                draws.pooled[static_cast<std::size_t>(c) * permutations + b] += draw[c];
                if (!per_sample) continue;
                NullDraws& part = partial[worker];
                part.sum[sample][c] += draw[c];
                part.sum_sq[sample][c] += static_cast<long long>(draw[c]) * draw[c];
                if (draw[c] >= observed_by_sample[sample][c]) ++part.at_least[sample][c];
            }
        }
    });
    if (per_sample) {
        draws.sum = partial[0].sum;
        draws.sum_sq = partial[0].sum_sq;
        draws.at_least = partial[0].at_least;
        for (int worker = 1; worker < pool.workers; ++worker) {
            for (int sample = 0; sample < n_samples; ++sample) {
                for (int c = 0; c < family.size; ++c) {
                    draws.sum[sample][c] += partial[worker].sum[sample][c];
                    draws.sum_sq[sample][c] += partial[worker].sum_sq[sample][c];
                    draws.at_least[sample][c] += partial[worker].at_least[sample][c];
                }
            }
        }
    }
    return draws;
}

// ---------------------------------------------------------------------------
// Permutation test

MotifTest summarize_motif_test(const MotifCandidate& candidate,
                               std::vector<int> null_counts,
                               const int rank,
                               const double alpha) {
    MotifTest test;
    test.rank = rank;
    test.observed = candidate.count;
    test.pattern = candidate.pattern;
    test.null_counts = std::move(null_counts);

    int ge = 0;
    double sum = 0.0;
    for (const int v : test.null_counts) {
        if (v >= test.observed) ++ge;
        sum += v;
    }
    const int n_permutations = static_cast<int>(test.null_counts.size());
    test.p_value = (1.0 + ge) / (n_permutations + 1.0);
    test.null_mean = sum / n_permutations;

    double var = 0.0;
    for (const int v : test.null_counts) {
        const double diff = v - test.null_mean;
        var += diff * diff;
    }
    test.null_sd = std::sqrt(var / n_permutations);
    if (test.null_sd > 0.0) {
        test.z_score = (test.observed - test.null_mean) / test.null_sd;
    } else if (test.observed != test.null_mean) {
        const double inf = std::numeric_limits<double>::infinity();
        test.z_score = test.observed > test.null_mean ? inf : -inf;
    }
    test.significant = test.p_value <= alpha;
    return test;
}

void record_neighborhood_coverage(Result& result,
                                  const CandidateNeighborhoods& candidates,
                                  const int n_cells) {
    std::vector<int> memberships(n_cells, 0);
    for (const int row : candidates.rows) {
        if (row < 0 || row >= static_cast<int>(candidates.neighbors.size())) {
            throw std::runtime_error("Internal error: neighborhood row is out of range");
        }
        result.neighborhood_membership_count += candidates.neighbors[row].size();
        for (const int cell : candidates.neighbors[row]) {
            const int count = ++memberships[cell];
            if (count == 1) ++result.covered_cell_count;
            else if (count == 2) ++result.overlapped_cell_count;
        }
    }
    result.candidate_center_ids = candidates.centers;
}

void validate_analysis_config(const AnalysisConfig& config) {
    validate_null_model_config(config);
    if (config.rho < 0.0) throw std::runtime_error("rho must be non-negative");
    if (config.permutations < 1) throw std::runtime_error("permutations must be at least 1");
    if (config.alpha <= 0.0 || config.alpha >= 1.0) {
        throw std::runtime_error("alpha must be greater than 0 and less than 1");
    }
    if (config.max_motifs < 1) throw std::runtime_error("max-motifs must be at least 1");
    if (config.error_control != ErrorControl::FewRSFDR) return;

    if (config.fdr_failure_probability <= 0.0 ||
        config.fdr_failure_probability >= config.alpha) {
        throw std::runtime_error("FewRS FDR failure probability must be between 0 and alpha");
    }
    if (config.fdr_min_discoveries < 1 || config.fdr_min_discoveries > config.max_motifs) {
        throw std::runtime_error("FewRS FDR minimum discoveries must be between 1 and max-motifs");
    }
    if (config.fdr_order < 1 || config.fdr_order > config.fdr_min_discoveries) {
        throw std::runtime_error("FewRS FDR null order is out of range");
    }
    if (config.permutations < fewrs_resample_count(config.fdr_failure_probability)) {
        throw std::runtime_error("FewRS FDR has too few permutations for its failure probability");
    }
}

// Re-runs the motif search on every label permutation. Returns the null count
// of each observed rank per permutation; fills the FewRS order statistic.
std::vector<std::vector<int>> permutation_null_counts(Result& result,
                                                      const Dataset& data,
                                                      const AnalysisConfig& config,
                                                      const LabelShufflePlan& shuffle_plan,
                                                      const CandidateNeighborhoods& candidates,
                                                      const int n_tests,
                                                      const double motif_separation) {
    const bool fewrs = config.error_control == ErrorControl::FewRSFDR;
    const int n_labels = static_cast<int>(data.label_names.size());
    const int search_limit = fewrs ? std::max(n_tests, config.fdr_order) : n_tests;
    std::vector<std::vector<int>> null_by_rank(n_tests, std::vector<int>(config.permutations, 0));
    if (fewrs) result.fdr_null_order_statistics.assign(config.permutations, 0);

    const PermutationWorkers pool = permutation_workers(
        config, permutation_bytes(candidates.rows.size(), n_labels, data.labels.size()));
    for_each_permutation(data.labels, shuffle_plan, config, pool,
                         [&](const int b, int, const std::vector<int>& labels, const int threads) {
        const std::vector<MotifCandidate> null_motifs = top_motifs(
            label_histograms(candidates.neighbors, candidates.rows, labels, n_labels, threads),
            static_cast<int>(candidates.centers.size()), n_labels, config.rho, config.metric,
            threads, search_limit, motif_separation);
        const int found = static_cast<int>(null_motifs.size());
        for (int rank = 0; rank < std::min(n_tests, found); ++rank) {
            null_by_rank[rank][b] = null_motifs[rank].count;
        }
        if (fewrs && config.fdr_order <= found) {
            result.fdr_null_order_statistics[b] = null_motifs[config.fdr_order - 1].count;
        }
    });
    return null_by_rank;
}

// FewRS-FDR: motifs above the maximal null order statistic are discoveries when
// there are enough of them for the FDP bound to meet alpha.
void certify_fewrs(Result& result, const AnalysisConfig& config,
                   const CandidateNeighborhoods& candidates, const int n_labels, const int n_cells) {
    result.fdr_threshold = *std::max_element(result.fdr_null_order_statistics.begin(),
                                             result.fdr_null_order_statistics.end());
    result.n_fdr_candidates = static_cast<int>(std::count_if(
        result.motif_tests.begin(), result.motif_tests.end(),
        [&](const MotifTest& test) { return test.observed > result.fdr_threshold; }));
    if (result.n_fdr_candidates > 0) {
        result.fdr_bound = config.fdr_failure_probability +
            static_cast<double>(config.fdr_order - 1) / result.n_fdr_candidates;
    }
    result.fdr_certified = result.n_fdr_candidates >= config.fdr_min_discoveries &&
        result.fdr_bound <= config.alpha + 1e-12;

    for (int rank = 0; rank < static_cast<int>(result.motif_tests.size()); ++rank) {
        MotifTest& test = result.motif_tests[rank];
        test.significant = result.fdr_certified && test.observed > result.fdr_threshold;
        if (test.significant) ++result.n_significant_motifs;
        if (test.significant || rank == 0) {
            test.match_mask = motif_matches_for_centers(
                candidates.observed_histograms, n_labels, test.pattern, config.rho,
                config.metric, candidates.centers, n_cells);
        }
    }
    result.tested_motifs = static_cast<int>(result.motif_tests.size());
    result.selected_rank = 1;
}

// Diagnostics shared by both searches: permutation population, tile counts and
// neighbourhood coverage.
Result start_result(const Dataset& data,
                    const AnalysisConfig& config,
                    const std::vector<int>& permutable_cells,
                    const LabelShufflePlan& shuffle_plan,
                    const CandidateNeighborhoods& candidates) {
    if (candidates.rows.size() != candidates.centers.size()) {
        throw std::runtime_error("Internal error: center and neighbor-row counts differ");
    }
    Result result;
    result.permutable_cell_count = static_cast<int>(permutable_cells.size());
    result.null_block_count = shuffle_plan.block_count;
    result.null_shuffle_group_count = static_cast<int>(shuffle_plan.groups.size());
    result.null_mixed_group_count = shuffle_plan.mixed_group_count;
    result.null_exchangeable_cell_count = shuffle_plan.exchangeable_cell_count;
    result.candidate_center_count = static_cast<int>(candidates.centers.size());
    if (config.null_model == NullModel::Block) {
        std::cerr << "Block null: " << result.null_block_count << " occupied tiles, "
                  << result.null_shuffle_group_count << " shuffle groups, "
                  << result.null_mixed_group_count << " mixed groups, "
                  << result.null_exchangeable_cell_count << " exchangeable cells"
                  << (result.null_exchangeable_cell_count == 0 ? " (no labels can change)" : "")
                  << '\n';
    }
    record_neighborhood_coverage(result, candidates, static_cast<int>(data.labels.size()));
    return result;
}

Result run_pattern_test_prepared(const Dataset& data,
                                 const AnalysisConfig& config,
                                 const std::vector<int>& permutable_cells,
                                 const LabelShufflePlan& shuffle_plan,
                                 const CandidateNeighborhoods& candidates,
                                 const bool select_rank_one) {
    const int n_cells = static_cast<int>(data.labels.size());
    const int n_labels = static_cast<int>(data.label_names.size());
    Result result = start_result(data, config, permutable_cells, shuffle_plan, candidates);

    const bool fewrs = config.error_control == ErrorControl::FewRSFDR;
    const double motif_separation = fewrs ? 2.0 * config.rho : 0.0;
    const std::vector<MotifCandidate> observed = top_motifs(
        candidates.observed_histograms, result.candidate_center_count, n_labels, config.rho,
        config.metric, config.threads, config.max_motifs, motif_separation);
    if (observed.empty()) {
        result.motif_tests.push_back(summarize_motif_test(
            {0, 0, std::vector<int>(n_labels, 0)},
            std::vector<int>(config.permutations, 0), 1, config.alpha));
        return result;
    }

    const int n_tests = static_cast<int>(observed.size());
    std::vector<std::vector<int>> null_by_rank = permutation_null_counts(
        result, data, config, shuffle_plan, candidates, n_tests, motif_separation);
    for (int rank = 0; rank < n_tests; ++rank) {
        result.motif_tests.push_back(summarize_motif_test(
            observed[rank], std::move(null_by_rank[rank]), rank + 1, config.alpha));
    }

    if (fewrs) {
        certify_fewrs(result, config, candidates, n_labels, n_cells);
        return result;
    }

    // Pointwise: report rank 1, or (standalone runs) the first significant rank.
    int selected = 0;
    if (!select_rank_one) {
        const auto first_significant = std::find_if(
            result.motif_tests.begin(), result.motif_tests.end(),
            [](const MotifTest& test) { return test.significant; });
        selected = first_significant == result.motif_tests.end()
            ? n_tests - 1
            : static_cast<int>(first_significant - result.motif_tests.begin());
        result.motif_tests.resize(selected + 1);
    }
    result.tested_motifs = select_rank_one ? n_tests : selected + 1;
    result.selected_rank = selected + 1;
    MotifTest& motif = result.motif_tests[selected];
    motif.match_mask = motif_matches_for_centers(
        candidates.observed_histograms, n_labels, motif.pattern, config.rho, config.metric,
        candidates.centers, n_cells);
    result.n_significant_motifs = motif.significant ? 1 : 0;
    return result;
}

// Positions into the candidate arrays, split by the parity of the tile holding
// each center. A center is kept only when its whole neighbourhood fits inside
// its own tile, so no testing neighbourhood shares cells with a selecting one.
struct CheckerboardSplit {
    std::vector<int> selection;
    std::vector<int> test;
    int discarded = 0;
};

CheckerboardSplit checkerboard_split(const Dataset& data,
                                     const AnalysisConfig& config,
                                     const CandidateNeighborhoods& candidates) {
    const double size = config.split_size;
    const bool block = config.null_model == NullModel::Block;
    const double origin_x = block ? config.block_origin_x : 0.0;
    const double origin_y = block ? config.block_origin_y : 0.0;
    CheckerboardSplit split;
    for (int i = 0; i < static_cast<int>(candidates.centers.size()); ++i) {
        const std::array<double, 2>& point = data.coords[candidates.centers[i]];
        const std::array<std::int64_t, 2> tile = grid_tile(point, size, origin_x, origin_y);
        const double lo_x = origin_x + static_cast<double>(tile[0]) * size;
        const double lo_y = origin_y + static_cast<double>(tile[1]) * size;
        const bool fits =
            point[0] - lo_x >= config.radius && lo_x + size - point[0] >= config.radius &&
            point[1] - lo_y >= config.radius && lo_y + size - point[1] >= config.radius;
        if (!fits) {
            ++split.discarded;
            continue;
        }
        const bool selects = ((tile[0] + tile[1]) % 2 + 2) % 2 == 0;
        (selects ? split.selection : split.test).push_back(i);
    }
    return split;
}

std::vector<int> subset_rows(const std::vector<int>& histograms,
                             const std::vector<int>& positions,
                             const int n_labels) {
    std::vector<int> out;
    out.reserve(positions.size() * n_labels);
    for (const int position : positions) {
        const int* row = row_of(histograms, position, n_labels);
        out.insert(out.end(), row, row + n_labels);
    }
    return out;
}

// Occurrences sharing no cell, taken greedily in center order: a dense clump of
// overlapping neighbourhoods counts once, not once per neighbourhood.
// Stops early once the support reaches `limit`, since it can only grow.
int disjoint_support(const std::vector<unsigned char>& matched,
                     const std::vector<int>& positions,
                     const CandidateNeighborhoods& candidates,
                     const int n_cells,
                     const int limit = std::numeric_limits<int>::max()) {
    std::vector<unsigned char> used(n_cells, 0);
    int support = 0;
    for (int i = 0; i < static_cast<int>(positions.size()) && support < limit; ++i) {
        if (!matched[i]) continue;
        const std::vector<int>& cells = candidates.neighbors[candidates.rows[positions[i]]];
        if (std::any_of(cells.begin(), cells.end(),
                        [&](const int cell) { return used[cell] != 0; })) {
            continue;
        }
        for (const int cell : cells) used[cell] = 1;
        ++support;
    }
    return support;
}

// Whether the disjoint support reaches min_support. It never exceeds the
// number of matching rows, which settles most candidates without the scan.
bool has_disjoint_support(const std::vector<unsigned char>& matched,
                          const std::vector<int>& positions,
                          const CandidateNeighborhoods& candidates,
                          const int n_cells,
                          const int min_support) {
    if (std::count(matched.begin(), matched.end(), 1) < min_support) return false;
    return disjoint_support(matched, positions, candidates, n_cells, min_support) >= min_support;
}

// True when a cell sits in a tile that the checkerboard assigns to testing.
bool in_test_tile(const std::array<double, 2>& point, const AnalysisConfig& config) {
    const bool block = config.null_model == NullModel::Block;
    const std::array<std::int64_t, 2> tile = grid_tile(
        point, config.split_size, block ? config.block_origin_x : 0.0,
        block ? config.block_origin_y : 0.0);
    return ((tile[0] + tile[1]) % 2 + 2) % 2 != 0;
}

// The null has to condition on whatever chose the family. Only labels in the
// testing tiles are shuffled; the selection half keeps its observed labels, so
// the permutation distribution matches the conditional distribution the
// selected family was drawn from.
LabelShufflePlan restrict_to_test_tiles(const LabelShufflePlan& plan,
                                        const Dataset& data,
                                        const AnalysisConfig& config) {
    LabelShufflePlan restricted;
    restricted.block_count = plan.block_count;
    for (const std::vector<int>& group : plan.groups) {
        std::vector<int> kept;
        for (const int cell : group) {
            if (in_test_tile(data.coords[cell], config)) kept.push_back(cell);
        }
        if (!kept.empty()) restricted.groups.push_back(std::move(kept));
    }
    for (const std::vector<int>& group : restricted.groups) {
        const bool mixed = std::any_of(group.begin() + 1, group.end(),
            [&](const int cell) { return data.labels[cell] != data.labels[group.front()]; });
        if (mixed) {
            ++restricted.mixed_group_count;
            restricted.exchangeable_cell_count += static_cast<int>(group.size());
        }
    }
    return restricted;
}

// Reports the most significant candidates instead of the most frequent one.
//
// The family is chosen on one half of a checkerboard of tiles and everything
// inferential - observed counts, nulls, p-values - is computed on the other
// half. Selecting and testing on the same neighbourhoods would guarantee small
// p-values: a composition would be picked *because* it was frequent, and its
// null, which never selects it, sits far below.
Result run_significance_test(const Dataset& data,
                             const AnalysisConfig& config,
                             const std::vector<int>& permutable_cells,
                             const LabelShufflePlan& shuffle_plan,
                             const CandidateNeighborhoods& candidates) {
    const int n_cells = static_cast<int>(data.labels.size());
    const int n_labels = static_cast<int>(data.label_names.size());
    Result result = start_result(data, config, permutable_cells, shuffle_plan, candidates);

    const CheckerboardSplit split = checkerboard_split(data, config, candidates);
    result.selection_center_count = static_cast<int>(split.selection.size());
    result.test_center_count = static_cast<int>(split.test.size());
    result.split_discarded_count = split.discarded;
    std::cerr << "Checkerboard split at " << config.split_size << " units: "
              << result.selection_center_count << " centers select, "
              << result.test_center_count << " test, " << split.discarded
              << " dropped at tile edges\n";
    if (split.selection.empty() || split.test.empty()) {
        throw std::runtime_error("Checkerboard split left one half empty; raise --split-size");
    }

    // Selection half: compositions that recur in enough separate places.
    const std::vector<int> selection_histograms =
        subset_rows(candidates.observed_histograms, split.selection, n_labels);
    const int n_selection = static_cast<int>(split.selection.size());
    const ScoredClasses selection = score_classes(
        selection_histograms, n_selection, n_labels, config.rho, config.metric, config.threads);
    // Each class is judged on its own, so the classes are judged in parallel;
    // the family then keeps them in class order.
    const RowMatcher selection_rows(
        selection_histograms, n_selection, n_labels, config.rho, config.metric);
    std::vector<unsigned char> in_family(selection.size(), 0);
    parallel_for_each(selection.size(), config.threads, [&](const int i, int) {
        const int* pattern = row_of(selection.classes.rows, i, n_labels);
        if (!passes_family_filters(pattern, selection.counts[i], n_labels, config)) return;
        const std::vector<unsigned char> matched = selection_rows.match(pattern);
        in_family[i] = has_disjoint_support(matched, split.selection, candidates, n_cells,
                                            config.min_support);
    });
    CandidateFamily family;
    for (int i = 0; i < selection.size(); ++i) {
        if (!in_family[i]) continue;
        const int* pattern = row_of(selection.classes.rows, i, n_labels);
        family.patterns.insert(family.patterns.end(), pattern, pattern + n_labels);
        family.first_index.push_back(selection.classes.first_index[i]);
        ++family.size;
    }
    result.candidate_family_size = family.size;
    std::cerr << "Candidate family: " << family.size << " of " << selection.size()
              << " compositions recur in >= " << config.min_support
              << " disjoint places in the selection half";
    if (config.min_types > 0) {
        std::cerr << " with >= " << config.min_types << " types of >= "
                  << config.min_type_cells << " cells";
    }
    std::cerr << "\n";
    // The family minimum hits its floor in up to family/(B+1) of the draws, so
    // that ratio bounds the smallest adjusted p-value a run of independent
    // candidates can produce. Correlated candidates, which near-duplicate
    // compositions are, reach lower.
    if (config.permutations + 1 < 20 * family.size) {
        std::cerr << "  warning: " << config.permutations << " permutations for "
                  << family.size << " candidates; the smallest attainable adjusted "
                     "p-value is about "
                  << static_cast<double>(family.size) / (config.permutations + 1)
                  << " for independent candidates. Raise --permutations or "
                     "--min-support if nothing reaches significance.\n";
    }
    if (family.size == 0) {
        result.motif_tests.push_back(summarize_motif_test(
            {0, 0, std::vector<int>(n_labels, 0)},
            std::vector<int>(config.permutations, 0), 1, config.alpha));
        return result;
    }
    if (config.rho != 0.0) {
        const bool js_family = config.metric == DistanceMetric::JensenShannon;
        family.points = js_family ? probability_rows(family.patterns, family.size, n_labels)
                                  : normalized_rows(family.patterns, family.size, n_labels);
        if (js_family) {
            family.entropies.resize(family.size);
            for (int i = 0; i < family.size; ++i) {
                family.entropies[i] =
                    shannon_entropy(row_of(family.points, i, n_labels), n_labels);
            }
        }
    }

    // Test half: observed counts, then the same counts under every permutation.
    // Each sample is counted separately; the family-wise test uses their sum.
    const int n_test = static_cast<int>(split.test.size());
    const int n_samples = data.sample_count();
    const std::vector<int> test_histograms =
        subset_rows(candidates.observed_histograms, split.test, n_labels);
    std::vector<SampleTestSet> sets(n_samples);
    for (const int position : split.test) {
        SampleTestSet& set = sets[data.sample_of(candidates.centers[position])];
        set.positions.push_back(position);
        set.rows.push_back(candidates.rows[position]);
    }
    std::vector<std::vector<int>> observed_by_sample(n_samples, std::vector<int>(family.size, 0));
    family.observed.assign(family.size, 0);
    for (int sample = 0; sample < n_samples; ++sample) {
        SampleTestSet& set = sets[sample];
        set.histograms = subset_rows(candidates.observed_histograms, set.positions, n_labels);
        if (set.rows.empty()) continue;
        count_family(set.histograms, static_cast<int>(set.rows.size()), n_labels, family,
                     config, config.threads, observed_by_sample[sample].data());
        for (int c = 0; c < family.size; ++c) family.observed[c] += observed_by_sample[sample][c];
    }
    if (n_samples > 1) {
        std::cerr << "Testing centers by sample:";
        for (int sample = 0; sample < n_samples; ++sample) {
            std::cerr << ' ' << data.sample_names[sample] << '=' << sets[sample].rows.size();
        }
        std::cerr << '\n';
    }

    const LabelShufflePlan test_plan = restrict_to_test_tiles(shuffle_plan, data, config);
    result.null_shuffle_group_count = static_cast<int>(test_plan.groups.size());
    result.null_mixed_group_count = test_plan.mixed_group_count;
    result.null_exchangeable_cell_count = test_plan.exchangeable_cell_count;
    const NullDraws draws = null_family_counts(
        data, config, test_plan, candidates.neighbors, sets, family, observed_by_sample);
    const std::vector<int>& null_counts = draws.pooled;
    const FamilyScores scores = score_family(
        family.observed, family.first_index, null_counts, config.permutations);

    // Report the most significant candidates, skipping near-duplicates (with a
    // positive rho, neighbouring classes are the same motif seen twice) and
    // candidates that do not recur in the test half. Both filters only ever
    // remove rejections, so the family-wise guarantee still holds.
    const double min_separation = 2.0 * config.rho;
    const bool js = config.metric == DistanceMetric::JensenShannon;
    const RowMatcher test_rows(test_histograms, n_test, n_labels, config.rho, config.metric);
    std::vector<int> reported;
    std::vector<std::vector<unsigned char>> reported_matches;
    std::vector<int> reported_support;
    for (const int c : scores.order) {
        if (static_cast<int>(reported.size()) >= config.max_motifs) break;
        const bool distinct = config.rho == 0.0 ||
            std::all_of(reported.begin(), reported.end(), [&](const int kept) {
                const double* a = row_of(family.points, c, n_labels);
                const double* b = row_of(family.points, kept, n_labels);
                double distance2 = 0.0;
                if (js) {
                    distance2 = jensen_shannon_divergence_from_entropy(
                        a, family.entropies[c], b, family.entropies[kept], n_labels);
                } else {
                    for (int d = 0; d < n_labels; ++d) {
                        const double diff = a[d] - b[d];
                        distance2 += diff * diff;
                    }
                }
                return distance2 > min_separation * min_separation + 1e-12;
            });
        if (!distinct) continue;
        const int* pattern = row_of(family.patterns, c, n_labels);
        std::vector<unsigned char> matched = test_rows.match(pattern);
        const int support = disjoint_support(matched, split.test, candidates, n_cells);
        if (support < config.min_support) continue;
        reported.push_back(c);
        reported_matches.push_back(std::move(matched));
        reported_support.push_back(support);
    }

    for (int rank = 0; rank < static_cast<int>(reported.size()); ++rank) {
        const int c = reported[rank];
        const int* pattern = row_of(family.patterns, c, n_labels);
        const int* column =
            null_counts.data() + static_cast<std::size_t>(c) * config.permutations;
        MotifTest test;
        test.rank = rank + 1;
        test.observed = family.observed[c];
        test.pattern.assign(pattern, pattern + n_labels);
        test.null_counts.assign(column, column + config.permutations);
        test.p_value = scores.p_raw[c];
        test.p_value_adjusted = scores.p_adjusted[c];
        test.null_mean = scores.null_mean[c];
        test.null_sd = scores.null_sd[c];
        test.lift = scores.lift[c];
        if (test.null_sd > 0.0) {
            test.z_score = (test.observed - test.null_mean) / test.null_sd;
        } else if (test.observed != test.null_mean) {
            const double inf = std::numeric_limits<double>::infinity();
            test.z_score = test.observed > test.null_mean ? inf : -inf;
        }
        test.significant = test.p_value_adjusted <= config.alpha;
        if (test.significant) ++result.n_significant_motifs;
        const std::vector<unsigned char>& matched = reported_matches[rank];
        test.disjoint_support = reported_support[rank];
        const MatchedShape shape =
            matched_shape(test_histograms, matched, n_test, n_labels);
        test.mean_cells = shape.mean_cells;
        test.mean_composition = shape.mean_composition;
        test.match_mask.assign(n_cells, 0);
        for (int i = 0; i < n_test; ++i) {
            if (matched[i]) test.match_mask[candidates.centers[split.test[i]]] = 1;
        }
        // Several samples: the same motif sample by sample. These per-sample
        // p-values are descriptive; only the pooled adjusted p-value controls
        // the family-wise error.
        for (int sample = 0; n_samples > 1 && sample < n_samples; ++sample) {
            const SampleTestSet& set = sets[sample];
            SampleEvidence evidence;
            evidence.observed = observed_by_sample[sample][c];
            evidence.null_mean = static_cast<double>(draws.sum[sample][c]) / config.permutations;
            const double mean_sq =
                static_cast<double>(draws.sum_sq[sample][c]) / config.permutations;
            evidence.null_sd =
                std::sqrt(std::max(0.0, mean_sq - evidence.null_mean * evidence.null_mean));
            evidence.lift = evidence.null_mean > 0.0
                ? evidence.observed / evidence.null_mean
                : (evidence.observed > 0 ? std::numeric_limits<double>::infinity() : 1.0);
            evidence.p_value = (1.0 + draws.at_least[sample][c]) / (config.permutations + 1.0);
            if (!set.rows.empty()) {
                const int sample_rows = static_cast<int>(set.rows.size());
                const std::vector<unsigned char> sample_matched =
                    RowMatcher(set.histograms, sample_rows, n_labels, config.rho, config.metric)
                        .match(test.pattern.data());
                evidence.disjoint_support =
                    disjoint_support(sample_matched, set.positions, candidates, n_cells);
                const MatchedShape sample_shape =
                    matched_shape(set.histograms, sample_matched, sample_rows, n_labels);
                evidence.mean_cells = sample_shape.mean_cells;
                evidence.mean_composition = sample_shape.mean_composition;
            }
            if (evidence.p_value <= config.alpha) ++test.replicated_in;
            test.per_sample.push_back(evidence);
        }
        // Spread of the per-sample compositions, over the samples that hold the
        // motif at all. It says whether every sample sees the same mixture.
        std::vector<const SampleEvidence*> seen;
        for (const SampleEvidence& evidence : test.per_sample) {
            if (evidence.observed > 0 && !evidence.mean_composition.empty()) seen.push_back(&evidence);
        }
        if (seen.size() > 1) {
            test.composition_sd.assign(n_labels, 0.0);
            for (int d = 0; d < n_labels; ++d) {
                double mean = 0.0;
                for (const SampleEvidence* evidence : seen) mean += evidence->mean_composition[d];
                mean /= seen.size();
                double variance = 0.0;
                for (const SampleEvidence* evidence : seen) {
                    const double diff = evidence->mean_composition[d] - mean;
                    variance += diff * diff;
                }
                test.composition_sd[d] = std::sqrt(variance / (seen.size() - 1));
            }
        }
        result.motif_tests.push_back(std::move(test));
    }
    if (result.motif_tests.empty()) {
        // The family was tested but nothing recurred often enough to report.
        // Downstream readers still expect one row.
        result.motif_tests.push_back(summarize_motif_test(
            {0, 0, std::vector<int>(n_labels, 0)},
            std::vector<int>(config.permutations, 0), 1, config.alpha));
    }
    result.tested_motifs = static_cast<int>(result.motif_tests.size());
    result.selected_rank = 1;
    return result;
}

void apply_stage_error_control(FreezeSweepResult& sweep, const AnalysisConfig& config) {
    std::vector<Result*> stages;
    for (FreezeStageResult& stage : sweep.stages) stages.push_back(&stage.result);
    if (config.stage_error_control == StageErrorControl::Holm) {
        // Holm step-down: sort raw p-values, scale by remaining count, keep monotone.
        std::vector<int> order(stages.size());
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](const int a, const int b) {
            return stages[a]->stage_p_value_raw < stages[b]->stage_p_value_raw;
        });
        double running = 0.0;
        for (int rank = 0; rank < static_cast<int>(order.size()); ++rank) {
            Result& stage = *stages[order[rank]];
            running = std::max(running, std::min(
                1.0, static_cast<double>(order.size() - rank) * stage.stage_p_value_raw));
            stage.stage_p_value_adjusted = running;
        }
    } else {
        for (Result* stage : stages) stage->stage_p_value_adjusted = stage->stage_p_value_raw;
    }
    for (Result* stage : stages) {
        stage->stage_significant = stage->stage_p_value_adjusted <= config.alpha;
    }
}

}  // namespace

double elapsed_seconds(const Clock::time_point start, const Clock::time_point end) {
    return std::chrono::duration<double>(end - start).count();
}

int normalize_thread_count(int requested, const int n) {
    if (requested <= 0) requested = static_cast<int>(std::thread::hardware_concurrency());
    if (requested <= 0) requested = 1;
    return std::min(requested, std::max(1, n));
}

std::string trim(std::string s) {
    const auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const auto last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

std::string lower_ascii(std::string s) {
    for (char& c : s) {
        if (c >= 'A' && c <= 'Z') c = static_cast<char>(c - 'A' + 'a');
    }
    return s;
}

int fewrs_resample_count(const double failure_probability) {
    if (failure_probability <= 0.0 || failure_probability >= 1.0) {
        throw std::runtime_error(
            "FewRS failure probability must be greater than 0 and less than 1");
    }
    const double numerator = std::log(1.0 / failure_probability);
    const double denominator = -std::log1p(-failure_probability);
    return static_cast<int>(std::ceil(numerator / denominator));
}

double shannon_entropy(const double* p, const int k) {
    double out = 0.0;
    for (int d = 0; d < k; ++d) {
        if (p[d] > 0.0) out -= p[d] * std::log(p[d]);
    }
    return out;
}

double jensen_shannon_divergence_from_entropy(const double* p,
                                              const double hp,
                                              const double* q,
                                              const double hq,
                                              const int k) {
    double hm = 0.0;
    for (int d = 0; d < k; ++d) {
        const double m = 0.5 * (p[d] + q[d]);
        if (m > 0.0) hm -= m * std::log(m);
    }
    const double div = hm - 0.5 * (hp + hq);
    return div > 0.0 ? div : 0.0;
}

void js_row_terms(const double* p, const int k, double* half, double* self) {
    for (int d = 0; d < k; ++d) {
        // The same m as the pairwise formula forms: 0.5 * (p + 0) and 0.5 * (p + p).
        const double m_half = 0.5 * (p[d] + 0.0);
        const double m_self = 0.5 * (p[d] + p[d]);
        half[d] = m_half > 0.0 ? m_half * std::log(m_half) : 0.0;
        self[d] = m_self > 0.0 ? m_self * std::log(m_self) : 0.0;
    }
}

Result most_frequent_pattern_test(const Dataset& data, const AnalysisConfig& config) {
    validate_analysis_config(config);
    const std::vector<int> permutable_cells =
        non_frozen_cells(data.labels, frozen_label_mask(data, config.frozen_cell_types));
    if (permutable_cells.empty()) {
        throw std::runtime_error("All cell types are frozen; no motif centers remain");
    }
    // Validate block coordinates before the spatial neighborhood index is built.
    const LabelShufflePlan shuffle_plan = make_label_shuffle_plan(data, permutable_cells, config);

    if (!config.frozen_cell_types.empty()) {
        std::cerr << "Frozen cell types:";
        for (const std::string& name : config.frozen_cell_types) std::cerr << ' ' << name;
        std::cerr << " (" << data.labels.size() - permutable_cells.size() << " fixed cells, "
                  << permutable_cells.size() << " permutable cells)\n";
    }

    const bool covering = config.neighborhood_mode == NeighborhoodMode::Covering;
    std::vector<int> centers = permutable_cells;
    std::vector<int> rows(centers.size());
    std::iota(rows.begin(), rows.end(), 0);
    Neighbors neighbors;
    if (covering) {
        // Cover all cells, so neighbourhoods are built around every cell first.
        std::vector<int> all_cells(data.labels.size());
        std::iota(all_cells.begin(), all_cells.end(), 0);
        neighbors = build_neighbors(data, all_cells, config.radius, config.threads);
        centers = covering_centers(neighbors, permutable_cells, static_cast<int>(data.labels.size()));
        rows = centers;
    } else {
        neighbors = build_neighbors(data, centers, config.radius, config.threads);
    }

    const std::vector<int> observed_histograms = label_histograms(
        neighbors, rows, data.labels, static_cast<int>(data.label_names.size()), config.threads);
    const CandidateNeighborhoods candidates{centers, neighbors, rows, observed_histograms};
    Result result = config.statistic == MotifStatistic::MinP
        ? run_significance_test(data, config, permutable_cells, shuffle_plan, candidates)
        : run_pattern_test_prepared(data, config, permutable_cells, shuffle_plan, candidates, false);
    result.active_label_count = static_cast<int>(data.label_names.size()) -
        static_cast<int>(config.frozen_cell_types.size());
    if (covering) {
        std::cerr << "Covering neighborhoods: " << centers.size() << " centers cover "
                  << result.covered_cell_count << '/' << data.labels.size() << " cells\n";
        if (result.covered_cell_count < static_cast<int>(data.labels.size())) {
            std::cerr << "  Uncovered cells have no eligible center within radius\n";
        }
    }
    return result;
}

FreezeSweepResult run_freeze_sweep(const Dataset& data, const AnalysisConfig& config) {
    validate_analysis_config(config);
    if (!config.freeze_by_abundance) {
        throw std::runtime_error("Internal error: run_freeze_sweep requires freeze-by-abundance");
    }
    if (!config.frozen_cell_types.empty()) {
        throw std::runtime_error("Cumulative freeze sweep cannot start with manually frozen types");
    }
    const int n = static_cast<int>(data.labels.size());
    const int n_labels = static_cast<int>(data.label_names.size());
    if (n_labels < 2) {
        throw std::runtime_error("Cumulative freeze sweep requires at least two observed cell types");
    }

    FreezeSweepResult sweep;
    sweep.label_counts.assign(n_labels, 0);
    for (const int label : data.labels) ++sweep.label_counts[label];
    sweep.freeze_order = abundance_order(data, sweep.label_counts);
    const int last_informative_stage = n_labels - 2;
    const int max_frozen = config.max_freeze_stages < 0
        ? last_informative_stage
        : std::min(config.max_freeze_stages, last_informative_stage);

    // Neighbourhoods and observed histograms of every cell are shared by all stages.
    std::vector<int> all_cells(n);
    std::iota(all_cells.begin(), all_cells.end(), 0);
    const LabelShufflePlan initial_shuffle_plan = make_label_shuffle_plan(data, all_cells, config);
    const auto neighbor_start = Clock::now();
    const Neighbors all_neighbors = build_neighbors(data, all_cells, config.radius, config.threads);
    sweep.neighbor_build_seconds = elapsed_seconds(neighbor_start);
    const auto histogram_start = Clock::now();
    const std::vector<int> all_histograms =
        label_histograms(all_neighbors, all_cells, data.labels, n_labels, config.threads);
    sweep.observed_histogram_seconds = elapsed_seconds(histogram_start);

    std::vector<unsigned char> frozen_labels(n_labels, 0);
    for (int stage = 0; stage <= max_frozen; ++stage) {
        AnalysisConfig stage_config = config;
        stage_config.freeze_by_abundance = false;
        stage_config.max_freeze_stages = -1;
        stage_config.frozen_cell_types.clear();
        for (int rank = 0; rank < stage; ++rank) {
            frozen_labels[sweep.freeze_order[rank]] = 1;
            stage_config.frozen_cell_types.push_back(data.label_names[sweep.freeze_order[rank]]);
        }
        const std::vector<int> permutable_cells = non_frozen_cells(data.labels, frozen_labels);
        const bool covering = config.neighborhood_mode == NeighborhoodMode::Covering;
        const std::vector<int> centers =
            covering ? covering_centers(all_neighbors, permutable_cells, n) : permutable_cells;

        std::cerr << "Freeze stage " << stage << '/' << max_frozen << ": " << stage
                  << " frozen types, " << centers.size() << " candidate centers";
        if (covering) std::cerr << " from " << permutable_cells.size() << " eligible cells";
        std::cerr << '\n';

        const auto stage_start = Clock::now();
        const LabelShufflePlan shuffle_plan = stage == 0
            ? initial_shuffle_plan
            : make_label_shuffle_plan(data, permutable_cells, stage_config);
        Result result = run_pattern_test_prepared(
            data, stage_config, permutable_cells, shuffle_plan,
            {centers, all_neighbors, centers, select_histogram_rows(all_histograms, centers, n_labels)},
            true);
        result.freeze_stage = stage;
        result.active_label_count = n_labels - stage;
        result.stage_p_value_raw = result.motif_tests.front().p_value;

        FreezeStageResult& stage_result = sweep.stages.emplace_back();
        stage_result.stage_index = stage;
        stage_result.frozen_label_ids.assign(sweep.freeze_order.begin(),
                                             sweep.freeze_order.begin() + stage);
        stage_result.result = std::move(result);
        stage_result.timings.analysis_seconds = elapsed_seconds(stage_start);
        stage_result.timings.total_seconds = stage_result.timings.analysis_seconds;
    }

    apply_stage_error_control(sweep, config);
    return sweep;
}

int count_matches(const Result& result) {
    const auto& mask = result.selected().match_mask;
    return static_cast<int>(std::count(mask.begin(), mask.end(), 1));
}
