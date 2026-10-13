#include "utils.h"
#include "KDTree.h"
#include "null_model.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
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

struct VectorHash {
    std::size_t operator()(const std::vector<int>& v) const {
        std::size_t h = 1469598103934665603ULL;
        for (const int x : v) {
            h ^= static_cast<std::uint64_t>(x) + 0x9E3779B97F4A7C15ULL;
            h *= 1099511628211ULL;
        }
        return h;
    }
};

using Neighbors = std::vector<std::vector<int>>;

// Distinct compositions (count vectors equal up to a common factor).
struct CompositionClasses {
    std::vector<int> rows;            // first count vector seen for each class, flattened
    std::vector<int> multiplicities;  // neighbourhoods in each class
    std::vector<int> class_of_row;    // class of every row, in row order
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

// ---------------------------------------------------------------------------
// Compositions and motif search

int composition_divisor(const int* row, const int k) {
    int divisor = 0;
    for (int d = 0; d < k; ++d) divisor = std::gcd(divisor, row[d]);
    return divisor;
}

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
    std::unordered_map<std::vector<int>, int, VectorHash> index_of;
    index_of.reserve(static_cast<std::size_t>(n * 1.3));
    CompositionClasses out;
    for (int i = 0; i < n; ++i) {
        const int* row = row_of(histograms, i, k);
        std::vector<int> key(row, row + k);
        if (const int divisor = composition_divisor(row, k); divisor > 1) {
            for (int& value : key) value /= divisor;
        }
        const auto [it, inserted] =
            index_of.try_emplace(std::move(key), static_cast<int>(out.multiplicities.size()));
        out.class_of_row.push_back(it->second);
        if (inserted) {
            out.rows.insert(out.rows.end(), row, row + k);
            out.multiplicities.push_back(1);
        } else {
            ++out.multiplicities[it->second];
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

// Every distinct composition of one histogram matrix, with the number of
// neighbourhoods within rho of it.
struct ScoredClasses {
    CompositionClasses classes;
    std::vector<int> counts;
    std::vector<double> points;  // query rows; empty when rho == 0
    int size() const { return static_cast<int>(counts.size()); }
};

ScoredClasses score_classes(const std::vector<int>& histograms,
                            const int n,
                            const int n_labels,
                            const double rho,
                            const int threads) {
    ScoredClasses out;
    out.classes = composition_classes(histograms, n, n_labels);
    out.counts = out.classes.multiplicities;
    if (rho == 0.0) return out;

    const int m = out.size();
    out.points = normalized_rows(out.classes.rows, m, n_labels);
    if (rho >= std::sqrt(2.0)) {
        out.counts.assign(m, n);
    } else {
        out.counts = KDTreeND(out.points, n_labels, out.classes.multiplicities)
                         .weighted_neighbor_counts(rho * rho, threads);
    }
    return out;
}

// Marks the histogram rows lying within rho of the pattern.
std::vector<unsigned char> motif_match_rows(const std::vector<int>& histograms,
                                            const int n,
                                            const int k,
                                            const std::vector<int>& pattern,
                                            const double rho) {
    std::vector<unsigned char> out(n, 0);
    if (rho == 0.0) {
        for (int i = 0; i < n; ++i) {
            out[i] = same_composition(row_of(histograms, i, k), pattern.data(), k);
        }
        return out;
    }

    const double rho2 = rho * rho;
    double pnorm2 = 0.0;
    for (const int v : pattern) pnorm2 += static_cast<double>(v) * v;
    std::vector<double> pn(k, 0.0);
    if (pnorm2 > 0.0) {
        const double inv = 1.0 / std::sqrt(pnorm2);
        for (int d = 0; d < k; ++d) pn[d] = pattern[d] * inv;
    }
    for (int i = 0; i < n; ++i) {
        const int* row = row_of(histograms, i, k);
        double norm2 = 0.0;
        for (int d = 0; d < k; ++d) norm2 += static_cast<double>(row[d]) * row[d];
        const double inv = norm2 > 0.0 ? 1.0 / std::sqrt(norm2) : 0.0;
        double dist2 = 0.0;
        for (int d = 0; d < k; ++d) {
            const double diff = row[d] * inv - pn[d];
            dist2 += diff * diff;
        }
        out[i] = dist2 <= rho2;
    }
    return out;
}

// A pattern needs at least min_types labels present. The recurrence floor is
// checked by the caller; here `count` is the frequency, which bounds it from
// above and so rules most classes out cheaply.
bool passes_family_filters(const int* pattern, const int count, const int n_labels,
                           const AnalysisConfig& config) {
    if (count < config.min_support) return false;
    if (config.min_types <= 0) return true;
    int types = 0;
    for (int d = 0; d < n_labels; ++d) {
        if (pattern[d] > 0) ++types;
    }
    return types >= config.min_types;
}

// ---------------------------------------------------------------------------
// Counting a fixed set of patterns on any labelling

struct CandidateFamily {
    std::vector<int> patterns;   // size x n_labels raw count vectors
    std::vector<double> points;  // size x n_labels query rows
    int size = 0;
};

// Counts every candidate in one (permuted) tissue.
void count_family(const std::vector<int>& histograms, const int n, const int n_labels,
                  const CandidateFamily& family, const AnalysisConfig& config, int* out) {
    const double rho = config.rho;
    const CompositionClasses classes = composition_classes(histograms, n, n_labels);
    const int m = static_cast<int>(classes.multiplicities.size());

    if (rho == 0.0) {
        // Exact proportional match: look each candidate up among the classes.
        auto reduced_key = [n_labels](const int* row) {
            std::vector<int> key(row, row + n_labels);
            if (const int divisor = composition_divisor(row, n_labels); divisor > 1) {
                for (int& value : key) value /= divisor;
            }
            return key;
        };
        std::unordered_map<std::vector<int>, int, VectorHash> multiplicity_of;
        multiplicity_of.reserve(static_cast<std::size_t>(m * 1.3));
        for (int i = 0; i < m; ++i) {
            multiplicity_of.emplace(reduced_key(row_of(classes.rows, i, n_labels)),
                                    classes.multiplicities[i]);
        }
        for (int c = 0; c < family.size; ++c) {
            const auto it = multiplicity_of.find(reduced_key(row_of(family.patterns, c, n_labels)));
            out[c] = it == multiplicity_of.end() ? 0 : it->second;
        }
        return;
    }

    if (rho >= std::sqrt(2.0)) {
        std::fill(out, out + family.size, n);
        return;
    }
    const std::vector<double> points = normalized_rows(classes.rows, m, n_labels);
    const KDTreeND tree(points, n_labels, classes.multiplicities);
    parallel_for(family.size, config.threads, [&](const int lo, const int hi) {
        for (int c = lo; c < hi; ++c) {
            out[c] = tree.weighted_count_within(row_of(family.points, c, n_labels), rho * rho);
        }
    });
}

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

// ---------------------------------------------------------------------------
// Permutation test

// How many standard deviations a count lies above the mean of a pattern's
// counts; zero when the counts do not vary, so the score is always finite.
double z_of(const double count, const double mean, const double sd) {
    return sd > 0.0 ? (count - mean) / sd : 0.0;
}

// The z-score of the observed count among all of a pattern's counts - the
// observed one together with its null counts - with their population mean and
// standard deviation.
double pooled_z_score(const int observed, const std::vector<int>& null_counts) {
    const double pool = static_cast<double>(null_counts.size()) + 1.0;
    double sum = observed;
    for (const int v : null_counts) sum += v;
    const double mean = sum / pool;
    double variance = (observed - mean) * (observed - mean);
    for (const int v : null_counts) variance += (v - mean) * (v - mean);
    return z_of(observed, mean, std::sqrt(variance / pool));
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
}

void validate_analysis_config(const AnalysisConfig& config) {
    validate_null_model_config(config);
    if (config.rho < 0.0) throw std::runtime_error("rho must be non-negative");
    if (config.permutations < 1) throw std::runtime_error("permutations must be at least 1");
    if (config.alpha <= 0.0 || config.alpha >= 1.0) {
        throw std::runtime_error("alpha must be greater than 0 and less than 1");
    }
    if (config.max_motifs < 1) throw std::runtime_error("max-motifs must be at least 1");
    if (config.min_support < 1) throw std::runtime_error("min-support must be at least 1");
    if (config.min_types < 0) throw std::runtime_error("min-types must be non-negative");
}

// Diagnostics: tile counts and neighbourhood coverage.
Result start_result(const Dataset& data,
                    const LabelShufflePlan& shuffle_plan,
                    const CandidateNeighborhoods& candidates) {
    Result result;
    result.null_block_count = static_cast<int>(shuffle_plan.groups.size());
    result.null_mixed_group_count = shuffle_plan.mixed_group_count;
    result.null_exchangeable_cell_count = shuffle_plan.exchangeable_cell_count;
    result.candidate_center_count = static_cast<int>(candidates.centers.size());
    std::cerr << "Block null: " << result.null_block_count << " occupied tiles, "
              << result.null_mixed_group_count << " mixed tiles, "
              << result.null_exchangeable_cell_count << " exchangeable cells"
              << (result.null_exchangeable_cell_count == 0 ? " (no labels can change)" : "")
              << '\n';
    record_neighborhood_coverage(result, candidates, static_cast<int>(data.labels.size()));
    return result;
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
int disjoint_support(const std::vector<unsigned char>& matched,
                     const std::vector<int>& positions,
                     const CandidateNeighborhoods& candidates,
                     const int n_cells) {
    std::vector<unsigned char> used(n_cells, 0);
    int support = 0;
    for (int i = 0; i < static_cast<int>(positions.size()); ++i) {
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

// ---------------------------------------------------------------------------
// Max-search test
//
// The whole search is the statistic. Every labelling - the observed one and
// each permutation - proposes its own candidates under the same rule and is
// scored against the same pooled reference, so the real tissue is simply one of
// the m+1 labellings.

// Patterns proposed by any labelling, keyed by their composition class so that
// the same mixture found in different labellings is one entry.
struct PatternUniverse {
    std::unordered_map<std::vector<int>, int, VectorHash> id_of;
    std::vector<int> rows;  // size x n_labels, the first raw vector seen
    int size = 0;

    int intern(const int* pattern, const int n_labels) {
        std::vector<int> key(pattern, pattern + n_labels);
        if (const int divisor = composition_divisor(pattern, n_labels); divisor > 1) {
            for (int& value : key) value /= divisor;
        }
        const auto [it, inserted] = id_of.try_emplace(std::move(key), size);
        if (inserted) {
            rows.insert(rows.end(), pattern, pattern + n_labels);
            ++size;
        }
        return it->second;
    }
};

// The candidates one labelling proposes: the same filters applied to every
// centre, identical for the observed tissue and for a permutation.
//
// The support filters are answered class by class rather than row by row: every
// row of a class has the same normalised vector, so the classes within rho of a
// candidate give its matching rows directly. Scanning all rows per candidate
// instead is what made this loop the cost of the whole test.
std::vector<int> search_candidates(const std::vector<int>& histograms,
                                   const std::vector<int>& positions,
                                   const CandidateNeighborhoods& candidates,
                                   const AnalysisConfig& config,
                                   const int n_labels,
                                   const int n_cells,
                                   PatternUniverse& universe) {
    const int n = static_cast<int>(positions.size());
    const ScoredClasses classes =
        score_classes(histograms, n, n_labels, config.rho, config.threads);
    const int class_count = classes.size();
    std::vector<int> ids;

    if (config.rho == 0.0) {
        // Exact matching keeps the straightforward path.
        for (int i = 0; i < class_count; ++i) {
            const int* pattern = row_of(classes.classes.rows, i, n_labels);
            if (!passes_family_filters(pattern, classes.counts[i], n_labels, config)) continue;
            const std::vector<unsigned char> matched = motif_match_rows(
                histograms, n, n_labels, std::vector<int>(pattern, pattern + n_labels),
                config.rho);
            if (disjoint_support(matched, positions, candidates, n_cells) < config.min_support) {
                continue;
            }
            ids.push_back(universe.intern(pattern, n_labels));
        }
        return ids;
    }

    // Rows of each class, in row order.
    std::vector<int> class_start(class_count + 1, 0);
    for (const int cls : classes.classes.class_of_row) ++class_start[cls + 1];
    for (int i = 0; i < class_count; ++i) class_start[i + 1] += class_start[i];
    std::vector<int> rows_by_class(n, 0);
    {
        std::vector<int> cursor(class_start.begin(), class_start.end() - 1);
        for (int row = 0; row < n; ++row) {
            rows_by_class[cursor[classes.classes.class_of_row[row]]++] = row;
        }
    }

    const KDTreeND tree(classes.points, n_labels, classes.classes.multiplicities);
    const double rho2 = config.rho * config.rho;
    // Generation stamps instead of a cleared mask, so each candidate costs only
    // the cells it actually touches.
    std::vector<int> stamp(n_cells, 0);
    int generation = 0;
    std::vector<int> hits;
    std::vector<int> rows;
    for (int i = 0; i < class_count; ++i) {
        const int* pattern = row_of(classes.classes.rows, i, n_labels);
        if (!passes_family_filters(pattern, classes.counts[i], n_labels, config)) continue;
        hits.clear();
        tree.neighbors_within(row_of(classes.points, i, n_labels), rho2, hits);
        rows.clear();
        for (const int cls : hits) {
            rows.insert(rows.end(), rows_by_class.begin() + class_start[cls],
                        rows_by_class.begin() + class_start[cls + 1]);
        }
        std::sort(rows.begin(), rows.end());
        ++generation;
        int support = 0;
        for (const int row : rows) {
            const std::vector<int>& cells = candidates.neighbors[candidates.rows[positions[row]]];
            if (std::any_of(cells.begin(), cells.end(),
                            [&](const int cell) { return stamp[cell] == generation; })) {
                continue;
            }
            for (const int cell : cells) stamp[cell] = generation;
            ++support;
        }
        if (support < config.min_support) continue;
        ids.push_back(universe.intern(pattern, n_labels));
    }
    return ids;
}

Result run_max_search_test(const Dataset& data,
                           const AnalysisConfig& config,
                           const LabelShufflePlan& shuffle_plan,
                           const CandidateNeighborhoods& candidates) {
    const int n_cells = static_cast<int>(data.labels.size());
    const int n_labels = static_cast<int>(data.label_names.size());
    const int n_samples = data.sample_count();
    const int m = config.permutations;
    const int pool = m + 1;
    Result result = start_result(data, shuffle_plan, candidates);

    const int step = std::max(1, pool / 20);  // progress lines, whatever the budget
    const int n_centers = static_cast<int>(candidates.rows.size());
    std::vector<int> all_positions(n_centers);
    std::iota(all_positions.begin(), all_positions.end(), 0);
    // Centres grouped by sample, for the per-sample evidence.
    std::vector<std::vector<int>> sample_positions(n_samples);
    for (int position = 0; position < n_centers; ++position) {
        sample_positions[data.sample_of(candidates.centers[position])].push_back(position);
    }

    // Pass 1: every labelling proposes its own candidates.
    PatternUniverse universe;
    std::vector<std::vector<int>> proposed(pool);
    {
        std::vector<int> labels = data.labels;
        std::mt19937 rng(config.seed);
        for (int j = 0; j < pool; ++j) {
            if (j > 0) shuffle_labels_within_groups(shuffle_plan, labels, rng);
            const std::vector<int> histograms =
                label_histograms(candidates.neighbors, candidates.rows, labels, n_labels,
                                 config.threads);
            proposed[j] = search_candidates(histograms, all_positions, candidates, config,
                                            n_labels, n_cells, universe);
            if ((j + 1) % step == 0 || j + 1 == pool) {
                std::cerr << "  search " << (j + 1) << '/' << pool << ", universe "
                          << universe.size << '\n';
            }
        }
    }
    result.universe_size = universe.size;
    result.candidate_family_size = static_cast<int>(proposed[0].size());
    std::cerr << "Observed tissue proposes " << proposed[0].size() << " candidates; "
              << universe.size << " patterns proposed by any labelling\n";

    // The universe as a query family, so every labelling can be counted on all
    // of it - a pattern needs its count everywhere, not only where it was
    // proposed.
    CandidateFamily family;
    family.size = universe.size;
    family.patterns = universe.rows;
    if (config.rho > 0.0) family.points = normalized_rows(family.patterns, family.size, n_labels);

    // Pass 2: replay the same labellings and count every pattern of the
    // universe on each, pooled over samples and kept per sample as well.
    std::vector<int> counts(static_cast<std::size_t>(family.size) * pool, 0);
    std::vector<std::vector<int>> observed_by_sample(n_samples, std::vector<int>(family.size, 0));
    std::vector<std::vector<long long>> sample_sum(n_samples, std::vector<long long>(family.size, 0));
    std::vector<std::vector<long long>> sample_sum_sq(n_samples, std::vector<long long>(family.size, 0));
    std::vector<std::vector<int>> sample_at_least(n_samples, std::vector<int>(family.size, 0));
    {
        std::vector<int> labels = data.labels;
        std::mt19937 rng(config.seed);
        std::vector<int> draw(family.size, 0);
        for (int j = 0; j < pool; ++j) {
            if (j > 0) shuffle_labels_within_groups(shuffle_plan, labels, rng);
            const std::vector<int> histograms =
                label_histograms(candidates.neighbors, candidates.rows, labels, n_labels,
                                 config.threads);
            for (int sample = 0; sample < n_samples; ++sample) {
                if (sample_positions[sample].empty()) continue;
                count_family(subset_rows(histograms, sample_positions[sample], n_labels),
                             static_cast<int>(sample_positions[sample].size()), n_labels,
                             family, config, draw.data());
                for (int c = 0; c < family.size; ++c) {
                    counts[static_cast<std::size_t>(c) * pool + j] += draw[c];
                    if (j == 0) {
                        observed_by_sample[sample][c] = draw[c];
                        continue;
                    }
                    sample_sum[sample][c] += draw[c];
                    sample_sum_sq[sample][c] += static_cast<long long>(draw[c]) * draw[c];
                    if (draw[c] >= observed_by_sample[sample][c]) ++sample_at_least[sample][c];
                }
            }
            if ((j + 1) % step == 0 || j + 1 == pool) {
                std::cerr << "  count " << (j + 1) << '/' << pool << '\n';
            }
        }
    }

    // Scores: the z-score of every labelling's count of every pattern among
    // that pattern's counts on all m+1 labellings, with their population mean
    // and standard deviation.
    std::vector<double> z(static_cast<std::size_t>(family.size) * pool);
    for (int c = 0; c < family.size; ++c) {
        const int* row = counts.data() + static_cast<std::size_t>(c) * pool;
        const double mean = static_cast<double>(std::accumulate(row, row + pool, 0LL)) / pool;
        double variance = 0.0;
        for (int j = 0; j < pool; ++j) variance += (row[j] - mean) * (row[j] - mean);
        const double sd = std::sqrt(variance / pool);
        for (int j = 0; j < pool; ++j) {
            z[static_cast<std::size_t>(c) * pool + j] = z_of(row[j], mean, sd);
        }
    }

    // Each labelling keeps only its best score T^(j); those m+1 values are the
    // reference the observed patterns are judged against. A labelling with no
    // patterns has no best score and never counts against an observed one.
    std::vector<double> best(pool, -std::numeric_limits<double>::infinity());
    for (int j = 0; j < pool; ++j) {
        for (const int c : proposed[j]) {
            best[j] = std::max(best[j], z[static_cast<std::size_t>(c) * pool + j]);
        }
    }

    // Report the observed tissue's candidates, highest z first, skipping
    // near-duplicates of something already reported.
    std::vector<int> order = proposed[0];
    std::sort(order.begin(), order.end(), [&](const int a, const int b) {
        const double za = z[static_cast<std::size_t>(a) * pool];
        const double zb = z[static_cast<std::size_t>(b) * pool];
        return za != zb ? za > zb : a < b;
    });
    const std::vector<int> observed_histograms = candidates.observed_histograms;
    const double min_separation = 2.0 * config.rho;
    std::vector<int> reported;
    for (const int c : order) {
        if (static_cast<int>(reported.size()) >= config.max_motifs) break;
        const bool distinct = config.rho == 0.0 ||
            std::all_of(reported.begin(), reported.end(), [&](const int kept) {
                const double* a = row_of(family.points, c, n_labels);
                const double* b = row_of(family.points, kept, n_labels);
                double distance2 = 0.0;
                for (int d = 0; d < n_labels; ++d) {
                    const double diff = a[d] - b[d];
                    distance2 += diff * diff;
                }
                return distance2 > min_separation * min_separation + 1e-12;
            });
        if (distinct) reported.push_back(c);
    }

    for (int rank = 0; rank < static_cast<int>(reported.size()); ++rank) {
        const int c = reported[rank];
        const int* row = counts.data() + static_cast<std::size_t>(c) * pool;
        const double observed_z = z[static_cast<std::size_t>(c) * pool];
        MotifTest test;
        test.rank = rank + 1;
        test.observed = row[0];
        test.pattern.assign(row_of(family.patterns, c, n_labels),
                            row_of(family.patterns, c, n_labels) + n_labels);
        test.null_counts.assign(row + 1, row + pool);
        // Pooled rank of the observed count among all m+1 labellings.
        test.p_value = static_cast<double>(
            std::count_if(row, row + pool,
                          [&](const int value) { return value >= row[0]; })) / pool;
        // How often any labelling's best score was at least this high.
        test.p_value_adjusted = static_cast<double>(std::count_if(
            best.begin(), best.end(), [&](const double t) { return t >= observed_z; })) / pool;
        const double mean = std::accumulate(row + 1, row + pool, 0LL) / static_cast<double>(m);
        double variance = 0.0;
        for (int j = 1; j < pool; ++j) variance += (row[j] - mean) * (row[j] - mean);
        test.null_mean = mean;
        test.null_sd = std::sqrt(variance / m);
        test.lift = mean > 0.0 ? test.observed / mean
                               : (test.observed > 0 ? std::numeric_limits<double>::infinity() : 1.0);
        test.z_score = pooled_z_score(test.observed, test.null_counts);
        test.significant = test.p_value_adjusted <= config.alpha;
        if (test.significant) ++result.n_significant_motifs;

        const std::vector<unsigned char> matched = motif_match_rows(
            observed_histograms, n_centers, n_labels, test.pattern, config.rho);
        test.disjoint_support = disjoint_support(matched, all_positions, candidates, n_cells);
        test.match_mask.assign(n_cells, 0);
        for (int i = 0; i < n_centers; ++i) {
            if (matched[i]) test.match_mask[candidates.centers[i]] = 1;
        }
        const MatchedShape shape =
            matched_shape(observed_histograms, matched, n_centers, n_labels);
        test.mean_cells = shape.mean_cells;
        test.mean_composition = shape.mean_composition;

        for (int sample = 0; n_samples > 1 && sample < n_samples; ++sample) {
            SampleEvidence evidence;
            evidence.observed = observed_by_sample[sample][c];
            evidence.null_mean = static_cast<double>(sample_sum[sample][c]) / m;
            const double mean_sq = static_cast<double>(sample_sum_sq[sample][c]) / m;
            evidence.null_sd =
                std::sqrt(std::max(0.0, mean_sq - evidence.null_mean * evidence.null_mean));
            evidence.lift = evidence.null_mean > 0.0
                ? evidence.observed / evidence.null_mean
                : (evidence.observed > 0 ? std::numeric_limits<double>::infinity() : 1.0);
            evidence.p_value = (1.0 + sample_at_least[sample][c]) / pool;
            if (!sample_positions[sample].empty()) {
                const std::vector<int> sample_histograms =
                    subset_rows(observed_histograms, sample_positions[sample], n_labels);
                const int rows = static_cast<int>(sample_positions[sample].size());
                const std::vector<unsigned char> sample_matched = motif_match_rows(
                    sample_histograms, rows, n_labels, test.pattern, config.rho);
                evidence.disjoint_support =
                    disjoint_support(sample_matched, sample_positions[sample], candidates, n_cells);
                const MatchedShape sample_shape =
                    matched_shape(sample_histograms, sample_matched, rows, n_labels);
                evidence.mean_cells = sample_shape.mean_cells;
                evidence.mean_composition = sample_shape.mean_composition;
            }
            if (evidence.p_value <= config.alpha) ++test.replicated_in;
            test.per_sample.push_back(evidence);
        }
        result.motif_tests.push_back(std::move(test));
    }
    if (result.motif_tests.empty()) {
        // Nothing proposed: one empty row, so the reports still have a line.
        MotifTest empty;
        empty.pattern.assign(n_labels, 0);
        empty.null_counts.assign(m, 0);
        result.motif_tests.push_back(std::move(empty));
    }
    result.tested_motifs = static_cast<int>(result.motif_tests.size());
    result.selected_rank = 1;
    return result;
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

Result most_frequent_pattern_test(const Dataset& data, const AnalysisConfig& config) {
    validate_analysis_config(config);
    // Validate block coordinates before the spatial neighborhood index is built.
    const LabelShufflePlan shuffle_plan = make_label_shuffle_plan(data, config);

    // Every spot is a centre.
    std::vector<int> centers(data.labels.size());
    std::iota(centers.begin(), centers.end(), 0);
    const std::vector<int>& rows = centers;
    const Neighbors neighbors = build_neighbors(data, centers, config.radius, config.threads);
    const std::vector<int> observed_histograms = label_histograms(
        neighbors, rows, data.labels, static_cast<int>(data.label_names.size()), config.threads);
    const CandidateNeighborhoods candidates{centers, neighbors, rows, observed_histograms};
    return run_max_search_test(data, config, shuffle_plan, candidates);
}

int count_matches(const Result& result) {
    const auto& mask = result.selected().match_mask;
    return static_cast<int>(std::count(mask.begin(), mask.end(), 1));
}
