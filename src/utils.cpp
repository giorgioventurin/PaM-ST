#include "utils.h"
#include "KDTree.h"
#include "null_model.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
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

const double* row_of(const std::vector<double>& flat, const int index, const int k) {
    return flat.data() + static_cast<std::size_t>(index) * k;
}

const int* row_of(const std::vector<int>& flat, const int index, const int k) {
    return flat.data() + static_cast<std::size_t>(index) * k;
}

// ---------------------------------------------------------------------------
// Neighbourhoods and count vectors

// neighbors[i] lists the spots of the same sample within `radius` of spot i,
// using a hash grid with bucket size `radius`. Every spot is a centre.
Neighbors build_neighbors(const Dataset& data, const double radius, const int threads) {
    const std::vector<std::array<double, 2>>& coords = data.coords;
    if (radius < 0.0) throw std::runtime_error("radius must be non-negative");
    const int n = static_cast<int>(coords.size());
    Neighbors neighbors(n);
    if (radius == 0.0) {
        for (int i = 0; i < n; ++i) neighbors[i].push_back(i);
        return neighbors;
    }

    auto bucket_of = [&](const int cell) {
        const std::array<double, 2>& p = coords[cell];
        return GridKey{data.sample_of(cell),
                       static_cast<std::int64_t>(std::floor(p[0] / radius)),
                       static_cast<std::int64_t>(std::floor(p[1] / radius))};
    };
    std::unordered_map<GridKey, std::vector<int>, GridKeyHash> grid;
    grid.reserve(static_cast<std::size_t>(n * 1.3));
    for (int i = 0; i < n; ++i) grid[bucket_of(i)].push_back(i);

    const double radius2 = radius * radius;
    parallel_for(n, threads, [&](const int lo, const int hi) {
        for (int i = lo; i < hi; ++i) {
            const auto& center = coords[i];
            const GridKey c = bucket_of(i);
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

// Row i of the flattened result is the count vector c_i of spot i's neighbourhood.
std::vector<int> count_vectors(const Neighbors& neighbors,
                               const std::vector<int>& labels,
                               const int n_labels,
                               const int threads) {
    std::vector<int> histograms(neighbors.size() * n_labels, 0);
    parallel_for(static_cast<int>(neighbors.size()), threads, [&](const int lo, const int hi) {
        for (int i = lo; i < hi; ++i) {
            int* row = histograms.data() + static_cast<std::size_t>(i) * n_labels;
            for (const int cell : neighbors[i]) ++row[labels[cell]];
        }
    });
    return histograms;
}

// ---------------------------------------------------------------------------
// Compositions

// Count vectors equal up to a common factor have the same normalised vector, so
// they are one composition. Its key is the vector divided by its gcd.
std::vector<int> composition_key(const int* row, const int k) {
    std::vector<int> key(row, row + k);
    int divisor = 0;
    for (int d = 0; d < k; ++d) divisor = std::gcd(divisor, row[d]);
    if (divisor > 1) {
        for (int& value : key) value /= divisor;
    }
    return key;
}

// Distinct compositions among n count vectors.
struct CompositionClasses {
    std::vector<int> rows;            // first count vector seen for each class, flattened
    std::vector<int> multiplicities;  // neighbourhoods in each class
    std::vector<int> class_of_row;    // class of every row, in row order
};

CompositionClasses composition_classes(const std::vector<int>& histograms,
                                       const int n,
                                       const int k) {
    std::unordered_map<std::vector<int>, int, VectorHash> index_of;
    index_of.reserve(static_cast<std::size_t>(n * 1.3));
    CompositionClasses out;
    for (int i = 0; i < n; ++i) {
        const int* row = row_of(histograms, i, k);
        const auto [it, inserted] = index_of.try_emplace(
            composition_key(row, k), static_cast<int>(out.multiplicities.size()));
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

// Each row divided by its Euclidean norm: the normalised count vectors c~.
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

// The compositions of one labelling, indexed so that the neighbourhoods within
// rho of any of them come back directly. Every row of a class has the same
// normalised vector, so a range query over the classes answers it for all
// their rows at once; scanning every row per class instead is what used to
// make the search the cost of the whole test.
class CompositionIndex {
public:
    CompositionIndex(const std::vector<int>& histograms, const int n, const int n_labels,
                     const double rho, const int threads)
        : n_labels_(n_labels), rho_(rho),
          classes_(composition_classes(histograms, n, n_labels)) {
        const int m = size();
        frequency_ = classes_.multiplicities;

        // Rows of each class, in row order.
        class_start_.assign(m + 1, 0);
        for (const int cls : classes_.class_of_row) ++class_start_[cls + 1];
        for (int i = 0; i < m; ++i) class_start_[i + 1] += class_start_[i];
        rows_by_class_.assign(n, 0);
        std::vector<int> cursor(class_start_.begin(), class_start_.end() - 1);
        for (int row = 0; row < n; ++row) {
            rows_by_class_[cursor[classes_.class_of_row[row]]++] = row;
        }

        if (rho == 0.0) return;
        points_ = normalized_rows(classes_.rows, m, n_labels);
        tree_ = std::make_unique<KDTreeND>(points_, n_labels, classes_.multiplicities);
        if (rho >= std::sqrt(2.0)) {
            frequency_.assign(m, n);  // no two normalised count vectors are farther apart
        } else {
            frequency_ = tree_->weighted_neighbor_counts(rho * rho, threads);
        }
    }

    // The tree points into this object, so it must stay where it was built.
    CompositionIndex(const CompositionIndex&) = delete;
    CompositionIndex& operator=(const CompositionIndex&) = delete;

    int size() const { return static_cast<int>(classes_.multiplicities.size()); }
    const int* pattern(const int cls) const { return row_of(classes_.rows, cls, n_labels_); }
    // f: the neighbourhoods within rho of the class.
    int frequency(const int cls) const { return frequency_[cls]; }

    // Those neighbourhoods themselves, in spot order.
    void rows_near(const int cls, std::vector<int>& out) const {
        out.clear();
        if (rho_ == 0.0) {
            out.assign(rows_by_class_.begin() + class_start_[cls],
                       rows_by_class_.begin() + class_start_[cls + 1]);
            return;
        }
        std::vector<int> hits;
        tree_->neighbors_within(row_of(points_, cls, n_labels_), rho_ * rho_, hits);
        for (const int hit : hits) {
            out.insert(out.end(), rows_by_class_.begin() + class_start_[hit],
                       rows_by_class_.begin() + class_start_[hit + 1]);
        }
        std::sort(out.begin(), out.end());
    }

private:
    int n_labels_ = 0;
    double rho_ = 0.0;
    CompositionClasses classes_;
    std::vector<int> frequency_;
    std::vector<int> class_start_;
    std::vector<int> rows_by_class_;
    std::vector<double> points_;  // empty when rho == 0
    std::unique_ptr<KDTreeND> tree_;
};

// The recurrence of a pattern: a set of its occurrences whose neighbourhoods
// share no spot, taken greedily in spot order. Generation stamps instead of a
// cleared mask, so each call costs only the spots it touches.
class SeparatedCounter {
public:
    explicit SeparatedCounter(const int n_cells) : stamp_(n_cells, 0) {}

    int count(const std::vector<int>& rows, const Neighbors& neighbors) {
        ++generation_;
        int kept = 0;
        for (const int row : rows) {
            const std::vector<int>& cells = neighbors[row];
            if (std::any_of(cells.begin(), cells.end(),
                            [&](const int cell) { return stamp_[cell] == generation_; })) {
                continue;
            }
            for (const int cell : cells) stamp_[cell] = generation_;
            ++kept;
        }
        return kept;
    }

private:
    std::vector<int> stamp_;
    int generation_ = 0;
};

// ---------------------------------------------------------------------------
// The patterns of one labelling

// Patterns proposed by any labelling, keyed by composition so that the same
// mixture found in different labellings is one entry: the universe U.
struct PatternUniverse {
    std::unordered_map<std::vector<int>, int, VectorHash> id_of;
    std::vector<int> rows;  // size x n_labels, the first count vector seen
    int size = 0;

    int intern(const int* pattern, const int n_labels) {
        const auto [it, inserted] = id_of.try_emplace(composition_key(pattern, n_labels), size);
        if (inserted) {
            rows.insert(rows.end(), pattern, pattern + n_labels);
            ++size;
        }
        return it->second;
    }
};

struct Proposal {
    int id = 0;          // in the universe
    int cls = 0;         // in the labelling's composition index
    int recurrence = 0;
};

// P^(j): the count vectors of labelling j holding at least min_types labels
// and recurring in at least min_support separate places. The frequency bounds
// the recurrence from above, so it rules most classes out before the scan.
std::vector<Proposal> propose(const CompositionIndex& index,
                              const Neighbors& neighbors,
                              const AnalysisConfig& config,
                              const int n_labels,
                              SeparatedCounter& counter,
                              PatternUniverse& universe) {
    std::vector<Proposal> out;
    std::vector<int> rows;
    for (int cls = 0; cls < index.size(); ++cls) {
        if (index.frequency(cls) < config.min_support) continue;
        const int* pattern = index.pattern(cls);
        const int labels = static_cast<int>(
            std::count_if(pattern, pattern + n_labels, [](const int v) { return v > 0; }));
        if (labels < config.min_types) continue;
        index.rows_near(cls, rows);
        const int recurrence = counter.count(rows, neighbors);
        if (recurrence < config.min_support) continue;
        out.push_back({universe.intern(pattern, n_labels), cls, recurrence});
    }
    return out;
}

// f^(j)_v for every v of the universe, on one labelling's count vectors.
void count_universe(const std::vector<int>& histograms, const int n, const int n_labels,
                    const PatternUniverse& universe, const std::vector<double>& points,
                    const AnalysisConfig& config, int* out) {
    const double rho = config.rho;
    const CompositionClasses classes = composition_classes(histograms, n, n_labels);
    const int m = static_cast<int>(classes.multiplicities.size());

    if (rho == 0.0) {
        // Exact proportional match: look each pattern up among the classes.
        std::unordered_map<std::vector<int>, int, VectorHash> multiplicity_of;
        multiplicity_of.reserve(static_cast<std::size_t>(m * 1.3));
        for (int i = 0; i < m; ++i) {
            multiplicity_of.emplace(composition_key(row_of(classes.rows, i, n_labels), n_labels),
                                    classes.multiplicities[i]);
        }
        for (int c = 0; c < universe.size; ++c) {
            const auto it = multiplicity_of.find(
                composition_key(row_of(universe.rows, c, n_labels), n_labels));
            out[c] = it == multiplicity_of.end() ? 0 : it->second;
        }
        return;
    }

    if (rho >= std::sqrt(2.0)) {
        std::fill(out, out + universe.size, n);
        return;
    }
    const std::vector<double> class_points = normalized_rows(classes.rows, m, n_labels);
    const KDTreeND tree(class_points, n_labels, classes.multiplicities);
    parallel_for(universe.size, config.threads, [&](const int lo, const int hi) {
        for (int c = lo; c < hi; ++c) {
            out[c] = tree.weighted_count_within(row_of(points, c, n_labels), rho * rho);
        }
    });
}

// How many standard deviations a count lies above the mean of a pattern's
// counts; zero when the counts do not vary, so the score is always finite.
double z_of(const double count, const double mean, const double sd) {
    return sd > 0.0 ? (count - mean) / sd : 0.0;
}

// The matching neighbourhoods as they actually look: their mean size and mean
// proportions per label. Every neighbourhood counts once, whatever its size,
// so one dense neighbourhood cannot outvote a dozen sparse ones.
void describe_occurrences(const std::vector<int>& histograms,
                          const std::vector<int>& rows,
                          const int n_labels,
                          Pattern& pattern) {
    pattern.mean_cells = 0.0;
    pattern.mean_composition.assign(n_labels, 0.0);
    int hits = 0;
    for (const int i : rows) {
        const int* row = row_of(histograms, i, n_labels);
        const int total = std::accumulate(row, row + n_labels, 0);
        if (total <= 0) continue;
        ++hits;
        pattern.mean_cells += total;
        for (int d = 0; d < n_labels; ++d) {
            pattern.mean_composition[d] += static_cast<double>(row[d]) / total;
        }
    }
    if (hits > 0) {
        pattern.mean_cells /= hits;
        for (double& value : pattern.mean_composition) value /= hits;
    }
}

void validate_analysis_config(const AnalysisConfig& config) {
    validate_null_model_config(config);
    if (config.rho < 0.0) throw std::runtime_error("rho must be non-negative");
    if (config.permutations < 1) throw std::runtime_error("permutations must be at least 1");
    if (config.alpha <= 0.0 || config.alpha >= 1.0) {
        throw std::runtime_error("alpha must be greater than 0 and less than 1");
    }
    if (config.min_support < 1) throw std::runtime_error("min-support must be at least 1");
    if (config.min_types < 0) throw std::runtime_error("min-types must be non-negative");
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

// The test of docs/pattern_test.tex. Every labelling - the observed one and
// each tile shuffle - proposes its own patterns under the same rule, so the
// real tissue is simply one of the m+1 labellings, and the p-value of an
// observed pattern is the share of labellings whose best z-score reaches its
// own.
Result run_pattern_test(const Dataset& data, const AnalysisConfig& config) {
    validate_analysis_config(config);
    // Validate tile coordinates before the neighbourhood index is built.
    const LabelShufflePlan shuffle_plan = make_label_shuffle_plan(data, config);
    std::cerr << "Tiles: " << shuffle_plan.groups.size() << " occupied, "
              << shuffle_plan.mixed_group_count << " with mixed labels ("
              << shuffle_plan.exchangeable_cell_count << " exchangeable cells)"
              << (shuffle_plan.exchangeable_cell_count == 0 ? "; no label can move" : "")
              << '\n';

    const int n_cells = static_cast<int>(data.labels.size());
    const int n_labels = static_cast<int>(data.label_names.size());
    const int m = config.permutations;
    const int pool = m + 1;
    const int step = std::max(1, pool / 20);  // progress lines, whatever the budget
    const Neighbors neighbors = build_neighbors(data, config.radius, config.threads);

    // Pass 1: every labelling proposes its own patterns P^(j).
    PatternUniverse universe;
    std::vector<std::vector<Proposal>> proposed(pool);
    std::vector<int> observed_histograms;
    std::unique_ptr<CompositionIndex> observed;
    {
        SeparatedCounter counter(n_cells);
        std::vector<int> labels = data.labels;
        std::mt19937 rng(config.seed);
        for (int j = 0; j < pool; ++j) {
            if (j > 0) shuffle_labels_within_groups(shuffle_plan, labels, rng);
            std::vector<int> histograms = count_vectors(neighbors, labels, n_labels, config.threads);
            auto index = std::make_unique<CompositionIndex>(
                histograms, n_cells, n_labels, config.rho, config.threads);
            proposed[j] = propose(*index, neighbors, config, n_labels, counter, universe);
            if (j == 0) {
                observed = std::move(index);
                observed_histograms = std::move(histograms);
            }
            if ((j + 1) % step == 0 || j + 1 == pool) {
                std::cerr << "  search " << (j + 1) << '/' << pool << ", universe "
                          << universe.size << '\n';
            }
        }
    }
    std::cerr << "Observed tissue proposes " << proposed[0].size() << " patterns; "
              << universe.size << " proposed by any labelling\n";

    // Pass 2: replay the same labellings and count every pattern of the
    // universe on each - a pattern needs its count everywhere, not only where
    // it was proposed. counts[v * pool + j] is f^(j)_v.
    const std::vector<double> points = config.rho > 0.0
        ? normalized_rows(universe.rows, universe.size, n_labels)
        : std::vector<double>{};
    std::vector<int> counts(static_cast<std::size_t>(universe.size) * pool, 0);
    {
        std::vector<int> labels = data.labels;
        std::mt19937 rng(config.seed);
        std::vector<int> draw(universe.size, 0);
        for (int j = 0; j < pool; ++j) {
            if (j > 0) shuffle_labels_within_groups(shuffle_plan, labels, rng);
            count_universe(count_vectors(neighbors, labels, n_labels, config.threads),
                           n_cells, n_labels, universe, points, config, draw.data());
            for (int c = 0; c < universe.size; ++c) {
                counts[static_cast<std::size_t>(c) * pool + j] = draw[c];
            }
            if ((j + 1) % step == 0 || j + 1 == pool) {
                std::cerr << "  count " << (j + 1) << '/' << pool << '\n';
            }
        }
    }

    // z^(j)_v: each count against the pattern's counts on all m+1 labellings,
    // with their population mean and standard deviation.
    std::vector<double> z(static_cast<std::size_t>(universe.size) * pool);
    for (int c = 0; c < universe.size; ++c) {
        const int* row = counts.data() + static_cast<std::size_t>(c) * pool;
        const double mean = static_cast<double>(std::accumulate(row, row + pool, 0LL)) / pool;
        double variance = 0.0;
        for (int j = 0; j < pool; ++j) variance += (row[j] - mean) * (row[j] - mean);
        const double sd = std::sqrt(variance / pool);
        for (int j = 0; j < pool; ++j) {
            z[static_cast<std::size_t>(c) * pool + j] = z_of(row[j], mean, sd);
        }
    }

    // T^(j): each labelling's best score over its own patterns. A labelling
    // with no patterns has no best score and never counts against one.
    std::vector<double> best(pool, -std::numeric_limits<double>::infinity());
    for (int j = 0; j < pool; ++j) {
        for (const Proposal& p : proposed[j]) {
            best[j] = std::max(best[j], z[static_cast<std::size_t>(p.id) * pool + j]);
        }
    }

    // Report every pattern of P^(0), highest z first.
    std::vector<Proposal> order = proposed[0];
    std::sort(order.begin(), order.end(), [&](const Proposal& a, const Proposal& b) {
        const double za = z[static_cast<std::size_t>(a.id) * pool];
        const double zb = z[static_cast<std::size_t>(b.id) * pool];
        return za != zb ? za > zb : a.id < b.id;
    });
    Result result;
    result.universe_size = universe.size;
    std::vector<int> rows;
    for (const Proposal& p : order) {
        Pattern pattern;
        pattern.counts.assign(observed->pattern(p.cls), observed->pattern(p.cls) + n_labels);
        pattern.frequency = counts[static_cast<std::size_t>(p.id) * pool];
        pattern.recurrence = p.recurrence;
        pattern.z_score = z[static_cast<std::size_t>(p.id) * pool];
        pattern.p_value = static_cast<double>(std::count_if(
            best.begin(), best.end(), [&](const double t) { return t >= pattern.z_score; })) / pool;
        pattern.significant = pattern.p_value <= config.alpha;
        observed->rows_near(p.cls, rows);
        describe_occurrences(observed_histograms, rows, n_labels, pattern);
        if (pattern.significant) {
            pattern.occurrences = rows;
            ++result.n_significant;
        }
        result.patterns.push_back(std::move(pattern));
    }
    return result;
}
