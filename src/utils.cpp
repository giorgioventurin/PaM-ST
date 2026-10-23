#include "utils.h"
#include "KDTree.h"
#include "null_model.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
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

// Each row divided by its L1 norm, the number of spots it counts: the
// normalised count vectors c~, which hold the proportion of each label.
// Dividing rather than multiplying by the inverse keeps proportional rows
// bit-identical, since each proportion is then the correctly rounded quotient.
std::vector<double> normalized_rows(const std::vector<int>& rows, const int n, const int k) {
    std::vector<double> out(static_cast<std::size_t>(n) * k, 0.0);
    for (int i = 0; i < n; ++i) {
        const int* row = row_of(rows, i, k);
        const double total = std::accumulate(row, row + k, 0.0);
        if (total == 0.0) continue;
        double* dst = out.data() + static_cast<std::size_t>(i) * k;
        for (int d = 0; d < k; ++d) dst[d] = row[d] / total;
    }
    return out;
}

// The compositions of one labelling, indexed so that the occurrences of any
// count vector come back directly. Every row of a class has the same
// normalised vector, so a range query over the classes answers it for all
// their rows at once.
class CompositionIndex {
public:
    CompositionIndex(const std::vector<int>& histograms, const int n, const int n_labels,
                     const double rho)
        : n_labels_(n_labels), rho_(rho),
          classes_(composition_classes(histograms, n, n_labels)) {
        const int m = size();
        // Rows of each class, in row order.
        class_start_.assign(m + 1, 0);
        for (const int cls : classes_.class_of_row) ++class_start_[cls + 1];
        for (int i = 0; i < m; ++i) class_start_[i + 1] += class_start_[i];
        rows_by_class_.assign(n, 0);
        std::vector<int> cursor(class_start_.begin(), class_start_.end() - 1);
        for (int row = 0; row < n; ++row) {
            rows_by_class_[cursor[classes_.class_of_row[row]]++] = row;
        }
        if (rho == 0.0) {
            class_of_key_.reserve(static_cast<std::size_t>(m * 1.3));
            for (int c = 0; c < m; ++c) class_of_key_.emplace(composition_key(pattern(c), n_labels), c);
            return;
        }
        points_ = normalized_rows(classes_.rows, m, n_labels);
        tree_ = std::make_unique<KDTreeND>(points_, n_labels, classes_.multiplicities);
    }

    // The tree points into this object, so it must stay where it was built.
    CompositionIndex(const CompositionIndex&) = delete;
    CompositionIndex& operator=(const CompositionIndex&) = delete;

    int size() const { return static_cast<int>(classes_.multiplicities.size()); }
    const int* pattern(const int cls) const { return row_of(classes_.rows, cls, n_labels_); }

    // occ(v): the rows whose count vector lies within rho of v, in spot order.
    // `key` is the composition key of v, used when rho == 0, and `point` its
    // normalised vector, used otherwise.
    void occurrences(const std::vector<int>& key, const double* point,
                     std::vector<int>& hits, std::vector<int>& out) const {
        out.clear();
        if (rho_ > 0.0) {
            within(point, hits, out);
            return;
        }
        const auto it = class_of_key_.find(key);
        if (it != class_of_key_.end()) rows_of_class(it->second, out);
    }

    // occ(v) for v one of this labelling's own compositions.
    void class_occurrences(const int cls, std::vector<int>& hits, std::vector<int>& out) const {
        out.clear();
        if (rho_ > 0.0) {
            within(row_of(points_, cls, n_labels_), hits, out);
        } else {
            rows_of_class(cls, out);
        }
    }

private:
    void rows_of_class(const int cls, std::vector<int>& out) const {
        out.insert(out.end(), rows_by_class_.begin() + class_start_[cls],
                   rows_by_class_.begin() + class_start_[cls + 1]);
    }

    void within(const double* point, std::vector<int>& hits, std::vector<int>& out) const {
        hits.clear();
        tree_->neighbors_within(point, rho_ * rho_, hits);
        for (const int hit : hits) rows_of_class(hit, out);
        std::sort(out.begin(), out.end());
    }

    int n_labels_ = 0;
    double rho_ = 0.0;
    CompositionClasses classes_;
    std::vector<int> class_start_;
    std::vector<int> rows_by_class_;
    std::unordered_map<std::vector<int>, int, VectorHash> class_of_key_;  // rho == 0 only
    std::vector<double> points_;                                          // rho > 0 only
    std::unique_ptr<KDTreeND> tree_;
};

// Working space for the greedy recurrence, one per thread. The per-spot
// arrays are left zeroed after every call.
struct GreedyScratch {
    explicit GreedyScratch(const int n_cells) : count(n_cells, 0), start(n_cells, 0) {}
    std::vector<int> count, start;        // per spot
    std::vector<int> touched, members;    // spots met, and who holds each
    std::vector<int> mark, degree, removed;
    std::vector<std::size_t> adj_start;
    std::vector<int> adj;
    std::vector<unsigned char> alive;
    std::vector<std::pair<int, int>> heap;
    std::vector<int> rows, hits;          // occurrences and class hits
};

// The greedy recurrence of Algorithm 1. The conflict graph joins two
// occurrences whose neighbourhoods share a spot; the greedy repeatedly takes a
// vertex of minimum degree in the remaining graph and deletes it with its
// neighbours, degrees being recomputed as vertices go. Ties go to the lowest
// index, a fixed rule, so the value depends on the occurrences alone.
int greedy_recurrence(const std::vector<int>& occ, const Neighbors& neighbors, GreedyScratch& w) {
    const int f = static_cast<int>(occ.size());
    if (f <= 1) return f;

    // For every spot, the occurrences whose neighbourhood holds it.
    w.touched.clear();
    for (const int i : occ) {
        for (const int s : neighbors[i]) {
            if (w.count[s]++ == 0) w.touched.push_back(s);
        }
    }
    int total = 0;
    for (const int s : w.touched) {
        w.start[s] = total;
        total += w.count[s];
        w.count[s] = w.start[s];  // becomes the fill cursor, then the end
    }
    w.members.resize(total);
    for (int a = 0; a < f; ++a) {
        for (const int s : neighbors[occ[a]]) w.members[w.count[s]++] = a;
    }

    // The conflict graph, as adjacency lists.
    w.adj_start.assign(f + 1, 0);
    w.adj.clear();
    w.mark.assign(f, -1);
    for (int a = 0; a < f; ++a) {
        w.mark[a] = a;
        for (const int s : neighbors[occ[a]]) {
            for (int k = w.start[s]; k < w.count[s]; ++k) {
                const int b = w.members[k];
                if (w.mark[b] != a) {
                    w.mark[b] = a;
                    w.adj.push_back(b);
                }
            }
        }
        w.adj_start[a + 1] = w.adj.size();
    }
    for (const int s : w.touched) w.count[s] = 0;

    // Minimum-degree greedy, with a lazy heap: a popped entry whose degree is
    // out of date is skipped, and every decrement pushes a fresh one.
    w.degree.resize(f);
    w.alive.assign(f, 1);
    w.heap.clear();
    for (int a = 0; a < f; ++a) {
        w.degree[a] = static_cast<int>(w.adj_start[a + 1] - w.adj_start[a]);
        w.heap.emplace_back(w.degree[a], a);
    }
    const auto later = std::greater<std::pair<int, int>>();
    std::make_heap(w.heap.begin(), w.heap.end(), later);
    int kept = 0;
    while (!w.heap.empty()) {
        std::pop_heap(w.heap.begin(), w.heap.end(), later);
        const auto [d, u] = w.heap.back();
        w.heap.pop_back();
        if (!w.alive[u] || d != w.degree[u]) continue;
        ++kept;
        // Delete u and its neighbours, then lower the degrees around them.
        w.removed.clear();
        w.alive[u] = 0;
        w.removed.push_back(u);
        for (std::size_t k = w.adj_start[u]; k < w.adj_start[u + 1]; ++k) {
            const int b = w.adj[k];
            if (w.alive[b]) {
                w.alive[b] = 0;
                w.removed.push_back(b);
            }
        }
        for (const int x : w.removed) {
            for (std::size_t k = w.adj_start[x]; k < w.adj_start[x + 1]; ++k) {
                const int y = w.adj[k];
                if (!w.alive[y]) continue;
                --w.degree[y];
                w.heap.emplace_back(w.degree[y], y);
                std::push_heap(w.heap.begin(), w.heap.end(), later);
            }
        }
    }
    return kept;
}

// Every count vector proposed by any labelling, keyed by composition so that
// proportional vectors, which share occurrences and recurrence, are one entry:
// the universe U.
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
    int find(const int* pattern, const int n_labels) const {
        const auto it = id_of.find(composition_key(pattern, n_labels));
        if (it == id_of.end()) throw std::runtime_error("Internal error: pattern missing from the universe");
        return it->second;
    }
};

// The score z of a pattern on one labelling, how many standard deviations its
// recurrence R lies above the mean of its m+1 recurrences, kept exactly as
// (pool R - S1) / sqrt(pool S2 - S1^2), with S1 and S2 the sums of the
// recurrences and of their squares. It is zero when they do not vary.
// Exactness matters because ties are common: a pattern seen on one labelling
// alone scores sqrt(m) whatever its recurrence, and m is the largest a score
// can be, so rounding would decide p-values.
struct Score {
    long long num = 0;    // pool R - S1
    __int128 spread = 0;  // pool S2 - S1^2

    int sign() const { return spread == 0 || num == 0 ? 0 : (num > 0 ? 1 : -1); }
    double value() const {
        return sign() == 0 ? 0.0 : static_cast<double>(num) / std::sqrt(static_cast<double>(spread));
    }
};

using Wide = std::pair<unsigned __int128, unsigned __int128>;  // high and low halves

Wide multiply_wide(const unsigned __int128 a, const unsigned __int128 b) {
    const unsigned __int128 mask = ~std::uint64_t{0};
    const unsigned __int128 a0 = a & mask, a1 = a >> 64, b0 = b & mask, b1 = b >> 64;
    const unsigned __int128 p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
    const unsigned __int128 mid = (p00 >> 64) + (p01 & mask) + (p10 & mask);
    return {p11 + (p01 >> 64) + (p10 >> 64) + (mid >> 64), (p00 & mask) | (mid << 64)};
}

// a >= b, decided exactly: for scores of one sign, |a| >= |b| exactly when
// num_a^2 spread_b >= num_b^2 spread_a.
bool at_least(const Score& a, const Score& b) {
    const int sa = a.sign(), sb = b.sign();
    if (sa != sb) return sa > sb;
    if (sa == 0) return true;
    const auto square = [](const long long x) {
        const auto u = static_cast<unsigned __int128>(x < 0 ? -x : x);
        return u * u;
    };
    const Wide lhs = multiply_wide(square(a.num), static_cast<unsigned __int128>(b.spread));
    const Wide rhs = multiply_wide(square(b.num), static_cast<unsigned __int128>(a.spread));
    return sa > 0 ? lhs >= rhs : lhs <= rhs;
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
    if (config.min_recurrence < 1) throw std::runtime_error("min-recurrence must be at least 1");
    if (config.permutations < 1) throw std::runtime_error("permutations must be at least 1");
    if (config.alpha <= 0.0 || config.alpha >= 1.0) {
        throw std::runtime_error("alpha must be greater than 0 and less than 1");
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

// The test of docs/pam-st.tex. The patterns of a labelling are its count
// vectors that recur at least kappa times there (kappa = 1 keeps every count
// vector), the statistic is the greedy recurrence, and the p-value of an
// observed pattern is the share of the m+1 labellings whose best score reaches
// its own.
// The labellings are replayed in three passes from one seeded stream, so that
// only per-pattern sums are kept: memory grows with |U|, not with |U| (m+1).
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
    const int kappa = config.min_recurrence;
    const int step = std::max(1, pool / 20);  // progress lines, whatever the budget
    const Neighbors neighbors = build_neighbors(data, config.radius, config.threads);

    // Labelling j > 0 is an independent tile shuffle of the observed labelling.
    // The shuffles come from one seeded stream, whose use depends only on the
    // tile sizes, so every pass sees the same m+1 labellings.
    auto for_each_labelling = [&](const char* what, auto&& visit) {
        std::mt19937 rng(config.seed);
        std::vector<int> labels;
        for (int j = 0; j < pool; ++j) {
            labels = data.labels;
            if (j > 0) shuffle_labels_within_groups(shuffle_plan, labels, rng);
            visit(j, count_vectors(neighbors, labels, n_labels, config.threads));
            if ((j + 1) % step == 0 || j + 1 == pool) {
                std::cerr << "  " << what << ' ' << (j + 1) << '/' << pool << '\n';
            }
        }
    };

    // The patterns of one labelling: the classes of its index that recur at
    // least kappa times there, as (class, recurrence). Recurrence never exceeds
    // frequency, so rarer classes skip the greedy.
    auto patterns_of = [&](const CompositionIndex& index) {
        const int k = index.size();
        std::vector<int> r(k, 0);
        parallel_for(k, config.threads, [&](const int lo, const int hi) {
            GreedyScratch w(n_cells);
            for (int cls = lo; cls < hi; ++cls) {
                index.class_occurrences(cls, w.hits, w.rows);
                if (static_cast<int>(w.rows.size()) >= kappa) r[cls] = greedy_recurrence(w.rows, neighbors, w);
            }
        });
        std::vector<std::pair<int, int>> patterns;
        for (int cls = 0; cls < k; ++cls) {
            if (r[cls] >= kappa) patterns.emplace_back(cls, r[cls]);
        }
        return patterns;
    };

    // Pass 1: the universe, every pattern of every labelling. Every count
    // vector recurs at least once on its own labelling, so kappa = 1 needs no
    // recurrences here.
    PatternUniverse universe;
    for_each_labelling("patterns", [&](int, const std::vector<int>& histograms) {
        if (kappa == 1) {
            const CompositionClasses classes = composition_classes(histograms, n_cells, n_labels);
            for (int c = 0; c < static_cast<int>(classes.multiplicities.size()); ++c) {
                universe.intern(row_of(classes.rows, c, n_labels), n_labels);
            }
            return;
        }
        const CompositionIndex index(histograms, n_cells, n_labels, config.rho);
        for (const auto& [cls, r] : patterns_of(index)) universe.intern(index.pattern(cls), n_labels);
    });
    const int n_universe = universe.size;
    std::cerr << "Universe: " << n_universe << " patterns\n";
    std::vector<std::vector<int>> keys(n_universe);
    for (int c = 0; c < n_universe; ++c) keys[c] = composition_key(row_of(universe.rows, c, n_labels), n_labels);
    const std::vector<double> points = config.rho > 0.0
        ? normalized_rows(universe.rows, n_universe, n_labels)
        : std::vector<double>(static_cast<std::size_t>(n_universe) * n_labels, 0.0);
    auto recurrence = [&](const CompositionIndex& index, const int c, GreedyScratch& w) {
        index.occurrences(keys[c], row_of(points, c, n_labels), w.hits, w.rows);
        return greedy_recurrence(w.rows, neighbors, w);
    };

    // Pass 2: the recurrence of every pattern of U on every labelling, kept as
    // sums for the mean and standard deviation over the m+1 labellings.
    std::vector<long long> sum(n_universe, 0), sum_sq(n_universe, 0);
    for_each_labelling("recurrences", [&](int, const std::vector<int>& histograms) {
        const CompositionIndex index(histograms, n_cells, n_labels, config.rho);
        parallel_for(n_universe, config.threads, [&](const int lo, const int hi) {
            GreedyScratch w(n_cells);
            for (int c = lo; c < hi; ++c) {
                const long long r = recurrence(index, c, w);
                sum[c] += r;
                sum_sq[c] += r * r;
            }
        });
    });
    auto score = [&](const int c, const long long r) {
        return Score{pool * r - sum[c],
                     static_cast<__int128>(pool) * sum_sq[c] - static_cast<__int128>(sum[c]) * sum[c]};
    };

    // Pass 3: each labelling's best score over its own patterns, and the
    // observed patterns with their scores.
    // best[j] is T^(j); a labelling without patterns, which only an empty
    // tissue has, would have T = -infinity and never count.
    std::vector<Score> best(pool);
    std::vector<unsigned char> has_best(pool, 0);
    std::vector<Pattern> observed;
    std::vector<Score> observed_score;
    std::vector<int> observed_id;
    for_each_labelling("scores", [&](const int j, const std::vector<int>& histograms) {
        const CompositionIndex index(histograms, n_cells, n_labels, config.rho);
        const std::vector<std::pair<int, int>> patterns = patterns_of(index);
        const int k = static_cast<int>(patterns.size());
        std::vector<int> id(k);
        std::vector<Score> z(k);
        for (int a = 0; a < k; ++a) {
            id[a] = universe.find(index.pattern(patterns[a].first), n_labels);
            z[a] = score(id[a], patterns[a].second);
            if (!has_best[j] || !at_least(best[j], z[a])) best[j] = z[a];
            has_best[j] = 1;
        }
        if (j != 0) return;
        observed.assign(k, Pattern{});
        observed_id = id;
        observed_score = z;
        parallel_for(k, config.threads, [&](const int lo, const int hi) {
            std::vector<int> hits, rows;
            for (int a = lo; a < hi; ++a) {
                const int cls = patterns[a].first;
                index.class_occurrences(cls, hits, rows);
                Pattern& p = observed[a];
                p.counts.assign(index.pattern(cls), index.pattern(cls) + n_labels);
                p.frequency = static_cast<int>(rows.size());
                p.recurrence = patterns[a].second;
                p.z_score = z[a].value();
                describe_occurrences(histograms, rows, n_labels, p);
            }
        });
    });
    int at_ceiling = 0;
    for (int j = 0; j < pool; ++j) {
        // z = sqrt(m) exactly when num^2 == m spread.
        const auto magnitude = static_cast<__int128>(best[j].num) * best[j].num;
        if (has_best[j] && best[j].sign() > 0 && magnitude == static_cast<__int128>(m) * best[j].spread) {
            ++at_ceiling;
        }
    }
    std::cerr << "Labellings whose best score is the largest possible, sqrt(m): " << at_ceiling
              << " of " << pool << '\n';

    // p-values, then the observed patterns highest z first.
    Result result;
    result.universe_size = n_universe;
    std::vector<int> order(observed.size());
    std::iota(order.begin(), order.end(), 0);
    for (int a = 0; a < static_cast<int>(observed.size()); ++a) {
        int reached = 0;
        for (int j = 0; j < pool; ++j) reached += has_best[j] && at_least(best[j], observed_score[a]);
        observed[a].p_value = static_cast<double>(reached) / pool;
        observed[a].significant = observed[a].p_value <= config.alpha;
    }
    std::sort(order.begin(), order.end(), [&](const int a, const int b) {
        const Score& za = observed_score[a];
        const Score& zb = observed_score[b];
        if (!at_least(zb, za)) return true;   // za > zb
        if (!at_least(za, zb)) return false;  // za < zb
        return observed_id[a] < observed_id[b];
    });

    // The occurrences of the significant patterns, on the observed labelling.
    const std::vector<int> observed_histograms =
        count_vectors(neighbors, data.labels, n_labels, config.threads);
    const CompositionIndex observed_index(observed_histograms, n_cells, n_labels, config.rho);
    std::vector<int> hits;
    for (const int a : order) {
        Pattern& p = observed[a];
        if (p.significant) {
            observed_index.occurrences(keys[observed_id[a]], row_of(points, observed_id[a], n_labels),
                                       hits, p.occurrences);
            ++result.n_significant;
        }
        result.patterns.push_back(std::move(p));
    }
    return result;
}
