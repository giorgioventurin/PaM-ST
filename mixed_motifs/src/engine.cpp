#include "mixed.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

namespace mixed {
namespace {

struct VectorHash {
    std::size_t operator()(const std::vector<int>& values) const noexcept {
        std::size_t hash = 0x9e3779b9;
        for (int value : values)
            hash ^= static_cast<std::size_t>(value) + 0x9e3779b9 +
                    (hash << 6) + (hash >> 2);
        return hash;
    }
};

// A separate tree per (sample, component) prevents cross-tissue neighbors.
class KDTree {
    struct Node {
        int cell = -1, left = -1, right = -1, smallest_cell = -1;
        double xmin = 0, xmax = 0, ymin = 0, ymax = 0;
    };
    const Dataset& data_;
    std::vector<int> order_;
    std::vector<Node> nodes_;
    int root_ = -1;

    int build(int begin, int end, int depth) {
        if (begin == end) return -1;
        const int middle = begin + (end - begin) / 2;
        const bool split_x = depth % 2 == 0;
        std::nth_element(order_.begin() + begin, order_.begin() + middle,
                         order_.begin() + end, [&](int a, int b) {
            const double u = split_x ? data_.cells[a].x : data_.cells[a].y;
            const double v = split_x ? data_.cells[b].x : data_.cells[b].y;
            return u < v || (u == v && a < b);
        });
        const int index = static_cast<int>(nodes_.size());
        nodes_.push_back(Node{});
        const int cell = order_[middle];
        const int left = build(begin, middle, depth + 1);
        const int right = build(middle + 1, end, depth + 1);
        Node& node = nodes_[index];
        node.cell = node.smallest_cell = cell;
        node.left = left;
        node.right = right;
        node.xmin = node.xmax = data_.cells[cell].x;
        node.ymin = node.ymax = data_.cells[cell].y;
        for (int child : {left, right}) {
            if (child < 0) continue;
            const Node& other = nodes_[child];
            node.xmin = std::min(node.xmin, other.xmin);
            node.xmax = std::max(node.xmax, other.xmax);
            node.ymin = std::min(node.ymin, other.ymin);
            node.ymax = std::max(node.ymax, other.ymax);
            node.smallest_cell = std::min(node.smallest_cell, other.smallest_cell);
        }
        return index;
    }

    double lower_distance(int index, const Cell& query) const {
        if (index < 0) return std::numeric_limits<double>::infinity();
        const Node& node = nodes_[index];
        const double dx = std::max({node.xmin - query.x, 0.0, query.x - node.xmax});
        const double dy = std::max({node.ymin - query.y, 0.0, query.y - node.ymax});
        return std::hypot(dx, dy);
    }

    using Neighbor = std::pair<double, int>;
    using Heap = std::priority_queue<Neighbor>;

    void search(int index, int center, int wanted, double cutoff, Heap& heap) const {
        if (index < 0) return;
        const Node& node = nodes_[index];
        const Cell& query = data_.cells[center];
        const double lower = lower_distance(index, query);
        if (lower > cutoff) return;
        // Index is the secondary distance tie-breaker, including duplicate coordinates.
        if (static_cast<int>(heap.size()) == wanted &&
            Neighbor(lower, node.smallest_cell) > heap.top()) return;
        if (node.cell != center) {
            const Cell& candidate = data_.cells[node.cell];
            const Neighbor neighbor(std::hypot(candidate.x - query.x,
                                                candidate.y - query.y), node.cell);
            if (neighbor.first <= cutoff) {
                if (static_cast<int>(heap.size()) < wanted) heap.push(neighbor);
                else if (neighbor < heap.top()) {
                    heap.pop();
                    heap.push(neighbor);
                }
            }
        }
        int first = node.left, second = node.right;
        const auto child_key = [&](int child) {
            return Neighbor(lower_distance(child, query),
                            child < 0 ? std::numeric_limits<int>::max()
                                      : nodes_[child].smallest_cell);
        };
        if (child_key(second) < child_key(first)) std::swap(first, second);
        search(first, center, wanted, cutoff, heap);
        search(second, center, wanted, cutoff, heap);
    }

public:
    KDTree(const Dataset& data, const std::vector<int>& cells) : data_(data), order_(cells) {
        nodes_.reserve(cells.size());
        root_ = build(0, static_cast<int>(order_.size()), 0);
    }

    bool neighborhood(int center, int k, double cutoff, Neighborhood& result) const {
        if (static_cast<int>(order_.size()) < k) return false;
        Heap heap;
        if (k > 1) search(root_, center, k - 1, cutoff, heap);
        if (static_cast<int>(heap.size()) != k - 1) return false;
        result.center = center;
        result.radius = heap.empty() ? 0.0 : heap.top().first;
        result.members.reserve(k);
        result.members.push_back(center);  // Always include this cell, even at distance ties.
        while (!heap.empty()) {
            result.members.push_back(heap.top().second);
            heap.pop();
        }
        std::sort(result.members.begin(), result.members.end());
        return true;
    }
};

bool same_type_set(const std::vector<int>& a, const std::vector<int>& b) {
    for (std::size_t i = 0; i < a.size(); ++i)
        if ((a[i] != 0) != (b[i] != 0)) return false;
    return true;
}

// Deleting r members gives an inverted index for half-L1 distance <= r.
// Nondecreasing deletion positions enumerate each subhistogram exactly once.
void visit_deletions(std::vector<int>& counts, int remaining, int begin,
                     const std::function<void(const std::vector<int>&)>& visit) {
    if (remaining == 0) {
        visit(counts);
        return;
    }
    for (int type = begin; type < static_cast<int>(counts.size()); ++type) {
        if (counts[type] == 0) continue;
        --counts[type];
        visit_deletions(counts, remaining - 1, type, visit);
        ++counts[type];
    }
}

void next_stamp(std::vector<std::uint32_t>& marks, std::uint32_t& stamp) {
    if (++stamp == 0) {
        std::fill(marks.begin(), marks.end(), 0);
        stamp = 1;
    }
}

struct Histogram {
    std::vector<int> counts, neighborhoods;
};

class HistogramIndex {
    using SignatureIndex = std::unordered_map<std::vector<int>, std::vector<int>, VectorHash>;
    SignatureIndex signatures_;
    std::vector<std::uint32_t> seen_;
    std::uint32_t stamp_ = 0;
    const Config& config_;
    int deletions_;

public:
    std::vector<Histogram> histograms;

    HistogramIndex(const std::vector<NeighborhoodStats>& stats, const Config& config)
        : config_(config), deletions_(std::min(config.max_relabels, config.k)) {
        std::unordered_map<std::vector<int>, int, VectorHash> histogram_ids;
        histogram_ids.reserve(stats.size());
        for (int index = 0; index < static_cast<int>(stats.size()); ++index) {
            if (!stats[index].eligible) continue;
            const auto inserted = histogram_ids.emplace(stats[index].counts,
                                                        static_cast<int>(histograms.size()));
            if (inserted.second) histograms.push_back({stats[index].counts, {}});
            histograms[inserted.first->second].neighborhoods.push_back(index);
        }
        for (int id = 0; id < static_cast<int>(histograms.size()); ++id) {
            std::vector<int> reduced = histograms[id].counts;
            visit_deletions(reduced, deletions_, 0, [&](const std::vector<int>& signature) {
                signatures_[signature].push_back(id);
            });
        }
        seen_.assign(histograms.size(), 0);
    }

    std::vector<int> matches(const std::vector<int>& counts) {
        next_stamp(seen_, stamp_);
        std::vector<int> result, reduced = counts;
        visit_deletions(reduced, deletions_, 0, [&](const std::vector<int>& signature) {
            const auto found = signatures_.find(signature);
            if (found == signatures_.end()) return;
            for (int id : found->second) {
                if (seen_[id] == stamp_) continue;
                seen_[id] = stamp_;
                const Histogram& match = histograms[id];
                if (config_.same_types && !same_type_set(counts, match.counts)) continue;
                result.insert(result.end(), match.neighborhoods.begin(), match.neighborhoods.end());
            }
        });
        std::sort(result.begin(), result.end());
        return result;
    }
};

void pack_occurrences(const Dataset& data, const Geometry& geometry, const Config& config,
                      Motif& motif, std::vector<std::uint32_t>& occupied,
                      std::uint32_t stamp) {
    std::set<int> samples;
    std::set<std::tuple<int, double, double>> regions;
    motif.raw_matches = static_cast<int>(motif.matches.size());
    for (int index : motif.matches) {
        const Neighborhood& group = geometry.neighborhoods[index];
        bool overlaps = false;
        for (int member : group.members) {
            if (occupied[member] == stamp) {
                overlaps = true;
                break;
            }
        }
        if (overlaps) continue;
        for (int member : group.members) occupied[member] = stamp;
        motif.occurrences.push_back(index);
        const Cell& center = data.cells[group.center];
        samples.insert(center.sample);
        regions.emplace(center.sample, std::floor(center.x / config.region_size),
                        std::floor(center.y / config.region_size));
    }
    motif.support = static_cast<int>(motif.occurrences.size());
    motif.samples = static_cast<int>(samples.size());
    motif.regions = static_cast<int>(regions.size());
}

void validate_scoring(const Dataset& data, const Geometry& geometry,
                      const std::vector<NeighborhoodStats>& stats, const Config& config) {
    if (stats.size() != geometry.neighborhoods.size())
        throw std::invalid_argument("Neighborhood statistics and geometry have different sizes");
    if (config.max_relabels < 0 || config.max_relabels > 2)
        throw std::invalid_argument("max_relabels must be 0, 1, or 2");
    if (!(config.region_size > 0) || !std::isfinite(config.region_size))
        throw std::invalid_argument("region_size must be finite and positive");
    for (const NeighborhoodStats& item : stats)
        if (item.counts.size() != data.types.size())
            throw std::invalid_argument("A neighborhood has the wrong number of cell types");
}

}  // namespace

Geometry build_geometry(const Dataset& data, const Config& config) {
    if (config.k < 1) throw std::invalid_argument("k must be positive");
    if (!(config.max_radius >= 0) || !std::isfinite(config.max_radius))
        throw std::invalid_argument("max_radius must be finite and nonnegative");
    if (!(config.contact_radius >= 0) || !std::isfinite(config.contact_radius))
        throw std::invalid_argument("contact_radius must be finite and nonnegative");
    if (data.cells.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("Too many cells for 32-bit cell indices");
    std::map<std::pair<int, int>, int> partition_ids;
    std::vector<std::vector<int>> partitions;
    std::vector<int> partition_of(data.cells.size());
    for (int index = 0; index < static_cast<int>(data.cells.size()); ++index) {
        const Cell& cell = data.cells[index];
        if (!std::isfinite(cell.x) || !std::isfinite(cell.y))
            throw std::invalid_argument("Cell coordinates must be finite");
        const auto key = std::make_pair(cell.sample, cell.component);
        const auto inserted = partition_ids.emplace(key, static_cast<int>(partitions.size()));
        if (inserted.second) partitions.emplace_back();
        const int partition = inserted.first->second;
        partitions[partition].push_back(index);
        partition_of[index] = partition;
    }
    std::vector<KDTree> trees;
    trees.reserve(partitions.size());
    for (const auto& partition : partitions) trees.emplace_back(data, partition);
    Geometry geometry;
    geometry.neighborhoods.reserve(data.cells.size());
    std::unordered_map<std::vector<int>, int, VectorHash> unique_sets;
    unique_sets.reserve(data.cells.size());
    for (int center = 0; center < static_cast<int>(data.cells.size()); ++center) {
        Neighborhood group;
        if (!trees[partition_of[center]].neighborhood(center, config.k, config.max_radius, group)) {
            ++geometry.rejected_radius;
            continue;
        }
        // Centers are visited by index: the first copy has the lowest center index.
        if (!unique_sets.emplace(group.members, center).second) {
            ++geometry.duplicate_sets;
            continue;
        }
        if (config.contact_radius > 0) {
            for (int a = 0; a < config.k; ++a) {
                const Cell& first = data.cells[group.members[a]];
                for (int b = a + 1; b < config.k; ++b) {
                    const Cell& second = data.cells[group.members[b]];
                    if (std::hypot(first.x - second.x, first.y - second.y) <= config.contact_radius)
                        group.edges.emplace_back(group.members[a], group.members[b]);
                }
            }
        }
        geometry.neighborhoods.push_back(std::move(group));
    }
    std::sort(geometry.neighborhoods.begin(), geometry.neighborhoods.end(),
              [](const Neighborhood& a, const Neighborhood& b) {
        return a.radius < b.radius || (a.radius == b.radius && a.center < b.center);
    });
    return geometry;
}

bool is_eligible_counts(const std::vector<int>& counts, const Config& config) {
    int total = 0, types = 0, largest = 0;
    for (int count : counts) {
        if (count < 0) return false;
        total += count;
        types += count > 0;
        largest = std::max(largest, count);
    }
    return total == config.k && total > 0 && types >= config.min_types &&
           static_cast<double>(largest) / total <= config.max_dominance + 1e-12;
}

int relabel_distance(const std::vector<int>& a, const std::vector<int>& b) {
    if (a.size() != b.size()) throw std::invalid_argument("Count vector sizes differ");
    int total = 0;
    for (std::size_t i = 0; i < a.size(); ++i) total += std::abs(a[i] - b[i]);
    return total / 2;
}

std::vector<NeighborhoodStats> characterize(const Dataset& data, const Geometry& geometry,
                                          const std::vector<int>& labels, const Config& config) {
    if (labels.size() != data.cells.size())
        throw std::invalid_argument("There must be exactly one label per cell");
    for (int label : labels)
        if (label < 0 || label >= static_cast<int>(data.types.size()))
            throw std::invalid_argument("Cell label is outside the type vocabulary");
    std::vector<NeighborhoodStats> result;
    result.reserve(geometry.neighborhoods.size());
    for (const Neighborhood& group : geometry.neighborhoods) {
        NeighborhoodStats stats;
        stats.counts.assign(data.types.size(), 0);
        for (int member : group.members) ++stats.counts[labels[member]];
        const double size = static_cast<double>(group.members.size());
        for (int count : stats.counts) {
            if (count == 0) continue;
            ++stats.n_types;
            const double fraction = count / size;
            stats.dominance = std::max(stats.dominance, fraction);
            stats.entropy -= fraction * std::log(fraction);
        }
        stats.effective_types = std::exp(stats.entropy);
        int cross_edges = 0;
        for (const auto& edge : group.edges) cross_edges += labels[edge.first] != labels[edge.second];
        stats.cross_edge_fraction = group.edges.empty() ? 0.0 :
            static_cast<double>(cross_edges) / group.edges.size();
        stats.eligible = is_eligible_counts(stats.counts, config) &&
                         stats.cross_edge_fraction + 1e-12 >= config.min_cross_edge_fraction;
        // An edgeless graph cannot meet any strictly positive edge-fraction requirement.
        if (group.edges.empty() && config.min_cross_edge_fraction > 0) stats.eligible = false;
        result.push_back(std::move(stats));
    }
    return result;
}

Motif score_template(const Dataset& data, const Geometry& geometry,
                     const std::vector<NeighborhoodStats>& stats, const std::vector<int>& counts,
                     const Config& config) {
    return score_templates(data, geometry, stats, {counts}, config).front();
}

std::vector<Motif> score_templates(const Dataset& data, const Geometry& geometry,
                                  const std::vector<NeighborhoodStats>& stats,
                                  const std::vector<std::vector<int>>& templates,
                                  const Config& config) {
    validate_scoring(data, geometry, stats, config);
    HistogramIndex index(stats, config);
    std::vector<Motif> result;
    result.reserve(templates.size());
    std::vector<std::uint32_t> occupied(data.cells.size(), 0);
    std::uint32_t stamp = 0;
    for (const std::vector<int>& counts : templates) {
        if (counts.size() != data.types.size())
            throw std::invalid_argument("Template has the wrong number of cell types");
        Motif motif;
        motif.counts = counts;
        if (is_eligible_counts(counts, config)) {
            motif.matches = index.matches(counts);
            next_stamp(occupied, stamp);
            pack_occurrences(data, geometry, config, motif, occupied, stamp);
        }
        result.push_back(std::move(motif));
    }
    return result;
}

std::vector<Motif> mine(const Dataset& data, const Geometry& geometry,
                        const std::vector<NeighborhoodStats>& stats, const Config& config) {
    validate_scoring(data, geometry, stats, config);
    HistogramIndex index(stats, config);
    std::vector<Motif> result;
    std::vector<std::uint32_t> occupied(data.cells.size(), 0);
    std::uint32_t occupied_stamp = 0;
    for (const Histogram& candidate : index.histograms) {
        Motif motif;
        motif.counts = candidate.counts;
        motif.matches = index.matches(candidate.counts);
        if (static_cast<int>(motif.matches.size()) < config.min_support) continue;
        next_stamp(occupied, occupied_stamp);
        pack_occurrences(data, geometry, config, motif, occupied, occupied_stamp);
        if (motif.support >= config.min_support) result.push_back(std::move(motif));
    }
    std::sort(result.begin(), result.end(), [](const Motif& a, const Motif& b) {
        if (a.support != b.support) return a.support > b.support;
        if (a.raw_matches != b.raw_matches) return a.raw_matches > b.raw_matches;
        return a.counts < b.counts;
    });
    // max_motifs is an output limit; the null search must retain the complete family.
    return result;
}

}  // namespace mixed
