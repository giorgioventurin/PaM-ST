// Differential checks deliberately use exhaustive neighbors and direct histogram scans.
// Checks throw rather than assert so they remain active in Release builds.
#include "mixed.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>

using namespace mixed;
namespace {
void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

Geometry exhaustive_geometry(const Dataset& data, const Config& config) {
    Geometry result;
    std::set<std::vector<int>> seen;
    for (int center = 0; center < static_cast<int>(data.cells.size()); ++center) {
        const auto& query = data.cells[center];
        std::vector<std::pair<double, int>> neighbors;
        for (int cell = 0; cell < static_cast<int>(data.cells.size()); ++cell) {
            const auto& other = data.cells[cell];
            if (cell == center || query.sample != other.sample || query.component != other.component) continue;
            const double distance = std::hypot(query.x - other.x, query.y - other.y);
            if (distance <= config.max_radius) neighbors.emplace_back(distance, cell);
        }
        std::sort(neighbors.begin(), neighbors.end());
        if (static_cast<int>(neighbors.size()) < config.k - 1) {
            ++result.rejected_radius;
            continue;
        }
        Neighborhood group;
        group.center = center;
        group.members.push_back(center);
        group.radius = config.k == 1 ? 0 : neighbors[config.k - 2].first;
        for (int index = 0; index < config.k - 1; ++index) group.members.push_back(neighbors[index].second);
        std::sort(group.members.begin(), group.members.end());
        if (!seen.insert(group.members).second) {
            ++result.duplicate_sets;
            continue;
        }
        if (config.contact_radius > 0) {
            for (int a = 0; a < config.k; ++a) for (int b = a + 1; b < config.k; ++b) {
                const auto& first = data.cells[group.members[a]];
                const auto& second = data.cells[group.members[b]];
                if (std::hypot(first.x - second.x, first.y - second.y) <= config.contact_radius)
                    group.edges.emplace_back(group.members[a], group.members[b]);
            }
        }
        result.neighborhoods.push_back(std::move(group));
    }
    std::sort(result.neighborhoods.begin(), result.neighborhoods.end(), [](const auto& a, const auto& b) {
        return std::tie(a.radius, a.center) < std::tie(b.radius, b.center);
    });
    return result;
}

bool expected_eligible(const std::vector<int>& counts, const Config& config) {
    int total = 0, types = 0, largest = 0;
    for (int count : counts) {
        if (count < 0) return false;
        total += count;
        if (count > 0) ++types;
        largest = std::max(largest, count);
    }
    return total == config.k && total > 0 && types >= config.min_types &&
           largest <= config.max_dominance * total + 1e-12;
}

void check_geometry(const Geometry& actual, const Geometry& expected) {
    require(actual.rejected_radius == expected.rejected_radius, "Radius rejection count differs");
    require(actual.duplicate_sets == expected.duplicate_sets, "Duplicate-set count differs");
    require(actual.neighborhoods.size() == expected.neighborhoods.size(), "Geometry size differs");
    for (std::size_t i = 0; i < actual.neighborhoods.size(); ++i) {
        const auto& a = actual.neighborhoods[i];
        const auto& e = expected.neighborhoods[i];
        require(a.center == e.center && a.members == e.members && a.radius == e.radius && a.edges == e.edges,
                "Geometry differs at neighborhood " + std::to_string(i));
    }
}

void check_characterization(const Geometry& geometry, const std::vector<int>& labels,
                            const std::vector<NeighborhoodStats>& stats, const Config& config,
                            std::size_t vocabulary_size) {
    require(stats.size() == geometry.neighborhoods.size(), "Statistics size differs");
    for (std::size_t i = 0; i < stats.size(); ++i) {
        const auto& group = geometry.neighborhoods[i];
        const auto& actual = stats[i];
        std::vector<int> counts(vocabulary_size, 0);
        for (int member : group.members) ++counts[labels[member]];
        int types = 0;
        double entropy = 0, dominance = 0;
        for (int count : counts) if (count) {
            ++types;
            const double fraction = count / static_cast<double>(config.k);
            entropy -= fraction * std::log(fraction);
            dominance = std::max(dominance, fraction);
        }
        int cross = 0;
        for (const auto& edge : group.edges) if (labels[edge.first] != labels[edge.second]) ++cross;
        const double fraction = group.edges.empty() ? 0 : cross / static_cast<double>(group.edges.size());
        const bool eligible = expected_eligible(counts, config) &&
            fraction >= config.min_cross_edge_fraction &&
            (!group.edges.empty() || config.min_cross_edge_fraction == 0);
        require(actual.counts == counts && actual.n_types == types && actual.eligible == eligible,
                "Composition eligibility differs");
        require(std::abs(actual.entropy - entropy) < 1e-12 &&
                std::abs(actual.dominance - dominance) < 1e-12 &&
                std::abs(actual.effective_types - std::exp(entropy)) < 1e-12 &&
                std::abs(actual.cross_edge_fraction - fraction) < 1e-12,
                "Composition/contact summary differs");
        require(is_eligible_counts(counts, config) == expected_eligible(counts, config),
                "Count eligibility differs from reference");
    }
}

Motif exhaustive_score(const Dataset& data, const Geometry& geometry,
                       const std::vector<NeighborhoodStats>& stats,
                       const std::vector<int>& counts, const Config& config) {
    Motif result;
    result.counts = counts;
    if (!expected_eligible(counts, config)) return result;
    std::set<int> occupied, samples;
    std::set<std::tuple<int, double, double>> regions;
    for (int index = 0; index < static_cast<int>(stats.size()); ++index) {
        if (!stats[index].eligible) continue;
        int distance = 0;
        bool same_types = true;
        for (std::size_t type = 0; type < counts.size(); ++type) {
            distance += std::abs(counts[type] - stats[index].counts[type]);
            if ((counts[type] > 0) != (stats[index].counts[type] > 0)) same_types = false;
        }
        if (distance > 2 * config.max_relabels || (config.same_types && !same_types)) continue;
        result.matches.push_back(index);
        const auto& group = geometry.neighborhoods[index];
        bool overlaps = false;
        for (int member : group.members) if (occupied.count(member)) overlaps = true;
        if (overlaps) continue;
        result.occurrences.push_back(index);
        for (int member : group.members) occupied.insert(member);
        const auto& center = data.cells[group.center];
        samples.insert(center.sample);
        regions.emplace(center.sample, std::floor(center.x / config.region_size),
                        std::floor(center.y / config.region_size));
    }
    result.raw_matches = static_cast<int>(result.matches.size());
    result.support = static_cast<int>(result.occurrences.size());
    result.samples = static_cast<int>(samples.size());
    result.regions = static_cast<int>(regions.size());
    return result;
}

void check_score(const Motif& actual, const Motif& expected) {
    require(actual.counts == expected.counts && actual.matches == expected.matches &&
            actual.occurrences == expected.occurrences && actual.raw_matches == expected.raw_matches &&
            actual.support == expected.support && actual.samples == expected.samples &&
            actual.regions == expected.regions, "Indexed score/packing differs from exhaustive reference");
}

void randomized_differential_checks() {
    std::mt19937 rng(40);
    for (int run = 0; run < 60; ++run) {
        Dataset data;
        data.types = {"A", "B", "C", "D", "E", "absent"};
        std::vector<int> labels;
        const int size = 50 + rng() % 200;
        for (int i = 0; i < size; ++i) {
            Cell cell;
            cell.x = static_cast<int>(rng() % 30) - 15;
            cell.y = static_cast<int>(rng() % 30) - 15;
            cell.sample = rng() % 2;
            cell.component = rng() % 2;
            cell.stratum = rng() % 2;
            cell.type = rng() % 5;
            data.cells.push_back(cell);
            labels.push_back(cell.type);
        }
        Config config;
        config.k = 1 + rng() % 10;
        config.max_radius = run % 2 ? 8 : 100;
        config.contact_radius = run % 3 ? 0 : 7;
        config.min_cross_edge_fraction = run % 3 ? 0 : .3;
        config.region_size = 7;
        config.min_types = run % 2 ? 2 : 3;
        config.max_dominance = run % 3 ? .7 : .5;
        config.max_motifs = 1;  // Mining must ignore the display limit.
        const auto geometry = build_geometry(data, config);
        check_geometry(geometry, exhaustive_geometry(data, config));
        const auto stats = characterize(data, geometry, labels, config);
        check_characterization(geometry, labels, stats, config, data.types.size());
        for (int tolerance = 0; tolerance <= 2; ++tolerance) for (bool same : {true, false}) {
            config.max_relabels = tolerance;
            config.same_types = same;
            config.min_support = 1 + run % 3;
            std::set<std::vector<int>> unique;
            for (const auto& item : stats) if (item.eligible) unique.insert(item.counts);
            std::vector<std::vector<int>> templates(unique.begin(), unique.end());
            templates.push_back({2, 2, 2, 2, 1, 1});  // Also test unobserved templates and absent labels.
            const auto scored = score_templates(data, geometry, stats, templates, config);
            require(scored.size() == templates.size(), "Batch scoring dropped templates");
            for (std::size_t i = 0; i < templates.size(); ++i)
                check_score(scored[i], exhaustive_score(data, geometry, stats, templates[i], config));
            std::vector<Motif> expected_mined;
            for (const auto& candidate : unique) {
                auto motif = exhaustive_score(data, geometry, stats, candidate, config);
                if (motif.support >= config.min_support) expected_mined.push_back(std::move(motif));
            }
            std::sort(expected_mined.begin(), expected_mined.end(), [](const Motif& a, const Motif& b) {
                if (a.support != b.support) return a.support > b.support;
                if (a.raw_matches != b.raw_matches) return a.raw_matches > b.raw_matches;
                return a.counts < b.counts;
            });
            const auto mined = mine(data, geometry, stats, config);
            require(mined.size() == expected_mined.size(), "Mining lost/added candidate templates");
            for (std::size_t i = 0; i < mined.size(); ++i) check_score(mined[i], expected_mined[i]);
        }
    }
}

void duplicate_coordinate_regression() {
    Dataset data;
    data.types = {"A", "B"};
    for (int i = 0; i < 20; ++i) {
        Cell cell;
        cell.type = i % 2;
        cell.stratum = i;  // Neighborhoods are allowed to cross strata.
        data.cells.push_back(cell);
    }
    Config config;
    config.k = 10;
    config.max_radius = 0;
    const auto geometry = build_geometry(data, config);
    check_geometry(geometry, exhaustive_geometry(data, config));
    require(geometry.neighborhoods.size() == 11 && geometry.duplicate_sets == 9,
            "Coincident points were not deterministically deduplicated");
    for (const auto& group : geometry.neighborhoods)
        require(std::binary_search(group.members.begin(), group.members.end(), group.center),
                "Center omitted when all neighbor distances tie");
}
}  // namespace

int main() {
    try {
        randomized_differential_checks();
        duplicate_coordinate_regression();
        std::cout << "Passed 60 exhaustive geometry comparisons, 360 scoring/mining comparisons, "
                     "and duplicate-coordinate regression.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Engine test failed: " << error.what() << '\n';
        return 1;
    }
}
