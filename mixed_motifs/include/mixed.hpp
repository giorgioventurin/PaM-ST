#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace mixed {
struct Cell {
    std::string id;
    double x = 0, y = 0;
    int type = 0, sample = 0, stratum = 0, component = 0;
};
struct Dataset {
    std::vector<Cell> cells;
    std::vector<std::string> types, samples, strata, components;
};
struct Config {
    int k = 10, min_types = 2, max_relabels = 1;
    double max_radius = 100, max_dominance = 0.7;
    double contact_radius = 0, min_cross_edge_fraction = 0;
    double region_size = 1000;
    int min_support = 5, max_motifs = 20;
    int permutations = 199, threads = 1;
    std::uint64_t seed = 37;
    double alpha = 0.05;
    bool same_types = true;
};
struct Neighborhood {
    int center = 0;
    std::vector<int> members;
    std::vector<std::pair<int, int>> edges;
    double radius = 0;
};
struct Geometry {
    std::vector<Neighborhood> neighborhoods;
    int rejected_radius = 0, duplicate_sets = 0;
};
struct NeighborhoodStats {
    std::vector<int> counts;
    int n_types = 0;
    double dominance = 0, entropy = 0, effective_types = 0;
    double cross_edge_fraction = 0;
    bool eligible = false;
};
struct Motif {
    std::vector<int> counts;
    int raw_matches = 0, support = 0, regions = 0, samples = 0;
    std::vector<int> matches, occurrences;
};
Geometry build_geometry(const Dataset&, const Config&);
std::vector<NeighborhoodStats> characterize(const Dataset&, const Geometry&,
                                          const std::vector<int>& labels, const Config&);
std::vector<Motif> mine(const Dataset&, const Geometry&,
                        const std::vector<NeighborhoodStats>&, const Config&);
Motif score_template(const Dataset&, const Geometry&,
                     const std::vector<NeighborhoodStats>&, const std::vector<int>&,
                     const Config&);
std::vector<Motif> score_templates(const Dataset&, const Geometry&,
                                  const std::vector<NeighborhoodStats>&,
                                  const std::vector<std::vector<int>>& templates,
                                  const Config&);
bool is_eligible_counts(const std::vector<int>&, const Config&);
int relabel_distance(const std::vector<int>&, const std::vector<int>&);
}  // namespace mixed
