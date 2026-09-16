#include "null_model.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

namespace {

void validate_grid_config(const AnalysisConfig& config) {
    if (!std::isfinite(config.block_size) || config.block_size <= 0.0) {
        throw std::runtime_error("block-size must be finite and greater than 0");
    }
    if (!std::isfinite(config.block_origin_x) ||
        !std::isfinite(config.block_origin_y)) {
        throw std::runtime_error("block origins must be finite");
    }
}

std::int64_t tile_index(const double coordinate, const double origin, const double size) {
    if (!std::isfinite(coordinate)) {
        throw std::runtime_error("Block null requires finite cell coordinates");
    }
    const long double index = std::floor(
        (static_cast<long double>(coordinate) - static_cast<long double>(origin)) /
        static_cast<long double>(size));
    // Exclusive power-of-two bound: INT64_MAX can round up to 2^63 as a float.
    const long double limit = std::ldexp(1.0L, 63);
    if (!std::isfinite(index) || index < -limit || index >= limit) {
        throw std::runtime_error("Block null tile index exceeds int64 range");
    }
    return static_cast<std::int64_t>(index);
}

// The neighbourhood grid divides coordinates by the radius in double precision
// and probes adjacent buckets; reject coordinates whose bucket would overflow.
void validate_neighbor_grid_range(const std::array<double, 2>& coordinate, const double radius) {
    if (radius <= 0.0) return;
    const double limit = std::ldexp(1.0, 63);
    for (const double value : coordinate) {
        const double bucket = std::floor(value / radius);
        if (!std::isfinite(bucket) || bucket <= -limit || bucket >= limit) {
            throw std::runtime_error("Cell coordinate exceeds safe neighborhood grid range");
        }
    }
}

}  // namespace

std::array<std::int64_t, 2> grid_tile(const std::array<double, 2>& coordinate,
                                      const double size, const double origin_x,
                                      const double origin_y) {
    return {tile_index(coordinate[0], origin_x, size), tile_index(coordinate[1], origin_y, size)};
}

void validate_null_model_config(const AnalysisConfig& config) {
    if (config.null_model != NullModel::Global && config.null_model != NullModel::Block) {
        throw std::runtime_error("Unknown null model");
    }
    if (config.null_model != NullModel::Block) return;
    validate_grid_config(config);
    if (!std::isfinite(config.radius) || !std::isfinite(config.rho) ||
        !std::isfinite(config.alpha) || !std::isfinite(config.fdr_failure_probability)) {
        throw std::runtime_error("Block null requires finite analysis parameters");
    }
}

LabelShufflePlan make_label_shuffle_plan(
    const Dataset& data,
    const std::vector<int>& permutable_cells,
    const AnalysisConfig& config) {
    validate_null_model_config(config);
    const std::size_t n = data.labels.size();
    if (n > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("Too many cells for the permutation engine");
    }
    std::vector<unsigned char> eligible(n, 0);
    for (const int cell : permutable_cells) {
        if (cell < 0 || static_cast<std::size_t>(cell) >= n) {
            throw std::runtime_error("Permutable cell index is out of range");
        }
        if (eligible[cell]) throw std::runtime_error("Duplicate permutable cell index");
        eligible[cell] = 1;
    }

    LabelShufflePlan plan;
    if (config.null_model == NullModel::Global) {
        if (!permutable_cells.empty()) plan.groups.push_back(permutable_cells);
    } else {
        if (data.coords.size() != n) {
            throw std::runtime_error("Coordinate and label counts differ");
        }
        using Tile = std::array<std::int64_t, 2>;
        std::map<Tile, std::vector<int>> tiles;
        std::vector<Tile> cell_tiles;
        cell_tiles.reserve(n);
        for (const auto& coordinate : data.coords) {
            cell_tiles.push_back(grid_tile(coordinate, config.block_size,
                                           config.block_origin_x, config.block_origin_y));
            tiles.try_emplace(cell_tiles.back());
            validate_neighbor_grid_range(coordinate, config.radius);
        }
        for (const int cell : permutable_cells) tiles[cell_tiles[cell]].push_back(cell);
        plan.block_count = static_cast<int>(tiles.size());
        for (auto& [tile, cells] : tiles) {
            if (!cells.empty()) plan.groups.push_back(std::move(cells));
        }
    }

    for (const auto& group : plan.groups) {
        const bool mixed = std::any_of(group.begin() + 1, group.end(),
            [&](const int cell) { return data.labels[cell] != data.labels[group.front()]; });
        if (mixed) {
            ++plan.mixed_group_count;
            plan.exchangeable_cell_count += static_cast<int>(group.size());
        }
    }
    return plan;
}

void shuffle_labels_within_groups(
    const LabelShufflePlan& plan,
    std::vector<int>& labels,
    std::mt19937& rng) {
    std::vector<int> shuffled;
    for (const auto& group : plan.groups) {
        shuffled.clear();
        for (const int cell : group) shuffled.push_back(labels[cell]);
        std::shuffle(shuffled.begin(), shuffled.end(), rng);
        for (std::size_t i = 0; i < group.size(); ++i) labels[group[i]] = shuffled[i];
    }
}
