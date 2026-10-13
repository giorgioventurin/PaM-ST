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

std::int64_t tile_index(const double coordinate, const double size) {
    if (!std::isfinite(coordinate)) {
        throw std::runtime_error("Block null requires finite cell coordinates");
    }
    const long double index =
        std::floor(static_cast<long double>(coordinate) / static_cast<long double>(size));
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

void validate_null_model_config(const AnalysisConfig& config) {
    if (!std::isfinite(config.block_size) || config.block_size <= 0.0) {
        throw std::runtime_error("block-size must be finite and greater than 0");
    }
    if (!std::isfinite(config.radius) || !std::isfinite(config.rho) ||
        !std::isfinite(config.alpha)) {
        throw std::runtime_error("Block null requires finite analysis parameters");
    }
}

LabelShufflePlan make_label_shuffle_plan(const Dataset& data, const AnalysisConfig& config) {
    validate_null_model_config(config);
    const std::size_t n = data.labels.size();
    if (n > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("Too many cells for the permutation engine");
    }
    if (data.coords.size() != n) {
        throw std::runtime_error("Coordinate and label counts differ");
    }

    // Tiles are keyed by sample as well, since each sample has its own
    // coordinate system and the same (x, y) tile means different places.
    using Tile = std::array<std::int64_t, 3>;
    std::map<Tile, std::vector<int>> tiles;
    for (std::size_t cell = 0; cell < n; ++cell) {
        const auto& coordinate = data.coords[cell];
        validate_neighbor_grid_range(coordinate, config.radius);
        const Tile tile{data.sample_of(static_cast<int>(cell)),
                        tile_index(coordinate[0], config.block_size),
                        tile_index(coordinate[1], config.block_size)};
        tiles[tile].push_back(static_cast<int>(cell));
    }

    LabelShufflePlan plan;
    for (auto& [tile, cells] : tiles) plan.groups.push_back(std::move(cells));
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
