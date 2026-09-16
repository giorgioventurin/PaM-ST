#ifndef PAM_ST_NULL_MODEL_H
#define PAM_ST_NULL_MODEL_H

#include "structs.h"

#include <array>
#include <cstdint>
#include <random>
#include <vector>

// Nonempty groups of non-frozen cells whose labels are exchanged together.
// Global: one group in the supplied order. Block: one group per occupied tile,
// ordered by tile (x, y), each keeping the supplied cell order.
struct LabelShufflePlan {
    std::vector<std::vector<int>> groups;
    int block_count = 0;              // occupied tiles, including frozen-only tiles
    int mixed_group_count = 0;        // groups holding at least two distinct labels
    int exchangeable_cell_count = 0;  // cells in mixed groups
};

// Half-open tile containing a coordinate: floor((coordinate - origin) / size).
std::array<std::int64_t, 2> grid_tile(const std::array<double, 2>& coordinate,
                                      double size, double origin_x, double origin_y);

// Throws for invalid block-grid or analysis parameters.
void validate_null_model_config(const AnalysisConfig& config);

// Tiles are half-open squares with index floor((coordinate - origin) / size).
// The caller supplies every non-frozen cell, including cells that are not
// neighbourhood centers. Frozen cell coordinates are validated in block mode.
LabelShufflePlan make_label_shuffle_plan(
    const Dataset& data,
    const std::vector<int>& permutable_cells,
    const AnalysisConfig& config);

// Shuffles the current labels independently within each group, leaving all
// other entries fixed. Successive calls continue from the previous shuffle.
void shuffle_labels_within_groups(
    const LabelShufflePlan& plan,
    std::vector<int>& labels,
    std::mt19937& rng);

#endif  // PAM_ST_NULL_MODEL_H
