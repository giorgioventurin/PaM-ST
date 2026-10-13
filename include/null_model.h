#ifndef PAM_ST_NULL_MODEL_H
#define PAM_ST_NULL_MODEL_H

#include "structs.h"

#include <random>
#include <vector>

// The tiles whose labels are exchanged together: one group per occupied tile,
// ordered by (sample, x, y), each keeping the cells in input order.
struct LabelShufflePlan {
    std::vector<std::vector<int>> groups;
    int mixed_group_count = 0;        // tiles holding at least two distinct labels
    int exchangeable_cell_count = 0;  // cells in those tiles
};

// Throws for invalid tile or analysis parameters.
void validate_null_model_config(const AnalysisConfig& config);

// Tiles are half-open squares of side block_size with index
// floor(coordinate / block_size), keyed by sample as well.
LabelShufflePlan make_label_shuffle_plan(const Dataset& data, const AnalysisConfig& config);

// Shuffles the current labels independently within each group, leaving all
// other entries fixed. Successive calls continue from the previous shuffle.
void shuffle_labels_within_groups(
    const LabelShufflePlan& plan,
    std::vector<int>& labels,
    std::mt19937& rng);

#endif  // PAM_ST_NULL_MODEL_H
