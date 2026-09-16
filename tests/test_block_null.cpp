#include "null_model.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

std::vector<std::vector<int>> canonical_groups(const LabelShufflePlan& plan) {
    auto groups = plan.groups;
    for (auto& group : groups) std::sort(group.begin(), group.end());
    std::sort(groups.begin(), groups.end());
    return groups;
}

std::map<int, int> histogram(const std::vector<int>& labels,
                             const std::vector<int>& cells) {
    std::map<int, int> counts;
    for (int cell : cells) ++counts[labels[cell]];
    return counts;
}

void test_grid_boundaries() {
    Dataset data;
    data.coords = {{-10.1, 0}, {-10, 0}, {-0.1, 0}, {0, 0},
                   {9.999, 0}, {10, 0}, {0, -0.1}, {0, 10}};
    data.labels = {0, 0, 1, 0, 1, 1, 0, 1};
    std::vector<int> cells(data.labels.size());
    std::iota(cells.begin(), cells.end(), 0);
    AnalysisConfig config;
    config.null_model = NullModel::Block;
    config.block_size = 10;
    const auto plan = make_label_shuffle_plan(data, cells, config);
    const std::vector<std::vector<int>> expected = {{0}, {1, 2}, {3, 4}, {5}, {6}, {7}};
    require(canonical_groups(plan) == expected,
            "Negative coordinates and exact x/y boundaries must use half-open tiles");

    data.coords = {{4.9, -5}, {5, -5}, {14.9, 4.9}, {15, 5}, {-5, -15}, {4.9, -5.1}};
    data.labels = {0, 1, 0, 1, 0, 1};
    cells.resize(data.labels.size());
    std::iota(cells.begin(), cells.end(), 0);
    config.block_origin_x = 5;
    config.block_origin_y = -5;
    const std::vector<std::vector<int>> shifted_expected = {{0}, {1, 2}, {3}, {4, 5}};
    require(canonical_groups(make_label_shuffle_plan(data, cells, config)) == shifted_expected,
            "Both grid origin coordinates must affect tile membership");
}

void test_conditional_shuffle_invariants() {
    Dataset data;
    data.coords = {{0, 0}, {1, 0}, {2, 0}, {3, 0}, {4, 0},
                   {10, 0}, {11, 0}, {12, 0}, {13, 0}, {14, 0},
                   {20, 0}, {21, 0}, {30, 0}, {40, 0}};
    data.labels = {0, 1, 0, 2, 9, 0, 2, 2, 1, 9, 3, 3, 4, 9};
    // Frozen cells 4 and 9 share mixed tiles, while 13 is a frozen-only tile.
    // All other cells must be shuffled, including cells not selected as centers.
    const std::vector<int> active = {0, 1, 2, 3, 5, 6, 7, 8, 10, 11, 12};
    const std::vector<std::vector<int>> tile_cells = {
        {0, 1, 2, 3, 4}, {5, 6, 7, 8, 9}, {10, 11}, {12}, {13}};
    AnalysisConfig config;
    config.null_model = NullModel::Block;
    config.block_size = 10;
    const auto plan = make_label_shuffle_plan(data, active, config);
    std::vector<int> planned;
    for (const auto& group : plan.groups) planned.insert(planned.end(), group.begin(), group.end());
    std::sort(planned.begin(), planned.end());
    require(planned == active, "Plan must contain every non-frozen cell exactly once");

    std::mt19937 rng(37), same_rng(37);
    auto labels = data.labels;
    auto repeated_labels = labels;
    std::vector<bool> changed(labels.size(), false);
    for (int draw = 0; draw < 300; ++draw) {
        shuffle_labels_within_groups(plan, labels, rng);
        shuffle_labels_within_groups(plan, repeated_labels, same_rng);
        require(labels == repeated_labels, "Equal seeds must reproduce every permutation");
        for (const auto& tile : tile_cells) {
            require(histogram(labels, tile) == histogram(data.labels, tile),
                    "Each tile must preserve its complete label histogram on every draw");
        }
        for (int cell : {4, 9, 13}) {
            require(labels[cell] == data.labels[cell], "Frozen cells must remain at their locations");
        }
        for (int cell : {10, 11, 12}) {
            require(labels[cell] == data.labels[cell], "Pure and singleton tiles cannot change labels");
        }
        for (std::size_t cell = 0; cell < labels.size(); ++cell) {
            changed[cell] = changed[cell] || labels[cell] != data.labels[cell];
        }
    }
    for (int cell : {0, 1, 2, 3, 5, 6, 7, 8}) {
        require(changed[cell], "Every non-frozen position in a mixed tile should participate");
    }
}

void test_one_tile_matches_legacy_shuffle() {
    Dataset data;
    for (int i = 0; i < 20; ++i) {
        data.coords.push_back({static_cast<double>(i), 0});
        data.labels.push_back(i % 4);
    }
    // Deliberately retain caller order rather than spatial or numeric index order.
    const std::vector<int> active = {13, 2, 17, 0, 9, 6, 19, 4, 7, 1, 12, 10};
    for (NullModel model : {NullModel::Global, NullModel::Block}) {
        AnalysisConfig config;
        config.null_model = model;
        config.block_size = 1000;
        config.block_origin_x = -100;
        config.block_origin_y = -100;
        const auto plan = make_label_shuffle_plan(data, active, config);
        require(plan.groups.size() == 1 && plan.groups.front() == active,
                "One tile must preserve legacy permutable-cell ordering");
        for (unsigned seed : {0U, 37U, 123456U}) {
            std::mt19937 rng(seed), legacy_rng(seed);
            auto labels = data.labels;
            auto expected = labels;
            std::vector<int> shuffled;
            for (int cell : active) shuffled.push_back(data.labels[cell]);
            for (int draw = 0; draw < 40; ++draw) {
                std::shuffle(shuffled.begin(), shuffled.end(), legacy_rng);
                for (std::size_t i = 0; i < active.size(); ++i) expected[active[i]] = shuffled[i];
                shuffle_labels_within_groups(plan, labels, rng);
                require(labels == expected,
                        "Global and one-tile block shuffles must exactly reproduce legacy RNG sequence");
            }
        }
    }
}

}  // namespace

int main() {
    try {
        test_grid_boundaries();
        test_conditional_shuffle_invariants();
        test_one_tile_matches_legacy_shuffle();
        std::cout << "Block null grid, conditional permutation, and legacy equivalence checks passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
