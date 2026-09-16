#include "significance.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(const bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

bool close(const double a, const double b) {
    return std::abs(a - b) < 1e-12;
}

// null_counts is candidate-major, so each row below is one candidate's draws.
std::vector<int> matrix(const std::vector<std::vector<int>>& rows) {
    std::vector<int> flat;
    for (const auto& row : rows) flat.insert(flat.end(), row.begin(), row.end());
    return flat;
}

void test_raw_p_values() {
    // Four draws. Candidate 0 is never matched, candidate 1 is beaten twice.
    const std::vector<int> observed = {10, 5};
    const std::vector<int> first_index = {0, 1};
    const auto counts = matrix({{1, 2, 3, 4}, {5, 6, 1, 2}});
    const FamilyScores scores = score_family(observed, first_index, counts, 4);

    require(close(scores.p_raw[0], 1.0 / 5.0), "unbeaten candidate must have p = 1/(B+1)");
    require(close(scores.p_raw[1], 3.0 / 5.0), "p counts draws greater than or equal to observed");
    require(close(scores.null_mean[0], 2.5), "null mean");
    require(close(scores.null_mean[1], 3.5), "null mean");
    require(close(scores.lift[0], 4.0), "lift is observed over null mean");
    require(scores.order[0] == 0, "the more extreme candidate ranks first");
}

void test_single_candidate_is_unadjusted() {
    // With one candidate there is nothing to correct for.
    const auto counts = matrix({{0, 1, 2, 7, 3, 1, 0, 2, 9}});
    const FamilyScores scores = score_family({5}, {0}, counts, 9);
    require(close(scores.p_adjusted[0], scores.p_raw[0]),
            "a family of one must leave the p-value unchanged");
}

void test_duplicate_candidates_cost_nothing() {
    // Identical candidates are the same test repeated; the family minimum is
    // the same as for one of them, so the adjustment must not inflate.
    const std::vector<int> draws = {4, 1, 0, 3, 2, 6, 1, 0};
    const FamilyScores one = score_family({5}, {0}, matrix({draws}), 8);
    const FamilyScores three =
        score_family({5, 5, 5}, {0, 1, 2}, matrix({draws, draws, draws}), 8);
    for (int c = 0; c < 3; ++c) {
        require(close(three.p_adjusted[c], one.p_adjusted[0]),
                "duplicated candidates must not change the adjusted p-value");
    }
}

// A draw that is extreme for one candidate has to be able to match an observed
// value that is extreme for another. Scoring draws against the draws alone puts
// their floor at 1/B while the observed reaches 1/(B+1), which would leave an
// observed floor unmatchable and reject far too often.
void test_extreme_draw_matches_extreme_observation() {
    const int B = 9;
    // Candidate 0: the observation beats every draw. Candidate 1: one draw does
    // the same within its own column.
    const std::vector<int> observed = {10, 0};
    const auto counts = matrix({{0, 0, 0, 0, 0, 0, 0, 0, 0},
                                {0, 0, 0, 99, 0, 0, 0, 0, 0}});
    const FamilyScores scores = score_family(observed, {0, 1}, counts, B);
    require(close(scores.p_raw[0], 1.0 / (B + 1.0)), "the observation is the most extreme value");
    require(scores.p_adjusted[0] > scores.p_raw[0] + 1e-12,
            "a rival candidate's extreme draw must cost the family something");
    require(close(scores.p_adjusted[0], 2.0 / (B + 1.0)), "exactly one draw matches it");
}

void test_adjustment_is_conservative_and_monotone() {
    // Twenty independent-looking candidates: adjusted >= raw everywhere, and
    // adjusted values never decrease along the ranking.
    std::vector<std::vector<int>> rows;
    std::vector<int> observed;
    std::vector<int> first_index;
    for (int c = 0; c < 20; ++c) {
        std::vector<int> draws;
        for (int b = 0; b < 50; ++b) draws.push_back((c * 7 + b * 13) % 11);
        rows.push_back(draws);
        observed.push_back(c % 5 + 6);
        first_index.push_back(c);
    }
    const FamilyScores scores = score_family(observed, first_index, matrix(rows), 50);
    for (int c = 0; c < 20; ++c) {
        require(scores.p_adjusted[c] >= scores.p_raw[c] - 1e-12,
                "an adjusted p-value can never be smaller than the raw one");
        require(scores.p_adjusted[c] <= 1.0 + 1e-12, "p-values stay within [0, 1]");
    }
    for (std::size_t rank = 1; rank < scores.order.size(); ++rank) {
        require(scores.p_adjusted[scores.order[rank]] >=
                    scores.p_adjusted[scores.order[rank - 1]] - 1e-12,
                "adjusted p-values must be monotone along the ranking");
        require(scores.p_raw[scores.order[rank]] >= scores.p_raw[scores.order[rank - 1]] - 1e-12,
                "the ranking must be by raw p-value");
    }
}

void test_family_filters() {
    AnalysisConfig config;
    config.min_support = 10;
    const std::vector<int> mixed = {4, 3, 0, 1};
    const std::vector<int> blob = {20, 1, 0, 0};

    require(!passes_family_filters(mixed.data(), 9, 4, config), "support floor excludes rare candidates");
    require(passes_family_filters(mixed.data(), 10, 4, config), "support floor is inclusive");
    require(passes_family_filters(blob.data(), 50, 4, config), "no composition filter by default");

    config.min_types = 3;
    config.min_type_cells = 2;
    require(!passes_family_filters(blob.data(), 50, 4, config),
            "a near-pure composition fails the multi-type filter");
    require(!passes_family_filters(mixed.data(), 50, 4, config),
            "types below min-type-cells do not count");
    config.min_type_cells = 1;
    require(passes_family_filters(mixed.data(), 50, 4, config),
            "three types with at least one cell pass");
}

}  // namespace

int main() {
    try {
        test_raw_p_values();
        test_single_candidate_is_unadjusted();
        test_duplicate_candidates_cost_nothing();
        test_extreme_draw_matches_extreme_observation();
        test_adjustment_is_conservative_and_monotone();
        test_family_filters();
        std::cout << "Significance scoring, Westfall-Young adjustment and family filters passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
