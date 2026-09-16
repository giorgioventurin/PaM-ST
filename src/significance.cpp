#include "significance.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

// Raw permutation p-values per candidate, then the Westfall-Young step-down
// adjustment on the minimum p-value over the family: a candidate is credited
// only if chance would beat it *anywhere* in the family that often.
//
// p-values are compared as exact integers. A raw p-value is k/(B+1) and a null
// draw's p-value is k/B, so p1 <= p2 becomes k1*(B+1) <= k2*B.
FamilyScores score_family(const std::vector<int>& observed,
                          const std::vector<int>& first_index,
                          const std::vector<int>& null_counts,
                          const int permutations) {
    const int n_candidates = static_cast<int>(observed.size());
    const long long b_total = permutations;
    FamilyScores scores;
    scores.null_mean.assign(n_candidates, 0.0);
    scores.null_sd.assign(n_candidates, 0.0);
    scores.lift.assign(n_candidates, 0.0);
    scores.p_raw.assign(n_candidates, 1.0);
    scores.p_adjusted.assign(n_candidates, 1.0);

    std::vector<long long> raw_numerator(n_candidates, 1);
    // p-value numerator of each draw against the other draws, per candidate.
    std::vector<std::vector<int>> draw_numerator(n_candidates);
    for (int c = 0; c < n_candidates; ++c) {
        const int* column = null_counts.data() + static_cast<std::size_t>(c) * permutations;
        double sum = 0.0;
        int at_least_observed = 0;
        for (int b = 0; b < permutations; ++b) {
            sum += column[b];
            if (column[b] >= observed[c]) ++at_least_observed;
        }
        scores.null_mean[c] = sum / permutations;
        double var = 0.0;
        for (int b = 0; b < permutations; ++b) {
            const double diff = column[b] - scores.null_mean[c];
            var += diff * diff;
        }
        scores.null_sd[c] = std::sqrt(var / permutations);
        scores.lift[c] = scores.null_mean[c] > 0.0
            ? observed[c] / scores.null_mean[c]
            : (observed[c] > 0 ? std::numeric_limits<double>::infinity() : 1.0);
        raw_numerator[c] = 1 + at_least_observed;
        scores.p_raw[c] = static_cast<double>(raw_numerator[c]) / (permutations + 1.0);

        std::vector<int> sorted(column, column + permutations);
        std::sort(sorted.begin(), sorted.end());
        draw_numerator[c].resize(permutations);
        for (int b = 0; b < permutations; ++b) {
            // #{draws >= this draw}, which already counts the draw itself.
            draw_numerator[c][b] = permutations - static_cast<int>(
                std::lower_bound(sorted.begin(), sorted.end(), column[b]) - sorted.begin());
        }
    }

    scores.order.resize(n_candidates);
    std::iota(scores.order.begin(), scores.order.end(), 0);
    std::sort(scores.order.begin(), scores.order.end(), [&](const int a, const int b) {
        if (raw_numerator[a] != raw_numerator[b]) return raw_numerator[a] < raw_numerator[b];
        if (scores.lift[a] != scores.lift[b]) return scores.lift[a] > scores.lift[b];
        return first_index[a] < first_index[b];
    });

    // Step down from the least significant candidate, shrinking the family the
    // minimum is taken over, then make the adjusted values monotone.
    std::vector<int> running_min(permutations, std::numeric_limits<int>::max());
    for (int rank = n_candidates - 1; rank >= 0; --rank) {
        const int c = scores.order[rank];
        int exceedances = 0;
        for (int b = 0; b < permutations; ++b) {
            running_min[b] = std::min(running_min[b], draw_numerator[c][b]);
            if (static_cast<long long>(running_min[b]) * (b_total + 1) <=
                raw_numerator[c] * b_total) {
                ++exceedances;
            }
        }
        scores.p_adjusted[c] = (1.0 + exceedances) / (permutations + 1.0);
    }
    for (int rank = 1; rank < n_candidates; ++rank) {
        scores.p_adjusted[scores.order[rank]] = std::max(
            scores.p_adjusted[scores.order[rank]], scores.p_adjusted[scores.order[rank - 1]]);
    }
    return scores;
}

// A candidate must be common enough to be testable and, optionally, mixed
// enough to be interesting. Both rules are fixed before any resampling, so the
// same family is scored in the observed data and in every permutation.
bool passes_family_filters(const int* pattern, const int count, const int n_labels,
                           const AnalysisConfig& config) {
    if (count < config.min_support) return false;
    if (config.min_types <= 0) return true;
    int types = 0;
    for (int d = 0; d < n_labels; ++d) {
        if (pattern[d] >= config.min_type_cells) ++types;
    }
    return types >= config.min_types;
}

