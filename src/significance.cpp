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
// Observed value and draws are one pooled set of B + 1 exchangeable values, and
// every p-value counts within that same set. Scoring the draws against the
// draws alone instead would put their smallest attainable p-value at 1/B while
// the observed reaches 1/(B + 1), so an observed value at its floor could never
// be matched and the adjustment would reject far too often. Sharing the
// denominator also makes every comparison an exact integer one.
//
// Note that the smallest attainable adjusted p-value is about
// family size / (B + 1): B has to grow with the family to keep any power.
FamilyScores score_family(const std::vector<int>& observed,
                          const std::vector<int>& first_index,
                          const std::vector<int>& null_counts,
                          const int permutations) {
    const int n_candidates = static_cast<int>(observed.size());
    FamilyScores scores;
    scores.null_mean.assign(n_candidates, 0.0);
    scores.null_sd.assign(n_candidates, 0.0);
    scores.lift.assign(n_candidates, 0.0);
    scores.p_raw.assign(n_candidates, 1.0);
    scores.p_adjusted.assign(n_candidates, 1.0);

    std::vector<int> raw_numerator(n_candidates, 1);
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
            // #{pooled values >= this draw}: the draws, which already count the
            // draw itself, plus the observed value.
            draw_numerator[c][b] = permutations - static_cast<int>(
                std::lower_bound(sorted.begin(), sorted.end(), column[b]) - sorted.begin()) +
                (observed[c] >= column[b] ? 1 : 0);
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
            if (running_min[b] <= raw_numerator[c]) ++exceedances;
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

