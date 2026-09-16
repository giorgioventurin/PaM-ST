#ifndef PAM_ST_SIGNIFICANCE_H
#define PAM_ST_SIGNIFICANCE_H

#include "structs.h"

#include <vector>

// Per-candidate enrichment and family-wise adjusted significance. Pure
// arithmetic over the null count matrix, so it is unit-tested directly.
struct FamilyScores {
    std::vector<double> null_mean;
    std::vector<double> null_sd;
    std::vector<double> lift;        // observed / null mean
    std::vector<double> p_raw;       // permutation p-value of each candidate
    std::vector<double> p_adjusted;  // Westfall-Young step-down, family-wise
    std::vector<int> order;          // candidates, most significant first
};

// null_counts is candidate-major: null_counts[c * permutations + b] is the count
// of candidate c in permutation b. first_index breaks ties reproducibly.
//
// The raw p-value of a candidate is (1 + #{draws >= observed}) / (B + 1). The
// adjustment asks how often the *minimum* p-value over the still-contending
// candidates would be as small by chance, which is what pays for having
// searched the whole family; it is then made monotone along the ranking.
FamilyScores score_family(const std::vector<int>& observed,
                          const std::vector<int>& first_index,
                          const std::vector<int>& null_counts,
                          int permutations);

// A candidate must be common enough to be testable and, optionally, mixed
// enough to be interesting. Both rules are fixed before any resampling.
bool passes_family_filters(const int* pattern,
                           int count,
                           int n_labels,
                           const AnalysisConfig& config);

#endif  // PAM_ST_SIGNIFICANCE_H
