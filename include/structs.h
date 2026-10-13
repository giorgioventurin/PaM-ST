#ifndef PAM_ST_STRUCTS_H
#define PAM_ST_STRUCTS_H

#include <array>
#include <chrono>
#include <string>
#include <vector>

struct Dataset {
    std::vector<std::array<double, 2>> coords;
    std::vector<int> labels;
    std::vector<std::string> label_names;

    // Several tissues analysed together. Cells are stored sample by sample;
    // samples[i] is the sample of cell i, sample_starts[s] the first cell of
    // sample s, and sample_starts.back() the total.
    // Each sample keeps its own coordinate system. Empty vectors mean one sample.
    std::vector<int> samples;
    std::vector<std::string> sample_names;
    std::vector<int> sample_starts;

    int sample_count() const {
        return sample_names.empty() ? 1 : static_cast<int>(sample_names.size());
    }
    int sample_of(const int cell) const { return samples.empty() ? 0 : samples[cell]; }
    // Row of the cell within its own sample's input file (after filtering).
    int local_index(const int cell) const {
        return sample_starts.empty() ? cell : cell - sample_starts[sample_of(cell)];
    }
};

struct AnalysisConfig {
    double radius = 300.0;
    double rho = 0.05;
    int permutations = 1000;
    int seed = 37;
    double block_size = 1000.0;  // side of the square tiles labels are shuffled within
    int threads = 0;
    int max_motifs = 1;
    double alpha = 0.05;
    int min_support = 10;  // smallest number of separate places a pattern needs
    int min_types = 0;     // fewest labels a pattern holds; 0 disables the filter
};

// One motif's evidence within a single sample of a multi-sample run.
struct SampleEvidence {
    int observed = 0;
    int disjoint_support = 0;
    double null_mean = 0.0;
    double null_sd = 0.0;
    double lift = 0.0;
    double p_value = 1.0;  // permutation p-value within this sample, unadjusted
    double mean_cells = 0.0;             // mean cells inside the matching circles
    std::vector<double> mean_composition;  // their mean proportions, per label
};

struct MotifTest {
    int rank = 1;
    int observed = 0;
    std::vector<int> pattern;
    std::vector<int> null_counts;
    double p_value = 1.0;
    double null_mean = 0.0;
    double null_sd = 0.0;
    double z_score = 0.0;
    bool significant = false;
    std::vector<unsigned char> match_mask;

    // p_value is the pattern's own permutation rank; p_value_adjusted is the
    // max-search p-value, which accounts for every labelling's search.
    double p_value_adjusted = 1.0;
    double lift = 0.0;              // observed / null mean
    int disjoint_support = 0;       // occurrences sharing no cell

    // What the matching circles actually look like: their mean size and their
    // mean proportions per label, each circle weighted equally so that large
    // neighbourhoods do not dominate. composition_sd is the spread of the
    // per-sample means, and stays empty for a single sample.
    double mean_cells = 0.0;
    std::vector<double> mean_composition;
    std::vector<double> composition_sd;

    // Multi-sample runs only: the same evidence split by sample, and the number
    // of samples whose own p-value reaches alpha (descriptive, not adjusted).
    std::vector<SampleEvidence> per_sample;
    int replicated_in = 0;
};

struct Result {
    std::vector<MotifTest> motif_tests;
    int selected_rank = 1;  // 1-based rank of the reported motif in motif_tests
    int tested_motifs = 1;
    int n_significant_motifs = 0;
    int candidate_family_size = 0;  // patterns the observed tissue proposes
    int universe_size = 0;          // patterns proposed by any labelling

    int null_block_count = 0;
    int null_mixed_group_count = 0;
    int null_exchangeable_cell_count = 0;

    int candidate_center_count = 0;
    int covered_cell_count = 0;
    int overlapped_cell_count = 0;
    long long neighborhood_membership_count = 0;

    const MotifTest& selected() const { return motif_tests[selected_rank - 1]; }
};

struct Timings {
    double load_seconds = 0.0;
    double analysis_seconds = 0.0;
    double output_seconds = 0.0;
    double total_seconds = 0.0;
};

using Clock = std::chrono::steady_clock;

#endif  // PAM_ST_STRUCTS_H
