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

enum class DistanceMetric { Euclidean, JensenShannon };
enum class ErrorControl { Pointwise, FewRSFDR };
enum class StageErrorControl { None, Holm };
enum class NeighborhoodMode { Overlapping, Covering };
enum class NullModel { Global, Block };
// How candidate motifs are ranked and tested. Count is the original
// most-frequent-motif test; MinP ranks by family-wise adjusted significance.
enum class MotifStatistic { Count, MinP };

// Command-line spelling of each enum value, indexed by the enum's integer value.
template <typename Enum> struct EnumNames;
template <> struct EnumNames<DistanceMetric> {
    static constexpr std::array<const char*, 2> names{"l2", "js"};
};
template <> struct EnumNames<ErrorControl> {
    static constexpr std::array<const char*, 2> names{"pointwise", "fewrs-fdr"};
};
template <> struct EnumNames<StageErrorControl> {
    static constexpr std::array<const char*, 2> names{"none", "holm"};
};
template <> struct EnumNames<NeighborhoodMode> {
    static constexpr std::array<const char*, 2> names{"overlapping", "covering"};
};
template <> struct EnumNames<NullModel> {
    static constexpr std::array<const char*, 2> names{"global", "block"};
};
template <> struct EnumNames<MotifStatistic> {
    static constexpr std::array<const char*, 2> names{"count", "minp"};
};

struct AnalysisConfig {
    double radius = 300.0;
    double rho = 0.05;
    DistanceMetric metric = DistanceMetric::Euclidean;
    int permutations = 1000;
    int seed = 37;
    NullModel null_model = NullModel::Global;
    double block_size = 1000.0;
    double block_origin_x = 0.0;
    double block_origin_y = 0.0;
    int threads = 0;
    int max_motifs = 1;
    double alpha = 0.05;
    MotifStatistic statistic = MotifStatistic::Count;
    int min_support = 10;      // smallest observed count a candidate may have
    int min_types = 0;         // 0 disables the composition filter
    int min_type_cells = 1;    // cells a type needs to count towards min_types
    double split_size = 1000.0;  // checkerboard tile side for the selection/test split
    ErrorControl error_control = ErrorControl::Pointwise;
    double fdr_failure_probability = -1.0;
    int fdr_min_discoveries = 0;
    int fdr_order = 1;
    std::vector<std::string> frozen_cell_types;
    NeighborhoodMode neighborhood_mode = NeighborhoodMode::Overlapping;
    bool freeze_by_abundance = false;
    int max_freeze_stages = -1;
    StageErrorControl stage_error_control = StageErrorControl::Holm;
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

    // Significance search only; p_value holds the raw permutation p-value.
    double p_value_adjusted = 1.0;  // family-wise adjusted (Westfall-Young)
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
    int candidate_family_size = 0;  // candidates tested by the significance search
    int selection_center_count = 0; // centers that chose the family
    int test_center_count = 0;      // centers the family was tested on
    int split_discarded_count = 0;  // centers dropped at tile edges

    double fdr_bound = 1.0;
    int fdr_threshold = 0;
    int n_fdr_candidates = 0;
    bool fdr_certified = false;
    std::vector<int> fdr_null_order_statistics;

    int permutable_cell_count = 0;
    int null_block_count = 0;
    int null_shuffle_group_count = 0;
    int null_mixed_group_count = 0;
    int null_exchangeable_cell_count = 0;

    int candidate_center_count = 0;
    int covered_cell_count = 0;
    int overlapped_cell_count = 0;
    long long neighborhood_membership_count = 0;
    std::vector<int> candidate_center_ids;

    int freeze_stage = -1;
    int active_label_count = 0;
    double stage_p_value_raw = 1.0;
    double stage_p_value_adjusted = 1.0;
    bool stage_significant = false;

    const MotifTest& selected() const { return motif_tests[selected_rank - 1]; }
};

struct Timings {
    double load_seconds = 0.0;
    double analysis_seconds = 0.0;
    double output_seconds = 0.0;
    double total_seconds = 0.0;
};

struct FreezeStageResult {
    int stage_index = 0;
    std::vector<int> frozen_label_ids;
    Result result;
    Timings timings;
};

struct FreezeSweepResult {
    std::vector<int> freeze_order;
    std::vector<int> label_counts;
    std::vector<FreezeStageResult> stages;
    double neighbor_build_seconds = 0.0;
    double observed_histogram_seconds = 0.0;
};

using Clock = std::chrono::steady_clock;

#endif  // PAM_ST_STRUCTS_H
