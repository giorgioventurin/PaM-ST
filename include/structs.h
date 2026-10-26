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
    double tile_size = 1000.0;  // L: side of the square tiles labels are shuffled within
    int threads = 0;
    double alpha = 0.05;
    int min_recurrence = 10;  // kappa: fewest separate places a pattern recurs in
};

// A pattern of the observed tissue, with its test result.
struct Pattern {
    std::vector<int> counts;   // the count vector v, indexed by label
    int frequency = 0;         // f(v): neighbourhoods within rho of v
    int recurrence = 0;        // R(v): occurrences sharing no spot, min-degree greedy
    double z_score = 0.0;      // z of R(v) among its recurrences on all m+1 labellings
    double p_value = 1.0;      // share of labellings whose best z reaches it
    bool significant = false;  // p_value <= alpha

    // The matching neighbourhoods as they actually look: their mean size and
    // mean proportions per label, each neighbourhood weighted equally.
    double mean_cells = 0.0;
    std::vector<double> mean_composition;

    std::vector<int> occurrences;  // centre spots of occ(v); significant patterns only
};

struct Result {
    std::vector<Pattern> patterns;  // P^(0), highest z first
    int universe_size = 0;          // |U|: patterns of any labelling, up to scale
    int n_significant = 0;

    // Every labelling j = 0..m, the observed one first: how many patterns it
    // has, |P^(j)|, and its best score T^(j), -infinity when it has none.
    std::vector<int> labelling_patterns;
    std::vector<double> labelling_best;

    // Seconds spent in each pass over the m+1 labellings.
    double universe_seconds = 0.0;
    double recurrence_seconds = 0.0;
    double score_seconds = 0.0;
};

struct Timings {
    double load_seconds = 0.0;
    double analysis_seconds = 0.0;
    double output_seconds = 0.0;
    double total_seconds = 0.0;
};

using Clock = std::chrono::steady_clock;

#endif  // PAM_ST_STRUCTS_H
