#include "io_utils.h"
#include "utils.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    bool in_quotes = false;
    for (const char c : line) {
        if (c == '"') {
            in_quotes = !in_quotes;
        } else if (c == ',' && !in_quotes) {
            fields.push_back(field);
            field.clear();
        } else {
            field.push_back(c);
        }
    }
    fields.push_back(field);
    return fields;
}

void write_pattern_string(std::ostream& os,
                          const Dataset& data,
                          const std::vector<int>& pattern) {
    bool first = true;
    for (int i = 0; i < static_cast<int>(pattern.size()); ++i) {
        if (pattern[i] == 0) continue;
        if (!first) os << ';';
        first = false;
        os << data.label_names[i] << ':' << pattern[i];
    }
}

// Mean proportions as percentages, largest first, in the same "label:value"
// shape as a pattern. Components below 0.05% are left out to keep it readable.
void write_composition_string(std::ostream& os,
                              const Dataset& data,
                              const std::vector<double>& composition) {
    std::vector<int> order;
    for (int i = 0; i < static_cast<int>(composition.size()); ++i) {
        if (composition[i] * 100.0 >= 0.05) order.push_back(i);
    }
    // Ties broken by label index, so the string is reproducible run to run.
    std::sort(order.begin(), order.end(), [&](const int a, const int b) {
        return composition[a] != composition[b] ? composition[a] > composition[b] : a < b;
    });
    std::ostringstream out;
    out << std::fixed << std::setprecision(2);
    for (int i = 0; i < static_cast<int>(order.size()); ++i) {
        if (i > 0) out << ';';
        out << data.label_names[order[i]] << ':' << composition[order[i]] * 100.0;
    }
    os << out.str();
}

std::string join_strings(const std::vector<std::string>& values,
                         const char* separator = ";") {
    std::ostringstream out;
    for (int i = 0; i < static_cast<int>(values.size()); ++i) {
        if (i > 0) out << separator;
        out << values[i];
    }
    return out.str();
}

std::ofstream open_output_file(const std::string& out_dir, const std::string& name) {
    std::ofstream file(std::filesystem::path(out_dir) / name);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot write " + name + " in: " + out_dir);
    }
    return file;
}

template <typename T>
void write_field(std::ostream& os, const char separator, const char* name, const T& value) {
    os << name << separator << value << '\n';
}

void write_run_metrics(std::ostream& os,
                       const char separator,
                       const bool include_pattern,
                       const std::string& input,
                       const Dataset& data,
                       const Result& result,
                       const AnalysisConfig& config,
                       const Timings& timings) {
    const MotifTest& selected = result.selected();
    const int cell_count = static_cast<int>(data.labels.size());
    const double coverage_fraction = cell_count == 0
        ? 0.0
        : static_cast<double>(result.covered_cell_count) / cell_count;
    const double overlap_fraction = result.covered_cell_count == 0
        ? 0.0
        : static_cast<double>(result.overlapped_cell_count) / result.covered_cell_count;
    const double memberships_per_covered_cell = result.covered_cell_count == 0
        ? 0.0
        : static_cast<double>(result.neighborhood_membership_count) / result.covered_cell_count;

    write_field(os, separator, "input", input);
    write_field(os, separator, "cells", cell_count);
    if (data.sample_count() > 1) {
        std::vector<std::string> sample_cells;
        for (int sample = 0; sample < data.sample_count(); ++sample) {
            sample_cells.push_back(std::to_string(
                data.sample_starts[sample + 1] - data.sample_starts[sample]));
        }
        write_field(os, separator, "samples", data.sample_count());
        write_field(os, separator, "sample_names", join_strings(data.sample_names));
        write_field(os, separator, "sample_cells", join_strings(sample_cells));
    }
    write_field(os, separator, "labels", data.label_names.size());
    write_field(os, separator, "candidate_centers", result.candidate_center_count);
    write_field(os, separator, "covered_cells", result.covered_cell_count);
    write_field(os, separator, "uncovered_cells", cell_count - result.covered_cell_count);
    write_field(os, separator, "coverage_fraction", coverage_fraction);
    write_field(os, separator, "overlapped_cells", result.overlapped_cell_count);
    write_field(os, separator, "overlap_fraction", overlap_fraction);
    write_field(os, separator, "neighborhood_memberships", result.neighborhood_membership_count);
    write_field(os, separator, "extra_overlap_memberships",
                result.neighborhood_membership_count - result.covered_cell_count);
    write_field(os, separator, "memberships_per_covered_cell", memberships_per_covered_cell);
    {
        const auto flags = os.flags();
        const auto precision = os.precision();
        os << std::defaultfloat << std::setprecision(17);
        write_field(os, separator, "block_size", config.block_size);
        os.flags(flags);
        os.precision(precision);
    }
    write_field(os, separator, "null_blocks", result.null_block_count);
    write_field(os, separator, "null_mixed_groups", result.null_mixed_group_count);
    write_field(os, separator, "null_exchangeable_cells", result.null_exchangeable_cell_count);
    write_field(os, separator, "radius", config.radius);
    write_field(os, separator, "rho", config.rho);
    write_field(os, separator, "permutations", config.permutations);
    write_field(os, separator, "seed", config.seed);
    write_field(os, separator, "threads", config.threads);
    write_field(os, separator, "alpha", config.alpha);
    write_field(os, separator, "max_motifs", config.max_motifs);
    write_field(os, separator, "candidate_family_size", result.candidate_family_size);
    write_field(os, separator, "min_support", config.min_support);
    write_field(os, separator, "min_types", config.min_types);
    write_field(os, separator, "pattern_universe_size", result.universe_size);
    write_field(os, separator, "test_centers", result.candidate_center_count);
    write_field(os, separator, "disjoint_support", result.selected().disjoint_support);
    write_field(os, separator, "p_value_adjusted", result.selected().p_value_adjusted);
    write_field(os, separator, "lift", result.selected().lift);
    write_field(os, separator, "tested_motifs", result.tested_motifs);
    write_field(os, separator, "n_significant_motifs", result.n_significant_motifs);

    write_field(os, separator, "selected_rank", result.selected_rank);
    write_field(os, separator, "significant", selected.significant);
    write_field(os, separator, "observed_frequency", selected.observed);
    write_field(os, separator, "p_value", selected.p_value);
    write_field(os, separator, "z_score", selected.z_score);
    write_field(os, separator, "null_mean", selected.null_mean);
    write_field(os, separator, "null_sd", selected.null_sd);
    write_field(os, separator, "n_matches", count_matches(result));
    if (include_pattern) {
        os << "pattern" << separator;
        write_pattern_string(os, data, selected.pattern);
        os << '\n';
    }
    write_field(os, separator, "load_seconds", timings.load_seconds);
    write_field(os, separator, "analysis_seconds", timings.analysis_seconds);
    write_field(os, separator, "output_seconds", timings.output_seconds);
    write_field(os, separator, "total_seconds", timings.total_seconds);
}

struct RawSample {
    std::vector<std::array<double, 2>> coords;
    std::vector<std::string> labels;
};

RawSample read_cells(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) throw std::runtime_error("Cannot open input CSV: " + path);

    std::string line;
    if (!std::getline(f, line)) throw std::runtime_error("Empty CSV: " + path);
    if (!line.empty() && line.back() == '\r') line.pop_back();

    const auto header = split_csv(line);
    int col_x = -1;
    int col_y = -1;
    int col_ct = -1;
    for (int i = 0; i < static_cast<int>(header.size()); ++i) {
        if (header[i] == "X_centroid") col_x = i;
        if (header[i] == "Y_centroid") col_y = i;
        if (header[i] == "Cell_Type") col_ct = i;
    }
    if (col_x < 0 || col_y < 0 || col_ct < 0) {
        throw std::runtime_error("CSV must contain X_centroid,Y_centroid,Cell_Type columns");
    }

    RawSample sample;
    const int required_col = std::max({col_x, col_y, col_ct});
    int row = 1;
    while (std::getline(f, line)) {
        ++row;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;

        const auto fields = split_csv(line);
        if (required_col >= static_cast<int>(fields.size())) {
            throw std::runtime_error(
                "Malformed row " + std::to_string(row) + ": not enough columns");
        }

        const std::string cell_type = trim(fields[col_ct]);
        if (lower_ascii(cell_type) == "unclassified") continue;

        sample.coords.push_back({std::stod(fields[col_x]), std::stod(fields[col_y])});
        sample.labels.push_back(cell_type);
    }
    return sample;
}

}  // namespace

Dataset load_samples(const std::vector<std::string>& paths) {
    std::vector<RawSample> raw;
    for (const std::string& path : paths) raw.push_back(read_cells(path));

    // One label index shared by every sample: the sorted union of their types.
    std::vector<std::string> names;
    for (const RawSample& sample : raw) {
        names.insert(names.end(), sample.labels.begin(), sample.labels.end());
    }
    std::sort(names.begin(), names.end());
    names.erase(std::unique(names.begin(), names.end()), names.end());
    std::unordered_map<std::string, int> label_of;
    label_of.reserve(names.size());
    for (int i = 0; i < static_cast<int>(names.size()); ++i) label_of[names[i]] = i;

    Dataset data;
    data.label_names = std::move(names);
    for (int s = 0; s < static_cast<int>(raw.size()); ++s) {
        data.sample_starts.push_back(static_cast<int>(data.labels.size()));
        std::string name = std::filesystem::path(paths[s]).stem().string();
        const std::string base_name = name;
        for (int copy = 2; std::find(data.sample_names.begin(), data.sample_names.end(), name) !=
                           data.sample_names.end(); ++copy) {
            name = base_name + "#" + std::to_string(copy);
        }
        data.sample_names.push_back(name);
        for (std::size_t i = 0; i < raw[s].labels.size(); ++i) {
            data.coords.push_back(raw[s].coords[i]);
            data.labels.push_back(label_of.at(raw[s].labels[i]));
            data.samples.push_back(s);
        }
    }
    data.sample_starts.push_back(static_cast<int>(data.labels.size()));
    return data;
}

void write_text_file(const std::string& path,
                     const std::function<void(std::ostream&)>& write) {
    if (path.empty()) return;
    const std::filesystem::path file(path);
    if (file.has_parent_path()) std::filesystem::create_directories(file.parent_path());
    std::ofstream out(file);
    if (!out.is_open()) throw std::runtime_error("Cannot write output file: " + path);
    write(out);
}

void write_summary_csv(const std::string& out_dir,
                       const std::string& input,
                       const Dataset& data,
                       const Result& result,
                       const AnalysisConfig& config,
                       const Timings& timings) {
    if (out_dir.empty()) return;
    std::filesystem::create_directories(out_dir);

    std::ofstream f = open_output_file(out_dir, "summary.csv");
    f << std::boolalpha << std::setprecision(12) << "metric,value\n";
    write_run_metrics(f, ',', false, input, data, result, config, timings);
}

void write_outputs(const std::string& out_dir,
                   const Dataset& data,
                   const Result& result,
                   const AnalysisConfig& config) {
    if (out_dir.empty()) return;
    std::filesystem::create_directories(out_dir);
    const MotifTest& selected = result.selected();
    const bool several_samples = data.sample_count() > 1;
    if (!several_samples) std::filesystem::remove(std::filesystem::path(out_dir) / "motif_samples.csv");

    {
        std::ofstream f = open_output_file(out_dir, "pattern.csv");
        f << "label,count\n";
        for (int i = 0; i < static_cast<int>(selected.pattern.size()); ++i) {
            f << data.label_names[i] << ',' << selected.pattern[i] << '\n';
        }
    }

    {
        std::ofstream f = open_output_file(out_dir, "matches.csv");
        // cell_id is the row within the cell's own input file; with several
        // samples the sample column says which file.
        f << "cell_id,x,y,label" << (several_samples ? ",sample" : "") << '\n';
        for (int i = 0; i < static_cast<int>(selected.match_mask.size()); ++i) {
            if (!selected.match_mask[i]) continue;
            f << data.local_index(i) << ',' << data.coords[i][0] << ',' << data.coords[i][1]
              << ',' << data.label_names[data.labels[i]];
            if (several_samples) f << ',' << data.sample_names[data.sample_of(i)];
            f << '\n';
        }
    }

    {
        std::ofstream f = open_output_file(out_dir, "null_max.csv");
        f << "permutation,null_max\n";
        for (int i = 0; i < static_cast<int>(selected.null_counts.size()); ++i) {
            f << i + 1 << ',' << selected.null_counts[i] << '\n';
        }
    }

    {
        std::ofstream f = open_output_file(out_dir, "motif_tests.csv");
        f << "rank,observed_frequency,p_value,z_score,null_mean,null_sd,significant,pattern\n";
        for (const MotifTest& test : result.motif_tests) {
            f << test.rank << ','
              << test.observed << ','
              << std::setprecision(12) << test.p_value << ','
              << std::setprecision(12) << test.z_score << ','
              << std::setprecision(12) << test.null_mean << ','
              << std::setprecision(12) << test.null_sd << ','
              << (test.significant ? "true" : "false") << ',';
            write_pattern_string(f, data, test.pattern);
            f << '\n';
        }
    }

    {
        std::ofstream f = open_output_file(out_dir, "null_rank_counts.csv");
        f << "permutation,rank,null_count\n";
        for (const MotifTest& test : result.motif_tests) {
            for (int i = 0; i < static_cast<int>(test.null_counts.size()); ++i) {
                f << i + 1 << ',' << test.rank << ',' << test.null_counts[i] << '\n';
            }
        }
    }

    {
        std::ofstream f = open_output_file(out_dir, "motif_significance.csv");
        f << "rank,observed_frequency,disjoint_support,null_mean,null_sd,lift,z_score,"
             "p_raw,p_adjusted,significant,"
          << (several_samples ? "replicated_in," : "")
          << "family_size,mean_cells,mean_composition,"
          << (several_samples ? "composition_sd," : "") << "pattern\n";
        f << std::setprecision(12);
        for (const MotifTest& test : result.motif_tests) {
            f << test.rank << ',' << test.observed << ',' << test.disjoint_support
              << ',' << test.null_mean << ',' << test.null_sd << ',' << test.lift << ',' << test.z_score << ','
              << test.p_value << ',' << test.p_value_adjusted << ','
              << (test.significant ? "true" : "false") << ',';
            if (several_samples) f << test.replicated_in << ',';
            f << result.candidate_family_size << ',' << test.mean_cells << ',';
            write_composition_string(f, data, test.mean_composition);
            f << ',';
            if (several_samples) {
                write_composition_string(f, data, test.composition_sd);
                f << ',';
            }
            write_pattern_string(f, data, test.pattern);
            f << '\n';
        }
    }

    if (several_samples) {
        // The same motifs sample by sample. Per-sample p-values are descriptive:
        // the family-wise guarantee is carried by p_adjusted in
        // motif_significance.csv, which pools all samples.
        std::ofstream f = open_output_file(out_dir, "motif_samples.csv");
        f << "rank,sample,observed_frequency,disjoint_support,null_mean,null_sd,lift,"
             "p_raw,p_raw_below_alpha,mean_cells,mean_composition\n";
        f << std::setprecision(12);
        for (const MotifTest& test : result.motif_tests) {
            for (int sample = 0; sample < static_cast<int>(test.per_sample.size()); ++sample) {
                const SampleEvidence& e = test.per_sample[sample];
                f << test.rank << ',' << data.sample_names[sample] << ',' << e.observed << ','
                  << e.disjoint_support << ',' << e.null_mean << ',' << e.null_sd << ','
                  << e.lift << ',' << e.p_value << ','
                  << (e.p_value <= config.alpha ? "true" : "false") << ','
                  << e.mean_cells << ',';
                write_composition_string(f, data, e.mean_composition);
                f << '\n';
            }
        }
    }
}

void write_report(std::ostream& os,
                  const std::string& input,
                  const Dataset& data,
                  const Result& result,
                  const AnalysisConfig& config,
                  const Timings& timings) {
    os << std::boolalpha << std::fixed << std::setprecision(6);
    write_run_metrics(os, '=', true, input, data, result, config, timings);
}
