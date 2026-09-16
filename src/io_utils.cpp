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

std::string pattern_as_string(const Dataset& data, const std::vector<int>& pattern) {
    std::ostringstream out;
    write_pattern_string(out, data, pattern);
    return out.str();
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

void write_csv_text(std::ostream& os, const std::string& value) {
    os << '"';
    for (const char c : value) {
        if (c == '"') os << '"';
        os << c;
    }
    os << '"';
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

void write_null_config(std::ostream& os, const char separator, const AnalysisConfig& config) {
    write_field(os, separator, "null_model", name_of(config.null_model));
    if (config.null_model != NullModel::Block) return;
    const auto flags = os.flags();
    const auto precision = os.precision();
    os << std::defaultfloat << std::setprecision(17);
    write_field(os, separator, "block_size", config.block_size);
    write_field(os, separator, "block_origin_x", config.block_origin_x);
    write_field(os, separator, "block_origin_y", config.block_origin_y);
    os.flags(flags);
    os.precision(precision);
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
    const bool has_frozen_type = !config.frozen_cell_types.empty();
    const int cell_count = static_cast<int>(data.labels.size());
    const double candidate_fraction = cell_count == 0
        ? 0.0
        : static_cast<double>(result.candidate_center_count) / cell_count;
    const double coverage_fraction = cell_count == 0
        ? 0.0
        : static_cast<double>(result.covered_cell_count) / cell_count;
    const double overlap_fraction = result.covered_cell_count == 0
        ? 0.0
        : static_cast<double>(result.overlapped_cell_count) / result.covered_cell_count;
    const double memberships_per_covered_cell = result.covered_cell_count == 0
        ? 0.0
        : static_cast<double>(result.neighborhood_membership_count) / result.covered_cell_count;
    const std::string frozen_names =
        has_frozen_type ? join_strings(config.frozen_cell_types) : "none";

    write_field(os, separator, "input", input);
    write_field(os, separator, "cells", cell_count);
    write_field(os, separator, "labels", data.label_names.size());
    write_field(os, separator, "frozen_cell_type", frozen_names);
    write_field(os, separator, "frozen_cell_types", frozen_names);
    write_field(os, separator, "frozen_type_count", config.frozen_cell_types.size());
    write_field(os, separator, "frozen_cells", cell_count - result.permutable_cell_count);
    write_field(os, separator, "permutable_cells", result.permutable_cell_count);
    write_field(os, separator, "neighborhood_mode", name_of(config.neighborhood_mode));
    write_field(os, separator, "candidate_centers", result.candidate_center_count);
    write_field(os, separator, "candidate_fraction", candidate_fraction);
    write_field(os, separator, "covered_cells", result.covered_cell_count);
    write_field(os, separator, "uncovered_cells", cell_count - result.covered_cell_count);
    write_field(os, separator, "coverage_fraction", coverage_fraction);
    write_field(os, separator, "overlapped_cells", result.overlapped_cell_count);
    write_field(os, separator, "overlap_fraction", overlap_fraction);
    write_field(os, separator, "neighborhood_memberships", result.neighborhood_membership_count);
    write_field(os, separator, "extra_overlap_memberships",
                result.neighborhood_membership_count - result.covered_cell_count);
    write_field(os, separator, "memberships_per_covered_cell", memberships_per_covered_cell);
    write_null_config(os, separator, config);
    write_field(os, separator, "null_blocks", result.null_block_count);
    write_field(os, separator, "null_shuffle_groups", result.null_shuffle_group_count);
    write_field(os, separator, "null_mixed_groups", result.null_mixed_group_count);
    write_field(os, separator, "null_exchangeable_cells", result.null_exchangeable_cell_count);
    if (config.null_model == NullModel::Block) {
        write_field(os, separator, "permutation_scheme", has_frozen_type
            ? "shuffle_non_frozen_labels_within_blocks" : "shuffle_labels_within_blocks");
        write_field(os, separator, "permutation_assumption", has_frozen_type
            ? "exchangeable_non_frozen_labels_within_blocks_conditional_on_frozen_locations_and_block_counts"
            : "exchangeable_labels_within_blocks_conditional_on_block_counts");
    } else {
        write_field(os, separator, "permutation_scheme",
                    has_frozen_type ? "shuffle_non_frozen_labels" : "shuffle_all_labels");
        write_field(os, separator, "permutation_assumption", has_frozen_type
            ? "exchangeable_non_frozen_labels_conditional_on_frozen_locations"
            : "exchangeable_labels");
    }
    write_field(os, separator, "radius", config.radius);
    write_field(os, separator, "rho", config.rho);
    write_field(os, separator, "distance_metric", name_of(config.metric));
    write_field(os, separator, "permutations", config.permutations);
    write_field(os, separator, "seed", config.seed);
    write_field(os, separator, "threads", config.threads);
    write_field(os, separator, "alpha", config.alpha);
    write_field(os, separator, "error_control", name_of(config.error_control));
    write_field(os, separator, "max_motifs", config.max_motifs);
    if (config.statistic == MotifStatistic::MinP) {
        write_field(os, separator, "statistic", name_of(config.statistic));
        write_field(os, separator, "candidate_family_size", result.candidate_family_size);
        write_field(os, separator, "min_support", config.min_support);
        write_field(os, separator, "min_types", config.min_types);
        write_field(os, separator, "min_type_cells", config.min_type_cells);
        write_field(os, separator, "split_size", config.split_size);
        write_field(os, separator, "selection_centers", result.selection_center_count);
        write_field(os, separator, "test_centers", result.test_center_count);
        write_field(os, separator, "split_edge_centers_dropped", result.split_discarded_count);
        write_field(os, separator, "disjoint_support", result.selected().disjoint_support);
        write_field(os, separator, "p_value_adjusted", result.selected().p_value_adjusted);
        write_field(os, separator, "lift", result.selected().lift);
        write_field(os, separator, "family_error_control", "westfall-young-fwer");
    }
    write_field(os, separator, "tested_motifs", result.tested_motifs);
    write_field(os, separator, "n_significant_motifs", result.n_significant_motifs);

    if (result.freeze_stage >= 0) {
        write_field(os, separator, "freeze_stage", result.freeze_stage);
        write_field(os, separator, "active_label_types", result.active_label_count);
        write_field(os, separator, "stage_error_control", name_of(config.stage_error_control));
        write_field(os, separator, "stage_p_value_raw", result.stage_p_value_raw);
        write_field(os, separator, "stage_p_value_adjusted", result.stage_p_value_adjusted);
        write_field(os, separator, "stage_significant", result.stage_significant);
    }

    if (config.error_control == ErrorControl::FewRSFDR) {
        write_field(os, separator, "fdr_failure_probability", config.fdr_failure_probability);
        write_field(os, separator, "fdr_fdp_cap", config.alpha - config.fdr_failure_probability);
        write_field(os, separator, "fdr_order", config.fdr_order);
        write_field(os, separator, "fdr_min_discoveries", config.fdr_min_discoveries);
        write_field(os, separator, "motif_separation", 2.0 * config.rho);
        write_field(os, separator, "fdr_threshold", result.fdr_threshold);
        write_field(os, separator, "n_fdr_candidates", result.n_fdr_candidates);
        write_field(os, separator, "fdr_bound", result.fdr_bound);
        write_field(os, separator, "fdr_certified", result.fdr_certified);
        write_field(os, separator, "fdr_assumption", "subset_pivotality");
    }

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

std::vector<std::string> frozen_names_for_stage(const Dataset& data,
                                                const FreezeStageResult& stage) {
    std::vector<std::string> names;
    names.reserve(stage.frozen_label_ids.size());
    for (const int label : stage.frozen_label_ids) names.push_back(data.label_names[label]);
    return names;
}

std::string stage_directory_name(const int stage) {
    std::ostringstream out;
    out << "stage_" << std::setw(2) << std::setfill('0') << stage;
    out << (stage == 0 ? "_none" : "_top_" + std::to_string(stage));
    return out.str();
}

// Per-stage configuration as it would have been given on the command line.
AnalysisConfig stage_config_of(const Dataset& data,
                               const FreezeStageResult& stage,
                               const AnalysisConfig& config) {
    AnalysisConfig stage_config = config;
    stage_config.freeze_by_abundance = false;
    stage_config.max_freeze_stages = -1;
    stage_config.frozen_cell_types = frozen_names_for_stage(data, stage);
    return stage_config;
}

}  // namespace

Dataset load_data(const std::string& path) {
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

    std::vector<std::array<double, 2>> coords;
    std::vector<std::string> raw_labels;
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

        coords.push_back({std::stod(fields[col_x]), std::stod(fields[col_y])});
        raw_labels.push_back(cell_type);
    }

    std::vector<std::string> names = raw_labels;
    std::sort(names.begin(), names.end());
    names.erase(std::unique(names.begin(), names.end()), names.end());

    std::unordered_map<std::string, int> label_of;
    label_of.reserve(names.size());
    for (int i = 0; i < static_cast<int>(names.size()); ++i) label_of[names[i]] = i;

    Dataset data;
    data.coords = std::move(coords);
    data.label_names = std::move(names);
    data.labels.resize(raw_labels.size());
    for (int i = 0; i < static_cast<int>(raw_labels.size()); ++i) {
        data.labels[i] = label_of.at(raw_labels[i]);
    }
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
    const bool fewrs = config.error_control == ErrorControl::FewRSFDR;

    // Drop files left by an earlier run in a mode that still wrote them.
    auto drop_stale = [&](const char* name) {
        std::filesystem::remove(std::filesystem::path(out_dir) / name);
    };
    if (!fewrs) {
        drop_stale("motif_patterns.csv");
        drop_stale("motif_matches.csv");
        drop_stale("fewrs_fdr_null_order.csv");
    }
    if (config.neighborhood_mode != NeighborhoodMode::Covering) {
        drop_stale("neighborhood_centers.csv");
    }

    {
        std::ofstream f = open_output_file(out_dir, "pattern.csv");
        f << "label,count\n";
        for (int i = 0; i < static_cast<int>(selected.pattern.size()); ++i) {
            f << data.label_names[i] << ',' << selected.pattern[i] << '\n';
        }
    }

    if (config.neighborhood_mode == NeighborhoodMode::Covering) {
        std::ofstream f = open_output_file(out_dir, "neighborhood_centers.csv");
        f << "center_index,cell_id,x,y,label\n";
        for (int i = 0; i < static_cast<int>(result.candidate_center_ids.size()); ++i) {
            const int cell = result.candidate_center_ids[i];
            f << i << ',' << cell << ','
              << data.coords[cell][0] << ',' << data.coords[cell][1] << ',';
            write_csv_text(f, data.label_names[data.labels[cell]]);
            f << '\n';
        }
    }

    {
        std::ofstream f = open_output_file(out_dir, "matches.csv");
        f << "cell_id,x,y,label\n";
        for (int i = 0; i < static_cast<int>(selected.match_mask.size()); ++i) {
            if (!selected.match_mask[i]) continue;
            f << i << ',' << data.coords[i][0] << ',' << data.coords[i][1] << ','
              << data.label_names[data.labels[i]] << '\n';
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

    if (config.statistic == MotifStatistic::MinP) {
        std::ofstream f = open_output_file(out_dir, "motif_significance.csv");
        f << "rank,observed_frequency,disjoint_support,null_mean,null_sd,lift,z_score,"
             "p_raw,p_adjusted,significant,family_size,pattern\n";
        f << std::setprecision(12);
        for (const MotifTest& test : result.motif_tests) {
            f << test.rank << ',' << test.observed << ',' << test.disjoint_support
              << ',' << test.null_mean << ',' << test.null_sd << ',' << test.lift << ',' << test.z_score << ','
              << test.p_value << ',' << test.p_value_adjusted << ','
              << (test.significant ? "true" : "false") << ','
              << result.candidate_family_size << ',';
            write_pattern_string(f, data, test.pattern);
            f << '\n';
        }
    }

    if (!fewrs) return;

    {
        std::ofstream f = open_output_file(out_dir, "motif_patterns.csv");
        f << "rank,observed_frequency,significant,label,count\n";
        for (const MotifTest& test : result.motif_tests) {
            for (int label = 0; label < static_cast<int>(test.pattern.size()); ++label) {
                f << test.rank << ','
                  << test.observed << ','
                  << (test.significant ? "true" : "false") << ','
                  << data.label_names[label] << ','
                  << test.pattern[label] << '\n';
            }
        }
    }

    {
        std::ofstream f = open_output_file(out_dir, "motif_matches.csv");
        f << "rank,cell_id,x,y,label\n";
        for (const MotifTest& test : result.motif_tests) {
            if (!test.significant) continue;
            for (int i = 0; i < static_cast<int>(test.match_mask.size()); ++i) {
                if (!test.match_mask[i]) continue;
                f << test.rank << ','
                  << i << ','
                  << data.coords[i][0] << ','
                  << data.coords[i][1] << ','
                  << data.label_names[data.labels[i]] << '\n';
            }
        }
    }

    {
        std::ofstream f = open_output_file(out_dir, "fewrs_fdr_null_order.csv");
        f << "resample,null_order,null_count\n";
        for (int i = 0; i < static_cast<int>(result.fdr_null_order_statistics.size()); ++i) {
            f << i + 1 << ',' << config.fdr_order << ','
              << result.fdr_null_order_statistics[i] << '\n';
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

void write_freeze_sweep_outputs(const std::string& out_dir,
                                const std::string& input,
                                const Dataset& data,
                                const FreezeSweepResult& sweep,
                                const AnalysisConfig& config) {
    if (out_dir.empty()) return;
    std::filesystem::create_directories(out_dir);

    const int max_frozen = sweep.stages.empty() ? 0 : sweep.stages.back().stage_index;
    {
        std::ofstream f = open_output_file(out_dir, "freeze_order.csv");
        f << "abundance_rank,label,count,fraction,first_frozen_stage,included_in_sweep\n";
        for (int rank = 0; rank < static_cast<int>(sweep.freeze_order.size()); ++rank) {
            const int label = sweep.freeze_order[rank];
            const double fraction = data.labels.empty()
                ? 0.0
                : static_cast<double>(sweep.label_counts[label]) / data.labels.size();
            f << rank + 1 << ',';
            write_csv_text(f, data.label_names[label]);
            f << ',' << sweep.label_counts[label]
              << ',' << std::setprecision(12) << fraction
              << ',' << rank + 1
              << ',' << std::boolalpha << (rank + 1 <= max_frozen) << '\n';
        }
    }

    {
        std::ofstream f = open_output_file(out_dir, "freeze_sweep_summary.csv");
        f << "stage,frozen_type_count,frozen_cell_types,frozen_cells,"
             "active_label_types,permutable_cells,candidate_centers,candidate_fraction,"
             "covered_cells,coverage_fraction,overlapped_cells,overlap_fraction,"
             "extra_overlap_memberships,"
             "observed_frequency,observed_fraction,null_mean,null_sd,"
             "observed_null_ratio,z_score,p_value_raw,p_value_adjusted,"
             "raw_significant,stage_significant,pattern,neighborhood_mode,"
             "analysis_seconds,null_model,block_size,block_origin_x,block_origin_y,"
             "null_blocks,null_shuffle_groups,null_mixed_groups,null_exchangeable_cells\n";
        f << std::boolalpha << std::setprecision(12);
        for (const FreezeStageResult& stage : sweep.stages) {
            const Result& result = stage.result;
            const MotifTest& primary = result.motif_tests.front();
            const std::vector<std::string> frozen_names = frozen_names_for_stage(data, stage);
            const double candidate_fraction = data.labels.empty()
                ? 0.0
                : static_cast<double>(result.candidate_center_count) / data.labels.size();
            const double observed_fraction = result.candidate_center_count == 0
                ? 0.0
                : static_cast<double>(primary.observed) / result.candidate_center_count;
            const double coverage_fraction = data.labels.empty()
                ? 0.0
                : static_cast<double>(result.covered_cell_count) / data.labels.size();
            const double overlap_fraction = result.covered_cell_count == 0
                ? 0.0
                : static_cast<double>(result.overlapped_cell_count) / result.covered_cell_count;
            const double observed_null_ratio = primary.null_mean == 0.0
                ? (primary.observed > 0 ? std::numeric_limits<double>::infinity() : 1.0)
                : primary.observed / primary.null_mean;

            f << stage.stage_index << ','
              << frozen_names.size() << ',';
            write_csv_text(f, frozen_names.empty() ? "none" : join_strings(frozen_names));
            f << ',' << static_cast<int>(data.labels.size()) - result.permutable_cell_count
              << ',' << result.active_label_count
              << ',' << result.permutable_cell_count
              << ',' << result.candidate_center_count
              << ',' << candidate_fraction
              << ',' << result.covered_cell_count
              << ',' << coverage_fraction
              << ',' << result.overlapped_cell_count
              << ',' << overlap_fraction
              << ',' << result.neighborhood_membership_count - result.covered_cell_count
              << ',' << primary.observed
              << ',' << observed_fraction
              << ',' << primary.null_mean
              << ',' << primary.null_sd
              << ',' << observed_null_ratio
              << ',' << primary.z_score
              << ',' << result.stage_p_value_raw
              << ',' << result.stage_p_value_adjusted
              << ',' << (primary.p_value <= config.alpha)
              << ',' << result.stage_significant << ',';
            write_csv_text(f, pattern_as_string(data, primary.pattern));
            f << ',' << name_of(config.neighborhood_mode)
              << ',' << stage.timings.analysis_seconds
              << ',' << name_of(config.null_model);
            if (config.null_model == NullModel::Block) {
                f << ',' << std::setprecision(17) << config.block_size
                  << ',' << config.block_origin_x << ',' << config.block_origin_y
                  << std::setprecision(12);
            } else {
                f << ",,,";
            }
            f << ',' << result.null_block_count
              << ',' << result.null_shuffle_group_count
              << ',' << result.null_mixed_group_count
              << ',' << result.null_exchangeable_cell_count << '\n';
        }
    }

    for (const FreezeStageResult& stage : sweep.stages) {
        const AnalysisConfig stage_config = stage_config_of(data, stage, config);
        const std::filesystem::path stage_dir =
            std::filesystem::path(out_dir) / stage_directory_name(stage.stage_index);
        write_outputs(stage_dir.string(), data, stage.result, stage_config);
        write_summary_csv(stage_dir.string(), input, data, stage.result, stage_config,
                          stage.timings);
        write_text_file((stage_dir / "report.txt").string(), [&](std::ostream& os) {
            write_report(os, input, data, stage.result, stage_config, stage.timings);
        });
    }
}

void write_freeze_sweep_report(std::ostream& os,
                               const std::string& input,
                               const Dataset& data,
                               const FreezeSweepResult& sweep,
                               const AnalysisConfig& config,
                               const Timings& timings) {
    os << std::boolalpha << std::fixed << std::setprecision(6);
    os << "input=" << input << '\n';
    os << "cells=" << data.labels.size() << '\n';
    os << "labels=" << data.label_names.size() << '\n';
    os << "freeze_by_abundance=true\n";
    os << "freeze_stages=" << sweep.stages.size() << '\n';
    write_null_config(os, '=', config);
    os << "neighborhood_mode=" << name_of(config.neighborhood_mode) << '\n';
    os << "stage_error_control=" << name_of(config.stage_error_control) << '\n';
    os << "alpha=" << config.alpha << '\n';
    os << "neighbor_build_seconds=" << sweep.neighbor_build_seconds << '\n';
    os << "observed_histogram_seconds=" << sweep.observed_histogram_seconds << '\n';
    os << "analysis_seconds=" << timings.analysis_seconds << '\n';
    os << "output_seconds=" << timings.output_seconds << '\n';
    os << "total_seconds=" << timings.total_seconds << '\n';
    os << "stage,frozen_types,candidate_centers,covered_cells,coverage_fraction,"
          "observed,null_mean,z_score,p_raw,p_adjusted,significant\n";
    for (const FreezeStageResult& stage : sweep.stages) {
        const MotifTest& primary = stage.result.motif_tests.front();
        const std::vector<std::string> names = frozen_names_for_stage(data, stage);
        os << stage.stage_index << ','
           << (names.empty() ? "none" : join_strings(names)) << ','
           << stage.result.candidate_center_count << ','
           << stage.result.covered_cell_count << ','
           << (data.labels.empty()
                ? 0.0
                : static_cast<double>(stage.result.covered_cell_count) / data.labels.size()) << ','
           << primary.observed << ','
           << primary.null_mean << ','
           << primary.z_score << ','
           << stage.result.stage_p_value_raw << ','
           << stage.result.stage_p_value_adjusted << ','
           << stage.result.stage_significant << '\n';
    }
}
