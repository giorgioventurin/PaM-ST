#include "io_utils.h"
#include "utils.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
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

// The run in "name<separator>value" lines: inputs, parameters, what was found
// and how long it took. Shared by summary.csv and the text report.
void write_run_metrics(std::ostream& os,
                       const char separator,
                       const std::string& input,
                       const Dataset& data,
                       const Result& result,
                       const AnalysisConfig& config,
                       const Timings& timings) {
    write_field(os, separator, "input", input);
    write_field(os, separator, "cells", data.labels.size());
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
    write_field(os, separator, "radius", config.radius);
    write_field(os, separator, "rho", config.rho);
    write_field(os, separator, "min_recurrence", config.min_recurrence);
    write_field(os, separator, "tile_size", config.tile_size);
    write_field(os, separator, "permutations", config.permutations);
    write_field(os, separator, "seed", config.seed);
    write_field(os, separator, "alpha", config.alpha);
    write_field(os, separator, "threads", config.threads);
    write_field(os, separator, "observed_patterns", result.patterns.size());
    write_field(os, separator, "universe_size", result.universe_size);
    write_field(os, separator, "significant_patterns", result.n_significant);
    write_field(os, separator, "load_seconds", timings.load_seconds);
    write_field(os, separator, "analysis_seconds", timings.analysis_seconds);
    write_field(os, separator, "universe_seconds", result.universe_seconds);
    write_field(os, separator, "recurrence_seconds", result.recurrence_seconds);
    write_field(os, separator, "score_seconds", result.score_seconds);
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
    write_run_metrics(f, ',', input, data, result, config, timings);
}

void write_outputs(const std::string& out_dir, const Dataset& data, const Result& result) {
    if (out_dir.empty()) return;
    std::filesystem::create_directories(out_dir);

    {
        // Every pattern of the observed tissue, highest z first.
        std::ofstream f = open_output_file(out_dir, "patterns.csv");
        f << "rank,pattern,frequency,recurrence,z_score,p_value,significant,"
             "mean_cells,mean_composition\n";
        f << std::setprecision(12);
        for (int i = 0; i < static_cast<int>(result.patterns.size()); ++i) {
            const Pattern& pattern = result.patterns[i];
            f << i + 1 << ',';
            write_pattern_string(f, data, pattern.counts);
            f << ',' << pattern.frequency << ',' << pattern.recurrence << ','
              << pattern.z_score << ',' << pattern.p_value << ','
              << (pattern.significant ? "true" : "false") << ',' << pattern.mean_cells << ',';
            write_composition_string(f, data, pattern.mean_composition);
            f << '\n';
        }
    }

    {
        // Every labelling, the observed one (0) first: how many patterns it
        // has and its best score T, against which observed scores are judged.
        std::ofstream f = open_output_file(out_dir, "labellings.csv");
        f << "labelling,patterns,best_z\n";
        f << std::setprecision(12);
        for (int j = 0; j < static_cast<int>(result.labelling_best.size()); ++j) {
            f << j << ',' << result.labelling_patterns[j] << ',' << result.labelling_best[j] << '\n';
        }
    }

    {
        // The centre of every neighbourhood matching a significant pattern.
        // spot is the row within the sample's own input file, after filtering.
        std::ofstream f = open_output_file(out_dir, "occurrences.csv");
        f << "rank,sample,spot,x,y\n";
        f << std::setprecision(12);
        for (int i = 0; i < static_cast<int>(result.patterns.size()); ++i) {
            for (const int spot : result.patterns[i].occurrences) {
                f << i + 1 << ',' << data.sample_names[data.sample_of(spot)] << ','
                  << data.local_index(spot) << ',' << data.coords[spot][0] << ','
                  << data.coords[spot][1] << '\n';
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
    write_run_metrics(os, '=', input, data, result, config, timings);
}
