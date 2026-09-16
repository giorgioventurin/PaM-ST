#include "cli.h"
#include "null_model.h"
#include "utils.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>

namespace {

double parse_block_number(const std::string& value, const std::string& option) {
    std::size_t consumed = 0;
    double number;
    try {
        number = std::stod(value, &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error(option + " requires a finite number");
    }
    if (consumed != value.size() || !std::isfinite(number)) {
        throw std::runtime_error(option + " requires a finite number");
    }
    return number;
}

}  // namespace

void usage(const char* program) {
    std::cerr
        << "Usage: " << program << " --input cells.csv [options]\n"
        << "Options:\n"
        << "  --radius R          Spatial radius (default: 300)\n"
        << "  --rho RHO           Fixed normalized-composition tolerance (default: 0.05)\n"
        << "  --metric NAME       Distance metric: l2 or js (default: l2)\n"
        << "  --permutations B    Number of label permutations, or auto for FewRS (default: 1000)\n"
        << "  --seed S            RNG seed (default: 37)\n"
        << "  --null-model M      global or block (default: global)\n"
        << "  --block-size L      Square tile side in input units (default: 1000)\n"
        << "  --block-origin-x X  Grid x origin (default: 0)\n"
        << "  --block-origin-y Y  Grid y origin (default: 0)\n"
        << "                      Block options require --null-model block\n"
        << "  --threads N         Query threads for nonzero rho (default: hardware)\n"
        << "  --alpha A           Pointwise level or target FDR (default: 0.05)\n"
        << "  --max-motifs K      Motifs to report "
           "(frequency search: also the number tested; default: 1; alias: --k)\n"
        << "  --statistic S       count (most frequent motif, default) or minp "
           "(most significant)\n"
        << "  --min-support N     Smallest observed count a candidate may have "
           "(minp only; default: 10)\n"
        << "  --min-types M       Require M cell types in a candidate "
           "(minp only; default: off)\n"
        << "  --min-type-cells C  Cells a type needs to count towards --min-types "
           "(default: 1)\n"
        << "  --split-size L      Checkerboard tile side splitting selection from "
           "testing (minp only; default: 1000)\n"
        << "  --error-control M   pointwise or fewrs-fdr (default: pointwise)\n"
        << "  --fdr-failure-probability BETA\n"
        << "                      FewRS FDP-tail failure budget (default: alpha/2)\n"
        << "  --fdr-min-discoveries R\n"
        << "                      Minimum count required for FDR certification "
           "(default: max-motifs)\n"
        << "  --freeze-cell-type NAME\n"
        << "                      Freeze NAME in one run; repeat to freeze multiple types\n"
        << "  --neighborhood-mode M\n"
        << "                      overlapping or covering (default: overlapping)\n"
        << "  --freeze-by-abundance\n"
        << "                      Run cumulative stages freezing types by decreasing count\n"
        << "  --max-freeze-stages N\n"
        << "                      Freeze at most N types in the cumulative sweep\n"
        << "  --stage-error-control M\n"
        << "                      none or holm across sweep stages (default: holm)\n"
        << "  --output-dir DIR    Write summary/pattern/matches/null CSVs\n"
        << "  --output-file FILE  Write the text report to FILE "
           "(default: DIR/report.txt when --output-dir is set)\n";
}

CliOptions parse_args(const int argc, char** argv) {
    CliOptions options;
    AnalysisConfig& config = options.analysis;
    config.threads = static_cast<int>(std::thread::hardware_concurrency());
    if (config.threads <= 0) config.threads = 1;
    bool permutations_auto = false;
    bool max_freeze_stages_set = false;
    bool stage_error_control_set = false;
    bool block_option_set = false;
    bool min_support_set = false;
    bool family_filter_set = false;

    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        auto need_value = [&](const std::string& option) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("Missing value for " + option);
            return argv[++i];
        };

        if (key == "--input") options.input = need_value(key);
        else if (key == "--output-dir") options.output_dir = need_value(key);
        else if (key == "--output-file") options.output_file = need_value(key);
        else if (key == "--radius") config.radius = std::stod(need_value(key));
        else if (key == "--rho") config.rho = std::stod(need_value(key));
        else if (key == "--metric") {
            config.metric = parse_enum<DistanceMetric>(need_value(key), "metric", "metrics");
        }
        else if (key == "--permutations") {
            const std::string value = lower_ascii(trim(need_value(key)));
            permutations_auto = value == "auto";
            if (!permutations_auto) config.permutations = std::stoi(value);
        }
        else if (key == "--seed") config.seed = std::stoi(need_value(key));
        else if (key == "--null-model") {
            config.null_model = parse_enum<NullModel>(need_value(key), "null model", "models");
        }
        else if (key == "--block-size" || key == "--block-origin-x" ||
                 key == "--block-origin-y") {
            const double value = parse_block_number(need_value(key), key);
            if (key == "--block-size") config.block_size = value;
            else if (key == "--block-origin-x") config.block_origin_x = value;
            else config.block_origin_y = value;
            block_option_set = true;
        }
        else if (key == "--threads") config.threads = std::stoi(need_value(key));
        else if (key == "--alpha") config.alpha = std::stod(need_value(key));
        else if (key == "--max-motifs" || key == "--k") {
            config.max_motifs = std::stoi(need_value(key));
        }
        else if (key == "--statistic") {
            config.statistic =
                parse_enum<MotifStatistic>(need_value(key), "statistic", "statistics");
        }
        else if (key == "--min-support") {
            config.min_support = std::stoi(need_value(key));
            min_support_set = true;
        }
        else if (key == "--min-types") {
            config.min_types = std::stoi(need_value(key));
            family_filter_set = true;
        }
        else if (key == "--min-type-cells") {
            config.min_type_cells = std::stoi(need_value(key));
            family_filter_set = true;
        }
        else if (key == "--split-size") {
            config.split_size = parse_block_number(need_value(key), key);
            family_filter_set = true;
        }
        else if (key == "--error-control") {
            config.error_control =
                parse_enum<ErrorControl>(need_value(key), "error-control mode", "modes");
        }
        else if (key == "--fdr-failure-probability") {
            config.fdr_failure_probability = std::stod(need_value(key));
        }
        else if (key == "--fdr-min-discoveries") {
            config.fdr_min_discoveries = std::stoi(need_value(key));
        }
        else if (key == "--freeze-cell-type") {
            const std::string name = trim(need_value(key));
            if (name.empty()) throw std::runtime_error("freeze-cell-type must not be empty");
            config.frozen_cell_types.push_back(name);
        }
        else if (key == "--neighborhood-mode") {
            config.neighborhood_mode =
                parse_enum<NeighborhoodMode>(need_value(key), "neighborhood mode", "modes");
        }
        else if (key == "--freeze-by-abundance") config.freeze_by_abundance = true;
        else if (key == "--max-freeze-stages") {
            config.max_freeze_stages = std::stoi(need_value(key));
            max_freeze_stages_set = true;
        }
        else if (key == "--stage-error-control") {
            config.stage_error_control =
                parse_enum<StageErrorControl>(need_value(key), "stage-error-control mode", "modes");
            stage_error_control_set = true;
        }
        else if (key == "--help" || key == "-h") {
            usage(argv[0]);
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown argument: " + key);
        }
    }

    if (options.input.empty()) throw std::runtime_error("--input is required");
    if (block_option_set && config.null_model != NullModel::Block) {
        throw std::runtime_error("Block options require --null-model block");
    }
    validate_null_model_config(config);
    if (config.threads <= 0) config.threads = 1;
    if (config.alpha <= 0.0 || config.alpha >= 1.0) {
        throw std::runtime_error("alpha must be greater than 0 and less than 1");
    }
    if (config.max_motifs < 1) throw std::runtime_error("max-motifs must be at least 1");
    if (config.statistic == MotifStatistic::MinP) {
        if (config.min_support < 1) throw std::runtime_error("min-support must be at least 1");
        if (config.min_types < 0) throw std::runtime_error("min-types must be non-negative");
        if (config.min_type_cells < 1) {
            throw std::runtime_error("min-type-cells must be at least 1");
        }
        if (!(config.split_size > 0.0)) {
            throw std::runtime_error("split-size must be greater than 0");
        }
        if (config.split_size <= 2.0 * config.radius) {
            throw std::runtime_error(
                "split-size must exceed twice the radius, or no neighbourhood fits "
                "inside a tile");
        }
        if (config.freeze_by_abundance) {
            throw std::runtime_error(
                "--statistic minp does not support --freeze-by-abundance yet");
        }
        if (config.error_control != ErrorControl::Pointwise) {
            throw std::runtime_error(
                "--statistic minp controls family-wise error itself; use "
                "--error-control pointwise");
        }
    } else if (min_support_set || family_filter_set) {
        throw std::runtime_error(
            "--min-support, --min-types and --min-type-cells require --statistic minp");
    }
    if (config.freeze_by_abundance && !config.frozen_cell_types.empty()) {
        throw std::runtime_error(
            "--freeze-by-abundance and --freeze-cell-type are mutually exclusive");
    }
    if (max_freeze_stages_set && !config.freeze_by_abundance) {
        throw std::runtime_error("--max-freeze-stages requires --freeze-by-abundance");
    }
    if (stage_error_control_set && !config.freeze_by_abundance) {
        throw std::runtime_error("--stage-error-control requires --freeze-by-abundance");
    }
    if (max_freeze_stages_set && config.max_freeze_stages < 0) {
        throw std::runtime_error("max-freeze-stages must be non-negative");
    }
    if (config.freeze_by_abundance && config.error_control != ErrorControl::Pointwise) {
        throw std::runtime_error(
            "--freeze-by-abundance currently requires --error-control pointwise; "
            "stage-level Holm correction supplies whole-sweep error control");
    }

    if (config.error_control != ErrorControl::FewRSFDR) {
        if (permutations_auto) {
            throw std::runtime_error(
                "--permutations auto is only available with --error-control fewrs-fdr");
        }
        if (config.permutations < 1) throw std::runtime_error("permutations must be at least 1");
        return options;
    }

    if (config.fdr_failure_probability < 0.0) config.fdr_failure_probability = config.alpha / 2.0;
    if (config.fdr_failure_probability <= 0.0 ||
        config.fdr_failure_probability >= config.alpha) {
        throw std::runtime_error(
            "fdr-failure-probability must be greater than 0 and less than alpha");
    }
    if (config.fdr_min_discoveries == 0) config.fdr_min_discoveries = config.max_motifs;
    if (config.fdr_min_discoveries < 1 || config.fdr_min_discoveries > config.max_motifs) {
        throw std::runtime_error("fdr-min-discoveries must be between 1 and max-motifs");
    }

    const double fdp_cap = config.alpha - config.fdr_failure_probability;
    config.fdr_order = 1 + static_cast<int>(std::floor(fdp_cap * config.fdr_min_discoveries));
    const int required = fewrs_resample_count(config.fdr_failure_probability);
    if (permutations_auto) {
        config.permutations = required;
    } else if (config.permutations < required) {
        throw std::runtime_error(
            "fewrs-fdr requires at least " + std::to_string(required) +
            " permutations for the requested failure probability; use "
            "--permutations auto");
    }
    return options;
}
