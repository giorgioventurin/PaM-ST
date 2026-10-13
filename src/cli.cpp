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

double parse_finite(const std::string& value, const std::string& option) {
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
        << "  --input FILE        Cell table; repeat to analyse several samples together\n"
        << "  --radius R          Spatial radius (default: 300)\n"
        << "  --rho RHO           Normalized-composition tolerance (default: 0.05)\n"
        << "  --permutations M    Number of label permutations (default: 1000)\n"
        << "  --seed S            RNG seed (default: 37)\n"
        << "  --block-size L      Side of the square tiles labels are shuffled within "
           "(default: 1000)\n"
        << "  --threads N         Worker threads (default: hardware)\n"
        << "  --alpha A           Significance level (default: 0.05)\n"
        << "  --min-support N     Smallest number of separate places a pattern needs "
           "(default: 10)\n"
        << "  --min-types M       Fewest labels a pattern holds (default: off)\n"
        << "  --output-dir DIR    Write the result CSVs\n"
        << "  --output-file FILE  Write the text report to FILE "
           "(default: DIR/report.txt when --output-dir is set)\n";
}

CliOptions parse_args(const int argc, char** argv) {
    CliOptions options;
    AnalysisConfig& config = options.analysis;
    config.threads = static_cast<int>(std::thread::hardware_concurrency());
    if (config.threads <= 0) config.threads = 1;

    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        auto need_value = [&](const std::string& option) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("Missing value for " + option);
            return argv[++i];
        };

        if (key == "--input") options.inputs.push_back(need_value(key));
        else if (key == "--output-dir") options.output_dir = need_value(key);
        else if (key == "--output-file") options.output_file = need_value(key);
        else if (key == "--radius") config.radius = std::stod(need_value(key));
        else if (key == "--rho") config.rho = std::stod(need_value(key));
        else if (key == "--permutations") config.permutations = std::stoi(need_value(key));
        else if (key == "--seed") config.seed = std::stoi(need_value(key));
        else if (key == "--block-size") config.block_size = parse_finite(need_value(key), key);
        else if (key == "--threads") config.threads = std::stoi(need_value(key));
        else if (key == "--alpha") config.alpha = std::stod(need_value(key));
        else if (key == "--min-support") config.min_support = std::stoi(need_value(key));
        else if (key == "--min-types") config.min_types = std::stoi(need_value(key));
        else if (key == "--help" || key == "-h") {
            usage(argv[0]);
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown argument: " + key);
        }
    }

    if (options.inputs.empty()) throw std::runtime_error("--input is required");
    for (std::size_t i = 0; i < options.inputs.size(); ++i) {
        options.input += (i ? ";" : "") + options.inputs[i];
    }
    validate_null_model_config(config);
    if (config.threads <= 0) config.threads = 1;
    if (config.alpha <= 0.0 || config.alpha >= 1.0) {
        throw std::runtime_error("alpha must be greater than 0 and less than 1");
    }
    if (config.min_support < 1) throw std::runtime_error("min-support must be at least 1");
    if (config.min_types < 0) throw std::runtime_error("min-types must be non-negative");
    if (config.permutations < 1) throw std::runtime_error("permutations must be at least 1");
    return options;
}
