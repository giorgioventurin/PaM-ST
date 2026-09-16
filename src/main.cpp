#include "cli.h"
#include "io_utils.h"
#include "utils.h"

#include <exception>
#include <filesystem>
#include <functional>
#include <iostream>
#include <ostream>
#include <string>

int main(const int argc, char** argv) {
    const auto program_start = Clock::now();
    CliOptions options;
    try {
        options = parse_args(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        usage(argv[0]);
        return 1;
    }

    try {
        const AnalysisConfig& config = options.analysis;
        std::string report_file = options.output_file;
        if (report_file.empty() && !options.output_dir.empty()) {
            report_file = (std::filesystem::path(options.output_dir) / "report.txt").string();
        }

        Timings timings;
        const auto load_start = Clock::now();
        const Dataset data = load_data(options.input);
        timings.load_seconds = elapsed_seconds(load_start);

        std::cerr << "Cells: " << data.labels.size() << '\n';
        std::cerr << "Labels: " << data.label_names.size() << '\n';

        // The text report goes to stdout and, when requested, to a file.
        auto emit = [&](const std::function<void(std::ostream&)>& report) {
            report(std::cout);
            write_text_file(report_file, report);
        };

        const auto analysis_start = Clock::now();
        if (config.freeze_by_abundance) {
            const FreezeSweepResult sweep = run_freeze_sweep(data, config);
            timings.analysis_seconds = elapsed_seconds(analysis_start);

            const auto output_start = Clock::now();
            write_freeze_sweep_outputs(options.output_dir, options.input, data, sweep, config);
            timings.output_seconds = elapsed_seconds(output_start);
            timings.total_seconds = elapsed_seconds(program_start);
            emit([&](std::ostream& os) {
                write_freeze_sweep_report(os, options.input, data, sweep, config, timings);
            });
            return 0;
        }

        const Result result = most_frequent_pattern_test(data, config);
        timings.analysis_seconds = elapsed_seconds(analysis_start);

        const auto output_start = Clock::now();
        write_outputs(options.output_dir, data, result, config);
        timings.output_seconds = elapsed_seconds(output_start);
        timings.total_seconds = elapsed_seconds(program_start);
        write_summary_csv(options.output_dir, options.input, data, result, config, timings);

        emit([&](std::ostream& os) {
            write_report(os, options.input, data, result, config, timings);
        });
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
