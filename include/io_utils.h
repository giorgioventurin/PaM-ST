#ifndef PAM_ST_IO_UTILS_H
#define PAM_ST_IO_UTILS_H

#include "structs.h"

#include <functional>
#include <iosfwd>
#include <string>
#include <vector>

Dataset load_data(const std::string& path);

// Several tissues analysed together. Cells are concatenated sample by sample,
// each keeping its own coordinates; cell types are aligned by name.
Dataset load_samples(const std::vector<std::string>& paths);

// Creates parent directories of `path` and calls write on the opened file.
// Does nothing for an empty path.
void write_text_file(const std::string& path, const std::function<void(std::ostream&)>& write);

void write_summary_csv(const std::string& out_dir,
                       const std::string& input,
                       const Dataset& data,
                       const Result& result,
                       const AnalysisConfig& config,
                       const Timings& timings);

void write_outputs(const std::string& out_dir,
                   const Dataset& data,
                   const Result& result,
                   const AnalysisConfig& config);

void write_report(std::ostream& os,
                  const std::string& input,
                  const Dataset& data,
                  const Result& result,
                  const AnalysisConfig& config,
                  const Timings& timings);

void write_freeze_sweep_outputs(const std::string& out_dir,
                                const std::string& input,
                                const Dataset& data,
                                const FreezeSweepResult& sweep,
                                const AnalysisConfig& config);

void write_freeze_sweep_report(std::ostream& os,
                               const std::string& input,
                               const Dataset& data,
                               const FreezeSweepResult& sweep,
                               const AnalysisConfig& config,
                               const Timings& timings);

#endif  // PAM_ST_IO_UTILS_H
