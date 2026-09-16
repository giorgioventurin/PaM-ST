#ifndef PAM_ST_UTILS_H
#define PAM_ST_UTILS_H

#include "structs.h"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

double elapsed_seconds(Clock::time_point start,
                       Clock::time_point end = Clock::now());

int normalize_thread_count(int requested, int n);

// Calls fn(lo, hi) on contiguous chunks covering [0, n), one chunk per thread.
template <typename Fn>
void parallel_for(const int n, const int threads, Fn fn) {
    if (n <= 0) return;
    const int workers = normalize_thread_count(threads, n);
    const int chunk = (n + workers - 1) / workers;
    std::vector<std::thread> pool;
    pool.reserve(workers);
    for (int lo = 0; lo < n; lo += chunk) {
        pool.emplace_back(fn, lo, std::min(n, lo + chunk));
    }
    for (std::thread& worker : pool) worker.join();
}

std::string trim(std::string s);
std::string lower_ascii(std::string s);

template <typename Enum>
std::string name_of(const Enum value) {
    return EnumNames<Enum>::names.at(static_cast<std::size_t>(value));
}

// Parses a case-insensitive enum name. `what` and `plural` word the error,
// e.g. "Unknown metric: cos. Available metrics: l2, js".
template <typename Enum>
Enum parse_enum(std::string value, const std::string& what, const std::string& plural) {
    value = lower_ascii(trim(std::move(value)));
    const auto& names = EnumNames<Enum>::names;
    std::string available;
    for (std::size_t i = 0; i < names.size(); ++i) {
        if (value == names[i]) return static_cast<Enum>(i);
        available += (i ? ", " : "") + std::string(names[i]);
    }
    throw std::runtime_error(
        "Unknown " + what + ": " + value + ". Available " + plural + ": " + available);
}

int fewrs_resample_count(double failure_probability);

double shannon_entropy(const double* p, int k);
double jensen_shannon_divergence_from_entropy(const double* p,
                                              double hp,
                                              const double* q,
                                              double hq,
                                              int k);

Result most_frequent_pattern_test(const Dataset& data,
                                  const AnalysisConfig& config);

FreezeSweepResult run_freeze_sweep(const Dataset& data,
                                   const AnalysisConfig& config);

int count_matches(const Result& result);

#endif  // PAM_ST_UTILS_H
