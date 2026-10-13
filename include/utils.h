#ifndef PAM_ST_UTILS_H
#define PAM_ST_UTILS_H

#include "structs.h"

#include <algorithm>
#include <cstddef>
#include <string>
#include <thread>
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

Result most_frequent_pattern_test(const Dataset& data,
                                  const AnalysisConfig& config);

int count_matches(const Result& result);

#endif  // PAM_ST_UTILS_H
