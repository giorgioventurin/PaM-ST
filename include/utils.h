#ifndef PAM_ST_UTILS_H
#define PAM_ST_UTILS_H

#include "structs.h"

#include <algorithm>
#include <atomic>
#include <cmath>
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
// A single chunk runs on the calling thread.
template <typename Fn>
void parallel_for(const int n, const int threads, Fn fn) {
    if (n <= 0) return;
    const int workers = normalize_thread_count(threads, n);
    if (workers == 1) {
        fn(0, n);
        return;
    }
    const int chunk = (n + workers - 1) / workers;
    std::vector<std::thread> pool;
    pool.reserve(workers);
    for (int lo = 0; lo < n; lo += chunk) {
        pool.emplace_back(fn, lo, std::min(n, lo + chunk));
    }
    for (std::thread& worker : pool) worker.join();
}

// Calls fn(i, worker) for every i in [0, n), handing the indices out one at a
// time in increasing order, so uneven items still keep every thread busy.
// worker is in [0, number of threads) and identifies the calling thread.
template <typename Fn>
void parallel_for_each(const int n, const int threads, Fn fn) {
    if (n <= 0) return;
    const int workers = normalize_thread_count(threads, n);
    std::atomic<int> next{0};
    auto run = [&](const int worker) {
        for (int i = next++; i < n; i = next++) fn(i, worker);
    };
    std::vector<std::thread> pool;
    pool.reserve(workers - 1);
    for (int worker = 1; worker < workers; ++worker) pool.emplace_back(run, worker);
    run(0);
    for (std::thread& thread : pool) thread.join();
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

// The terms m log m of the mixture entropy that one row fixes by itself: where
// the other row is zero, m = p / 2 (half), and where both rows are equal, m = p
// (self). Zero where p is zero.
void js_row_terms(const double* p, int k, double* half, double* self);

// jensen_shannon_divergence_from_entropy with the terms above looked up rather
// than recomputed. Each m is formed exactly as there and the sum runs in the
// same order, so the result is bit-identical; only the logarithms of mixtures
// of two different nonzero values are still evaluated.
inline double jensen_shannon_divergence_from_terms(const double* p,
                                                   const double* p_half,
                                                   const double* p_self,
                                                   const double hp,
                                                   const double* q,
                                                   const double* q_half,
                                                   const double hq,
                                                   const int k) {
    double hm = 0.0;
    for (int d = 0; d < k; ++d) {
        const double a = p[d];
        const double b = q[d];
        if (a == b) {
            if (a > 0.0) hm -= p_self[d];
        } else if (a == 0.0) {
            hm -= q_half[d];
        } else if (b == 0.0) {
            hm -= p_half[d];
        } else {
            const double m = 0.5 * (a + b);
            hm -= m * std::log(m);
        }
    }
    const double div = hm - 0.5 * (hp + hq);
    return div > 0.0 ? div : 0.0;
}

// An interval that certainly holds the value jensen_shannon_divergence_from_entropy
// returns for two rows (with their shannon_entropy), found without logarithms.
// Search trees compare it with their thresholds first and compute the
// divergence only when the interval straddles one, so they decide exactly as
// with the divergence itself.
struct DivergenceBounds {
    double lo = 0.0;
    double hi = 0.0;
};

// max(x, 0), exactly (2x and its half are exact, and -x + x is zero), without
// the branch compilers tend to emit for a comparison.
inline double positive_part(const double x) {
    return 0.5 * (std::fabs(x) + x);
}

// Widens bounds on the exact divergence into bounds on its rounded value:
// relative slack for rounding in the bounds themselves, absolute slack for the
// rounding of the divergence (below 1e-14 per component). Both exceed the
// errors they cover by several orders of magnitude.
inline DivergenceBounds widen_divergence_bounds(const double lower, const double upper,
                                                const int k) {
    const double slack = 1e-12 * (k + 16);
    return {positive_part(lower * (1.0 - 1e-9) - slack),
            positive_part(upper * (1.0 + 1e-9) + slack)};
}

// A cheaper, looser interval from the square roots of p and q (root_p[d] =
// sqrt(p[d])). Per component, with g and x as for jensen_shannon_bounds below,
// the divergence term over (sqrt(a) - sqrt(b))^2 is g(x) / (4 (1 - sqrt(1 - x^2))),
// a ratio of power series in x^2 whose coefficient ratios decrease; so it falls
// from 1/2 at x = 0 to ln(2) / 2 at |x| = 1, and the divergence lies between
// ln(2) / 2 and 1/2 times ||sqrt(p) - sqrt(q)||^2.
inline DivergenceBounds jensen_shannon_bounds_from_roots(const double* root_p,
                                                         const double* root_q,
                                                         const int k) {
    // Four partial sums let the additions overlap; their order is immaterial here.
    double sum[4] = {0.0, 0.0, 0.0, 0.0};
    int d = 0;
    for (; d + 4 <= k; d += 4) {
        for (int j = 0; j < 4; ++j) {
            const double diff = root_p[d + j] - root_q[d + j];
            sum[j] += diff * diff;
        }
    }
    for (; d < k; ++d) {
        const double diff = root_p[d] - root_q[d];
        sum[0] += diff * diff;
    }
    const double hellinger2 = (sum[0] + sum[1]) + (sum[2] + sum[3]);
    return widen_divergence_bounds(0.34657359027997264 * hellinger2, 0.5 * hellinger2, k);
}

// A tight interval from p and q themselves. Per component, with s = a + b and
// x = (a - b) / s, the divergence term is (s / 4) g(x), where
// g(x) = (1 + x) ln(1 + x) + (1 - x) ln(1 - x) = sum_{n >= 1} x^(2n) / (n (2n - 1)).
// All terms are positive, so three of them bound g from below, and adding x^8
// times the sum of the rest at |x| = 1 bounds it from above; where one side is
// zero, g is exactly 2 ln 2.
inline DivergenceBounds jensen_shannon_bounds(const double* p, const double* q, const int k) {
    constexpr double two_ln2 = 1.3862943611198906;
    // 2 ln 2 - 1 - 1/6 - 1/15 = 0.1529610277866..., rounded up.
    constexpr double tail = 0.15296103;
    // Written without branches: equal components give x = 0, components where
    // both sides are zero s = 0, and those where one side is zero x^2 = 1,
    // where g is exactly 2 ln 2.
    double lower = 0.0;  // sum of s g_lower(x)
    double upper = 0.0;  // sum of s g_upper(x)
    for (int d = 0; d < k; ++d) {
        const double a = p[d];
        const double b = q[d];
        const double s = a + b;
        // Arithmetic rather than branches, as the cases come unpredictably:
        // s + 1 where s is zero, and the exact 2 ln 2 where x^2 is one.
        const double x = (a - b) / (s + static_cast<double>(s == 0.0));
        const double x2 = x * x;
        const double x4 = x2 * x2;
        const double series = x2 * (1.0 + x2 * (1.0 / 6.0) + x4 * (1.0 / 15.0));
        const double one_sided = static_cast<double>(x2 == 1.0);
        lower += s * (series + one_sided * (two_ln2 - series));
        upper += s * (series + tail * (x4 * x4));
    }
    return widen_divergence_bounds(0.25 * lower, 0.25 * upper, k);
}

Result most_frequent_pattern_test(const Dataset& data,
                                  const AnalysisConfig& config);

FreezeSweepResult run_freeze_sweep(const Dataset& data,
                                   const AnalysisConfig& config);

int count_matches(const Result& result);

#endif  // PAM_ST_UTILS_H
