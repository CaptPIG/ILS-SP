#pragma once
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <utility>
#include <vector>
#if defined(__SSE2__)
#include <xmmintrin.h>
#endif
#include "ils_sp/core.hpp"

namespace ils_sp::detail {
// A certificate of this restricted 0 <= x <= 1, Ax = 1 model only.
struct SpDualCertificate {
  bool valid{};
  bool certified{};
  bool relative_gap_certified{};  // Observation only; never enables fast return.
  double upper_h{};
  double lower_h{};
  double gap_upper_h{};
  double relative_gap_upper{};
  double acceptance_gap_h{};
  bool fixing_valid{};
  std::size_t fixing_columns_considered{};
  double fixing_threshold_h{};
  std::vector<std::size_t> fixed_zero_columns;
};
inline bool sp_certificate_fp_supported() noexcept {
#if defined(__x86_64__) && defined(__SSE2__) && FLT_EVAL_METHOD == 0 && !defined(__FAST_MATH__) && (!defined(__FINITE_MATH_ONLY__) || __FINITE_MATH_ONLY__ == 0)
  if (!std::numeric_limits<double>::is_iec559 ||
      std::numeric_limits<double>::digits != 53) return false;
  unsigned short control;
  asm volatile("fnstcw %0" : "=m"(control));
  const unsigned csr = _mm_getcsr();
  return (control & 0x003fU) == 0x003fU &&
         (control & 0x0c00U) == 0U &&
         (csr & 0x1f80U) == 0x1f80U && (csr & 0xe040U) == 0U;
#else
  return false;
#endif
}
inline SpDualCertificate certify_sp_start(
    std::span<const double> costs,
    const std::vector<std::vector<std::size_t>>& incidence,
    std::span<const std::size_t> start,
    std::span<const double> duals, double abs_gap, double rel_gap,
    bool allow_reduced_cost_fixing = false) {
  SpDualCertificate result;
  if (!sp_certificate_fp_supported() || costs.empty() || incidence.empty() ||
      start.empty() || duals.size() != incidence.size() ||
      !std::isfinite(abs_gap) || abs_gap < 0.0 ||
      !std::isfinite(rel_gap) || rel_gap < 0.0) return result;
  bool finite = true;
  const auto down = [&](double rounded) {
    if (!std::isfinite(rounded)) { finite = false; return 0.0; }
    const double bound = std::nextafter(rounded, -std::numeric_limits<double>::infinity());
    if (!std::isfinite(bound)) finite = false;
    return bound;
  };
  const auto up = [&](double rounded) {
    if (!std::isfinite(rounded)) { finite = false; return 0.0; }
    const double bound = std::nextafter(rounded, std::numeric_limits<double>::infinity());
    if (!std::isfinite(bound)) finite = false;
    return bound;
  };
  std::vector<bool> selected(costs.size(), false);
  for (const auto index : start) {
    if (index >= costs.size() || selected[index]) return result;
    selected[index] = true;
  }
  double upper = 0.0;
  for (std::size_t j = 0; j < costs.size(); ++j) {
    if (!std::isfinite(costs[j])) return result;
    if (selected[j]) upper = up(upper + costs[j]);
  }
  double base = 0.0;
  std::vector<double> ay_upper(costs.size(), 0.0);
  std::vector<std::size_t> seen(costs.size(), incidence.size());
  for (std::size_t i = 0; i < incidence.size(); ++i) {
    if (!std::isfinite(duals[i]) || incidence[i].empty()) return result;
    base = down(base + duals[i]);
    std::size_t covered = 0;
    for (const auto j : incidence[i]) {
      if (j >= costs.size() || seen[j] == i) return result;
      seen[j] = i;
      covered += selected[j];
      ay_upper[j] = up(ay_upper[j] + duals[i]);
    }
    if (covered != 1) return result;
  }
  // For any finite y, sum(y) + sum(min(0,c-A^T y)) is a lower bound.
  double lower = base;
  std::vector<double> residual_lowers;
  if (allow_reduced_cost_fixing) residual_lowers.reserve(costs.size());
  for (std::size_t j = 0; j < costs.size(); ++j) {
    const double residual_lower = down(costs[j] - ay_upper[j]);
    lower = down(lower + std::min(0.0, residual_lower));
    if (allow_reduced_cost_fixing) residual_lowers.push_back(residual_lower);
  }
  if (!finite || lower > upper) return result;
  const double gap = up(upper - lower);
  if (!finite || gap < 0.0) return result;
  double relative = 0.0;
  bool relative_ok = false;
  if (lower >= 0.0 && upper > 0.0) {
    // Coupled monotonic bound: 1-L/U increases with U when L >= 0.
    relative = up(gap / upper);
    if (!finite) return result;
    relative_ok = relative <= rel_gap;
  }
  result.valid = true;
  result.upper_h = upper;
  result.lower_h = lower;
  result.gap_upper_h = gap;
  result.relative_gap_upper = relative;
  result.relative_gap_certified = relative_ok;
  result.acceptance_gap_h = std::min(abs_gap, kCostTolerance);
  result.certified = gap <= result.acceptance_gap_h;
  if (allow_reduced_cost_fixing && !result.certified) {
    // For x_j=1, replace min(0,r_j) in the original dual bound by r_j.
    // Both terms are conservative lower bounds. A strict positive margin
    // preserves every complete feasible solution no worse than this start.
    const double threshold = up(upper + kCostTolerance);
    if (!finite) return SpDualCertificate{};
    std::vector<std::size_t> fixed;
    for (std::size_t j = 0; j < costs.size(); ++j) {
      const double conditional_lower = down(lower + std::max(0.0, residual_lowers[j]));
      if (!finite) return SpDualCertificate{};
      if (conditional_lower > threshold) {
        if (selected[j]) return SpDualCertificate{};
        fixed.push_back(j);
      }
    }
    result.fixing_valid = true;
    result.fixing_columns_considered = costs.size();
    result.fixing_threshold_h = threshold;
    result.fixed_zero_columns = std::move(fixed);
  }
  return result;
}
}  // namespace ils_sp::detail
