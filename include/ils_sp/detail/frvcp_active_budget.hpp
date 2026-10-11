#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>

namespace ils_sp::detail {

// One solve-local active-wall allowance. Idle search/SP/LP time is not charged.
// This is a cooperative engineering budget, not an optimality certificate.
struct FrvcpActiveWallBudget {
  using Clock = std::chrono::steady_clock;
  using Deadline = std::optional<Clock::time_point>;
  static constexpr double guard_s = 0.01;

  explicit FrvcpActiveWallBudget(std::optional<double> limit) : limit_s(limit) {
    if (limit_s && (!std::isfinite(*limit_s) || *limit_s <= 0.0)) {
      throw std::invalid_argument("FRVCP wall budget must be finite and positive");
    }
  }
  [[nodiscard]] bool enabled() const noexcept { return limit_s.has_value(); }
  [[nodiscard]] bool exhausted() const noexcept {
    return limit_s && *limit_s - spent_s <= guard_s;
  }
  [[nodiscard]] Deadline deadline(Clock::time_point now,
                                  Deadline global) const {
    if (!limit_s) return global;
    // Compare in floating tick space before any integral duration conversion.
    // This also handles huge finite inputs and negative test clock origins.
    const long double ticks = static_cast<long double>(
        std::max(0.0, *limit_s - spent_s)) * Clock::period::den / Clock::period::num;
    const long double target =
        static_cast<long double>(now.time_since_epoch().count()) + ticks;
    const auto maximum = Clock::time_point::max();
    const auto local = target >= static_cast<long double>(
        maximum.time_since_epoch().count()) ? maximum : Clock::time_point(
            Clock::duration(static_cast<Clock::duration::rep>(target)));
    return !global || local < *global ? Deadline(local) : global;
  }
  void charge(Clock::time_point start, Clock::time_point end) {
    if (limit_s) spent_s += std::chrono::duration<double>(end - start).count();
  }

  const std::optional<double> limit_s;
  double spent_s{}; // Keep observed overruns rather than clamping the receipt.
  std::size_t skipped_routes{};
};

}  // namespace ils_sp::detail
