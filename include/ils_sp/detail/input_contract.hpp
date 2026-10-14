#pragma once

#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>

namespace ils_sp::detail {

template <typename Integer>
Integer parse_nonnegative_integer(std::string_view value, std::string_view option) {
  static_assert(std::is_integral_v<Integer> && !std::is_same_v<Integer, bool>);
  std::uintmax_t parsed = 0;
  if (value.empty() || value.front() < '0' || value.front() > '9') {
    throw std::invalid_argument("invalid nonnegative integer for " + std::string(option));
  }
  const auto result = std::from_chars(value.data(), value.data() + value.size(), parsed);
  if (result.ec != std::errc{} || result.ptr != value.data() + value.size() ||
      parsed > static_cast<std::uintmax_t>(std::numeric_limits<Integer>::max())) {
    throw std::invalid_argument("integer out of range for " + std::string(option));
  }
  return static_cast<Integer>(parsed);
}

inline bool roundable_objective(double value) {
  const double cents = value * 100.0;
  // LLONG_MAX rounds upward when represented as double. Strictly below that
  // boundary is safely inside llround's nonnegative domain, including its .5 rule.
  return std::isfinite(value) && value >= 0.0 && std::isfinite(cents) &&
      cents < static_cast<double>(std::numeric_limits<long long>::max());
}

inline void validate_rounded_target(double value) {
  if (!roundable_objective(value)) {
    throw std::invalid_argument("rounded target is outside the finite rounding range");
  }
}

template <typename Clock>
[[nodiscard]] typename Clock::time_point checked_deadline_after(
    typename Clock::time_point start, double seconds) {
  using Duration = typename Clock::duration;
  using Rep = typename Duration::rep;
  static_assert(std::is_integral_v<Rep>);
  static_assert(std::numeric_limits<long double>::digits >=
                std::numeric_limits<Rep>::digits,
                "safe clock bounds require an exact long-double integer range");
  if (!std::isfinite(seconds) || seconds <= 0.0) {
    throw std::invalid_argument("time limit must be finite and positive");
  }
  // Compute in exactly the floating tick type used by duration_cast<double>,
  // then check before the floating-to-integral conversion and time-point sum.
  const double ticks = std::chrono::duration<double, typename Clock::period>(
      std::chrono::duration<double>(seconds)).count();
  const long double maximum = static_cast<long double>(std::numeric_limits<Rep>::max());
  const long double available = maximum -
      static_cast<long double>(start.time_since_epoch().count());
  if (!std::isfinite(ticks) || static_cast<long double>(ticks) > maximum ||
      static_cast<long double>(ticks) > available) {
    throw std::invalid_argument("time limit exceeds the clock representation");
  }
  return start + Duration(static_cast<Rep>(ticks));
}

inline void checked_flush(std::ostream& output) {
  output.flush();
  if (!output) throw std::runtime_error("output write failed");
}

}  // namespace ils_sp::detail
