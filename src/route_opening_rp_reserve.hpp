#pragma once
#include <cstddef>
#include <utility>
namespace ils_sp::detail {
// Own one of the existing lifetime tokens. Only the sampler debits remaining;
// releasing ownership never creates or replenishes a proposal token.
struct RouteOpeningRpReserve {
  bool held{};
  void initialize(std::size_t total) noexcept { held = total != 0; }
  std::size_t ordinary_allowance(std::size_t remaining) const noexcept {
    return remaining - static_cast<std::size_t>(held && remaining != 0);
  }
  bool release() noexcept { return std::exchange(held, false); }
  bool release_on_opening(std::size_t before, std::size_t after) noexcept {
    return after > before && release();
  }
};
}  // namespace ils_sp::detail
