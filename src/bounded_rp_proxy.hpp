#pragma once
#include "ils_sp/core.hpp"
#include <algorithm>
#include <cmath>
#include <optional>
#include <span>
#include <vector>
namespace ils_sp::detail {
// Work means one decoded/compared station ID, physical edge, charging-curve
// interval, catalogue preflight entry, or whole-solution addition. This is a
// deterministic engineering bound, not a CPU-time or arithmetic-op bound.
struct RpProxyWork {
  std::size_t cap{}, used{};
  bool take(std::size_t n = 1) noexcept {
    if (n > cap - used) return false;
    used += n;
    return true;
  }
};
// Select against the FINAL exact minimum, avoiding chained near-tie bias.
// The caller charges this bounded metadata pass before invoking it.
inline std::vector<std::size_t> rp_proxy_minimum_band(std::span<const double> scores) {
  std::vector<std::size_t> indices;
  if (scores.empty()) return indices;
  const double minimum = *std::min_element(scores.begin(), scores.end());
  for (std::size_t i = 0; i < scores.size(); ++i)
    if (scores[i] <= minimum + kCostTolerance) indices.push_back(i);
  return indices;
}
// Deliberately optimistic proposal-ranking heuristic, never a pruning bound.
// It ignores charger reachability/SOC, marginal-rate order and virtual SoC.
// Stream a single gap replacement without constructing or evaluating a Plan.
inline std::optional<double> rp_route_proxy(
    const Instance& instance, const Plan& plan, std::size_t replaced_gap,
    std::span<const int> replacement, double lambda, RpProxyWork& work) {
  const auto& vehicle = instance.vehicle();
  if (!std::isfinite(lambda) || lambda < 0 ||
      !std::isfinite(vehicle.speed_kmph) || vehicle.speed_kmph <= 0 ||
      !std::isfinite(vehicle.consumption_wh_per_km) || vehicle.consumption_wh_per_km <= 0 ||
      !std::isfinite(vehicle.battery_capacity_wh) || vehicle.battery_capacity_wh <= 0 ||
      !std::isfinite(vehicle.max_route_duration_h) || vehicle.max_route_duration_h < 0 ||
      plan.gaps.size() != plan.customer_ids.size() + 1) return std::nullopt;
  double distance = 0, segment_energy = 0, geometric_excess = 0;
  double service = 0, rate = 0;
  int previous = instance.depot().id;
  const auto append = [&](int id) {
    if (!work.take()) return false;
    const auto& node = instance.node(id);
    const double d = instance.distance_km(previous, id);
    if (!std::isfinite(d) || d < 0 || !std::isfinite(node.service_time_h) ||
        node.service_time_h < 0) return false;
    distance += d;
    segment_energy += d * vehicle.consumption_wh_per_km;
    service += node.service_time_h;
    if (node.kind == NodeKind::Station || node.kind == NodeKind::Depot) {
      geometric_excess += std::max(0.0, segment_energy - vehicle.battery_capacity_wh);
      segment_energy = 0;
    }
    if (node.kind == NodeKind::Station) {
      const auto& points = instance.curve_for_station(id).points();
      for (std::size_t j = 1; j < points.size(); ++j) {
        if (!work.take()) return false;
        const double dt = points[j].time_h - points[j - 1].time_h;
        const double dq = points[j].energy_wh - points[j - 1].energy_wh;
        if (!std::isfinite(dt) || dt <= 0 || !std::isfinite(dq) || dq < 0)
          return false;
        const double slope = dq / dt;
        if (!std::isfinite(slope)) return false;
        rate = std::max(rate, slope);
      }
    }
    previous = id;
    return true;
  };
  for (std::size_t gap = 0; gap < plan.gaps.size(); ++gap) {
    const std::span<const int> stations = gap == replaced_gap
        ? replacement : std::span<const int>(plan.gaps[gap]);
    for (int id : stations) if (!append(id)) return std::nullopt;
    if (gap < plan.customer_ids.size() && !append(plan.customer_ids[gap]))
      return std::nullopt;
  }
  if (!append(instance.depot().id)) return std::nullopt;
  const double energy = distance * vehicle.consumption_wh_per_km;
  const double charge = rate > 0
      ? std::max(0.0, energy - vehicle.battery_capacity_wh - geometric_excess) / rate : 0;
  const double raw = distance / vehicle.speed_kmph + charge;
  const double excess = std::max(0.0, raw + service - vehicle.max_route_duration_h);
  const double proxy = raw + lambda * (geometric_excess / vehicle.consumption_wh_per_km + excess);
  if (!std::isfinite(proxy) || !std::isfinite(energy) ||
      !std::isfinite(geometric_excess)) return std::nullopt;
  return proxy;
}
}  // namespace ils_sp::detail
