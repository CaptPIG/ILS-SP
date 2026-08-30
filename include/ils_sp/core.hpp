#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace ils_sp {

constexpr double kEnergyToleranceWh = 1e-6;
constexpr double kTimeToleranceH = 1e-9;
constexpr double kCostTolerance = 1e-9;

template <typename T>
inline void hash_combine(std::size_t& seed, const T& value) noexcept {
  seed ^= std::hash<T>{}(value) + 0x9e3779b97f4a7c15ULL + (seed << 6U) +
          (seed >> 2U);
}

struct IntVectorHash {
  std::size_t operator()(const std::vector<int>& values) const noexcept {
    std::size_t seed = values.size();
    for (const int value : values) {
      hash_combine(seed, value);
    }
    return seed;
  }
};

enum class NodeKind : int { Depot = 0, Customer = 1, Station = 2 };

struct ChargingPoint {
  double time_h{};
  double energy_wh{};
};

class ChargingCurve {
 public:
  ChargingCurve() = default;
  ChargingCurve(std::string station_type, std::vector<ChargingPoint> points);

  [[nodiscard]] const std::string& station_type() const noexcept {
    return station_type_;
  }
  [[nodiscard]] const std::vector<ChargingPoint>& points() const noexcept {
    return points_;
  }
  [[nodiscard]] double capacity_wh() const noexcept {
    return points_.back().energy_wh;
  }
  [[nodiscard]] double time_at_energy(double energy_wh) const;
  [[nodiscard]] double charge_duration(double arrival_wh,
                                       double departure_wh) const;
  [[nodiscard]] double marginal_rate_at(double energy_wh) const;
  [[nodiscard]] double next_breakpoint_energy(double energy_wh) const;

 private:
  std::string station_type_;
  std::vector<ChargingPoint> points_;
};

// True iff charging on `left` takes no longer than charging on `right` for
// every energy interval.  This is the curve component of the TRC Section 4.3
// path-dominance rule and is independent of station names.
[[nodiscard]] bool charging_curve_no_slower(const ChargingCurve& left,
                                             const ChargingCurve& right);

struct Node {
  int id{};
  NodeKind kind{NodeKind::Customer};
  double x_km{};
  double y_km{};
  double service_time_h{};
  std::string station_type;
};

struct Vehicle {
  double speed_kmph{};
  double consumption_wh_per_km{};
  double battery_capacity_wh{};
  double max_route_duration_h{};

  [[nodiscard]] double range_km() const noexcept {
    return battery_capacity_wh / consumption_wh_per_km;
  }
};

class Instance {
 public:
  Instance(std::string name, std::string dataset, std::vector<Node> nodes,
           Vehicle vehicle,
           std::unordered_map<std::string, ChargingCurve> charging_curves);

  [[nodiscard]] const std::string& name() const noexcept { return name_; }
  [[nodiscard]] const std::string& dataset() const noexcept { return dataset_; }
  [[nodiscard]] const std::vector<Node>& nodes() const noexcept { return nodes_; }
  [[nodiscard]] const Vehicle& vehicle() const noexcept { return vehicle_; }
  [[nodiscard]] const Node& depot() const noexcept { return node(depot_id_); }
  [[nodiscard]] const std::vector<int>& customer_ids() const noexcept {
    return customer_ids_;
  }
  [[nodiscard]] const std::vector<int>& station_ids() const noexcept {
    return station_ids_;
  }
  [[nodiscard]] const Node& node(int id) const;
  [[nodiscard]] const ChargingCurve& curve_for_station(int station_id) const;
  [[nodiscard]] const std::unordered_map<std::string, ChargingCurve>&
  charging_curves() const noexcept {
    return charging_curves_;
  }
  [[nodiscard]] double distance_km(int origin_id, int target_id) const;
  [[nodiscard]] double travel_time_h(int origin_id, int target_id) const;
  [[nodiscard]] double energy_wh(int origin_id, int target_id) const;

 private:
  std::string name_;
  std::string dataset_;
  std::vector<Node> nodes_;
  Vehicle vehicle_;
  std::unordered_map<std::string, ChargingCurve> charging_curves_;
  std::unordered_map<int, std::size_t> index_by_id_;
  std::vector<std::vector<double>> distances_km_;
  int depot_id_{};
  std::vector<int> customer_ids_;
  std::vector<int> station_ids_;
};

struct Route {
  std::vector<int> visits;
  bool operator==(const Route&) const = default;
};

struct ChargingEvent {
  std::size_t visit_index{};
  int station_id{};
  double arrival_energy_wh{};
  double departure_energy_wh{};
  double duration_h{};
};

struct RouteEvaluation {
  Route route;
  bool feasible{};
  double raw_cost_h{};
  double travel_time_h{};
  double charging_time_h{};
  double service_time_h{};
  double route_duration_h{};
  double energy_shortfall_wh{};
  double duration_excess_h{};
  std::vector<double> arrival_energy_wh;
  std::vector<double> departure_energy_wh;
  std::vector<ChargingEvent> charging_events;
};

struct Plan {
  std::vector<int> customer_ids;
  // gaps[g] is the station sequence between anchor g and anchor g + 1.
  std::vector<std::vector<int>> gaps;
  Route route;
  std::shared_ptr<const RouteEvaluation> evaluation;
  bool exact_charging{false};
};

struct Solution {
  std::vector<Plan> plans;

  [[nodiscard]] bool feasible() const;
  [[nodiscard]] double raw_cost_h() const;
  [[nodiscard]] std::vector<std::vector<int>> customer_sequences() const;
  void validate_partition(const Instance& instance) const;
};

[[nodiscard]] Route expand_route(const std::vector<int>& customer_ids,
                                 const std::vector<std::vector<int>>& gaps);
[[nodiscard]] std::vector<std::vector<int>> direct_gaps(std::size_t customers);
[[nodiscard]] double generalized_cost(const RouteEvaluation& evaluation,
                                      const Instance& instance,
                                      double penalty_lambda);
[[nodiscard]] double generalized_cost(const Solution& solution,
                                      const Instance& instance,
                                      double penalty_lambda);
[[nodiscard]] std::vector<int> sorted_customer_set(
    const std::vector<int>& customer_ids);

}  // namespace ils_sp
