#include "ils_sp/core.hpp"

#include <limits>
#include <numeric>
#include <stdexcept>

namespace ils_sp {

ChargingCurve::ChargingCurve(std::string station_type,
                             std::vector<ChargingPoint> points)
    : station_type_(std::move(station_type)), points_(std::move(points)) {
  if (station_type_.empty() || points_.size() < 2) {
    throw std::invalid_argument("charging curve is incomplete");
  }
  double previous_time = -1.0;
  double previous_energy = -1.0;
  for (const ChargingPoint point : points_) {
    if (!std::isfinite(point.time_h) || !std::isfinite(point.energy_wh) ||
        point.time_h < previous_time || point.energy_wh < previous_energy) {
      throw std::invalid_argument("charging curve must be finite and monotone");
    }
    previous_time = point.time_h;
    previous_energy = point.energy_wh;
  }
  if (std::abs(points_.front().time_h) > 1e-12 ||
      std::abs(points_.front().energy_wh) > 1e-9) {
    throw std::invalid_argument("charging curve must begin at (0, 0)");
  }
}

double ChargingCurve::time_at_energy(double energy_wh) const {
  if (energy_wh < -1e-7 || energy_wh > capacity_wh() + 1e-7) {
    throw std::out_of_range("energy is outside the charging curve");
  }
  const double energy = std::clamp(energy_wh, 0.0, capacity_wh());
  for (std::size_t index = 1; index < points_.size(); ++index) {
    const ChargingPoint& left = points_[index - 1];
    const ChargingPoint& right = points_[index];
    if (energy <= right.energy_wh + 1e-9) {
      const double delta_energy = right.energy_wh - left.energy_wh;
      if (delta_energy <= 1e-12) {
        return left.time_h;
      }
      const double ratio = (energy - left.energy_wh) / delta_energy;
      return left.time_h + ratio * (right.time_h - left.time_h);
    }
  }
  return points_.back().time_h;
}

double ChargingCurve::charge_duration(double arrival_wh,
                                      double departure_wh) const {
  if (departure_wh < arrival_wh - 1e-7) {
    throw std::invalid_argument("a charging event cannot discharge the battery");
  }
  return time_at_energy(departure_wh) - time_at_energy(arrival_wh);
}

double ChargingCurve::marginal_rate_at(double energy_wh) const {
  if (energy_wh >= capacity_wh() - 1e-8) {
    return 0.0;
  }
  const double energy = std::max(0.0, energy_wh);
  for (std::size_t index = 1; index < points_.size(); ++index) {
    const ChargingPoint& left = points_[index - 1];
    const ChargingPoint& right = points_[index];
    if (energy < right.energy_wh - 1e-8) {
      const double delta_time = right.time_h - left.time_h;
      if (delta_time <= 1e-15) {
        return std::numeric_limits<double>::infinity();
      }
      return (right.energy_wh - left.energy_wh) / delta_time;
    }
  }
  return 0.0;
}

double ChargingCurve::next_breakpoint_energy(double energy_wh) const {
  for (std::size_t index = 1; index < points_.size(); ++index) {
    if (points_[index].energy_wh > energy_wh + 1e-8) {
      return points_[index].energy_wh;
    }
  }
  return capacity_wh();
}

bool charging_curve_no_slower(const ChargingCurve& left,
                              const ChargingCurve& right) {
  std::vector<double> breakpoints;
  breakpoints.reserve(left.points().size() + right.points().size());
  for (const ChargingPoint& point : left.points()) {
    breakpoints.push_back(point.energy_wh);
  }
  for (const ChargingPoint& point : right.points()) {
    breakpoints.push_back(point.energy_wh);
  }
  std::sort(breakpoints.begin(), breakpoints.end());
  breakpoints.erase(
      std::unique(breakpoints.begin(), breakpoints.end(),
                  [](double first, double second) {
                    return std::abs(first - second) <= kEnergyToleranceWh;
                  }),
      breakpoints.end());
  for (std::size_t index = 1; index < breakpoints.size(); ++index) {
    const double start = breakpoints[index - 1];
    const double finish = breakpoints[index];
    if (finish <= start + kEnergyToleranceWh) continue;
    if (left.charge_duration(start, finish) >
        right.charge_duration(start, finish) + kTimeToleranceH) {
      return false;
    }
  }
  return true;
}

Instance::Instance(
    std::string name, std::string dataset, std::vector<Node> nodes,
    Vehicle vehicle,
    std::unordered_map<std::string, ChargingCurve> charging_curves)
    : name_(std::move(name)),
      dataset_(std::move(dataset)),
      nodes_(std::move(nodes)),
      vehicle_(vehicle),
      charging_curves_(std::move(charging_curves)) {
  if (name_.empty() || nodes_.empty()) {
    throw std::invalid_argument("instance name and node set are required");
  }
  const double vehicle_values[] = {
      vehicle_.speed_kmph, vehicle_.consumption_wh_per_km,
      vehicle_.battery_capacity_wh, vehicle_.max_route_duration_h};
  for (const double value : vehicle_values) {
    if (!std::isfinite(value) || value <= 0.0) {
      throw std::invalid_argument("vehicle parameters must be positive and finite");
    }
  }

  int depot_count = 0;
  for (std::size_t index = 0; index < nodes_.size(); ++index) {
    const Node& current = nodes_[index];
    if (!index_by_id_.emplace(current.id, index).second) {
      throw std::invalid_argument("node ids must be unique");
    }
    if (!std::isfinite(current.x_km) || !std::isfinite(current.y_km) ||
        current.service_time_h < 0.0) {
      throw std::invalid_argument("invalid node geometry or service time");
    }
    switch (current.kind) {
      case NodeKind::Depot:
        ++depot_count;
        depot_id_ = current.id;
        break;
      case NodeKind::Customer:
        customer_ids_.push_back(current.id);
        break;
      case NodeKind::Station: {
        station_ids_.push_back(current.id);
        const auto curve = charging_curves_.find(current.station_type);
        if (curve == charging_curves_.end()) {
          throw std::invalid_argument("station references an unknown curve");
        }
        if (std::abs(curve->second.capacity_wh() -
                     vehicle_.battery_capacity_wh) > 1e-6) {
          throw std::invalid_argument("curve and vehicle capacities differ");
        }
        break;
      }
    }
  }
  if (depot_count != 1 || customer_ids_.empty()) {
    throw std::invalid_argument("instance needs one depot and at least one customer");
  }
  std::sort(customer_ids_.begin(), customer_ids_.end());
  std::sort(station_ids_.begin(), station_ids_.end());

  distances_km_.assign(nodes_.size(), std::vector<double>(nodes_.size(), 0.0));
  for (std::size_t left = 0; left < nodes_.size(); ++left) {
    for (std::size_t right = 0; right < nodes_.size(); ++right) {
      distances_km_[left][right] =
          std::hypot(nodes_[left].x_km - nodes_[right].x_km,
                     nodes_[left].y_km - nodes_[right].y_km);
    }
  }
}

const Node& Instance::node(int id) const {
  const auto found = index_by_id_.find(id);
  if (found == index_by_id_.end()) {
    throw std::out_of_range("unknown node id " + std::to_string(id));
  }
  return nodes_[found->second];
}

const ChargingCurve& Instance::curve_for_station(int station_id) const {
  const Node& station = node(station_id);
  if (station.kind != NodeKind::Station) {
    throw std::invalid_argument("node is not a charging station");
  }
  return charging_curves_.at(station.station_type);
}

double Instance::distance_km(int origin_id, int target_id) const {
  return distances_km_.at(index_by_id_.at(origin_id))
      .at(index_by_id_.at(target_id));
}

double Instance::travel_time_h(int origin_id, int target_id) const {
  return distance_km(origin_id, target_id) / vehicle_.speed_kmph;
}

double Instance::energy_wh(int origin_id, int target_id) const {
  return distance_km(origin_id, target_id) * vehicle_.consumption_wh_per_km;
}

bool Solution::feasible() const {
  return std::all_of(plans.begin(), plans.end(), [](const Plan& plan) {
    return plan.evaluation != nullptr && plan.evaluation->feasible;
  });
}

double Solution::raw_cost_h() const {
  return std::accumulate(plans.begin(), plans.end(), 0.0,
                         [](double total, const Plan& plan) {
                           return total + plan.evaluation->raw_cost_h;
                         });
}

std::vector<std::vector<int>> Solution::customer_sequences() const {
  std::vector<std::vector<int>> sequences;
  sequences.reserve(plans.size());
  for (const Plan& plan : plans) {
    sequences.push_back(plan.customer_ids);
  }
  return sequences;
}

void Solution::validate_partition(const Instance& instance) const {
  std::unordered_set<int> seen;
  for (const Plan& plan : plans) {
    if (plan.customer_ids.empty() || plan.gaps.size() != plan.customer_ids.size() + 1 ||
        plan.route != expand_route(plan.customer_ids, plan.gaps)) {
      throw std::invalid_argument("plan representation is inconsistent");
    }
    for (const int customer_id : plan.customer_ids) {
      if (instance.node(customer_id).kind != NodeKind::Customer ||
          !seen.insert(customer_id).second) {
        throw std::invalid_argument("solution repeats or misclassifies a customer");
      }
    }
  }
  if (seen.size() != instance.customer_ids().size()) {
    throw std::invalid_argument("solution does not cover every customer");
  }
  for (const int customer_id : instance.customer_ids()) {
    if (!seen.contains(customer_id)) {
      throw std::invalid_argument("solution misses customer " +
                                  std::to_string(customer_id));
    }
  }
}

Route expand_route(const std::vector<int>& customer_ids,
                   const std::vector<std::vector<int>>& gaps) {
  if (gaps.size() != customer_ids.size() + 1) {
    throw std::invalid_argument("a route needs one more gap than customers");
  }
  Route result;
  std::size_t station_count = 0;
  for (const auto& gap : gaps) {
    station_count += gap.size();
  }
  result.visits.reserve(customer_ids.size() + station_count);
  for (std::size_t index = 0; index < customer_ids.size(); ++index) {
    result.visits.insert(result.visits.end(), gaps[index].begin(), gaps[index].end());
    result.visits.push_back(customer_ids[index]);
  }
  result.visits.insert(result.visits.end(), gaps.back().begin(), gaps.back().end());
  return result;
}

std::vector<std::vector<int>> direct_gaps(std::size_t customers) {
  return std::vector<std::vector<int>>(customers + 1);
}

double generalized_cost(const RouteEvaluation& evaluation,
                        const Instance& instance, double penalty_lambda) {
  if (!std::isfinite(penalty_lambda) || penalty_lambda < 0.0) {
    throw std::invalid_argument("penalty lambda must be finite and nonnegative");
  }
  const double distance_violation_km =
      evaluation.energy_shortfall_wh /
      instance.vehicle().consumption_wh_per_km;
  return evaluation.raw_cost_h +
         penalty_lambda *
             (distance_violation_km + evaluation.duration_excess_h);
}

double generalized_cost(const Solution& solution, const Instance& instance,
                        double penalty_lambda) {
  double result = 0.0;
  for (const Plan& plan : solution.plans) {
    result += generalized_cost(*plan.evaluation, instance, penalty_lambda);
  }
  return result;
}

std::vector<int> sorted_customer_set(const std::vector<int>& customer_ids) {
  std::vector<int> result = customer_ids;
  std::sort(result.begin(), result.end());
  result.erase(std::unique(result.begin(), result.end()), result.end());
  return result;
}

}  // namespace ils_sp
