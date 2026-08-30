#include "ils_sp/evaluator.hpp"

#include <array>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace ils_sp {

ChargingSequenceLabel NonlinearSequenceEvaluator::edge_label(
    int origin_id, std::span<const int> station_ids, int target_id) const {
  (void)instance_.node(origin_id);
  const Node& target = instance_.node(target_id);
  ChargingSequenceLabel result{.initialized = true,
                               .origin_id = origin_id,
                               .target_id = target_id,
                               .station_legs = {},
                               .tail_energy_wh = 0.0,
                               .travel_time_h = 0.0,
                               .service_time_h = 0.0};
  result.station_legs.reserve(station_ids.size() +
                              (target.kind == NodeKind::Station ? 1U : 0U));
  int previous = origin_id;
  double energy_since_station = 0.0;
  const auto traverse = [&](int node_id, bool charging_opportunity,
                            ChargingSequenceLabel& label,
                            double& accumulated_energy, int& prior) {
    label.travel_time_h += instance_.travel_time_h(prior, node_id);
    accumulated_energy += instance_.energy_wh(prior, node_id);
    if (charging_opportunity) {
      label.station_legs.push_back(ChargingSequenceLeg{
          .station_id = node_id,
          .energy_from_previous_station_wh = accumulated_energy});
      accumulated_energy = 0.0;
    }
    prior = node_id;
  };
  for (const int station_id : station_ids) {
    if (instance_.node(station_id).kind != NodeKind::Station) {
      throw std::invalid_argument("sequence label path contains a non-station");
    }
    traverse(station_id, true, result, energy_since_station, previous);
  }
  traverse(target_id, target.kind == NodeKind::Station, result,
           energy_since_station, previous);
  result.tail_energy_wh = energy_since_station;
  if (target.kind == NodeKind::Customer) {
    result.service_time_h = target.service_time_h;
  }
  return result;
}

void NonlinearSequenceEvaluator::append(ChargingSequenceLabel& left,
                                        const ChargingSequenceLabel& right)
    const {
  if (!right.initialized) return;
  if (!left.initialized) {
    left = right;
    return;
  }
  if (left.target_id != right.origin_id) {
    throw std::invalid_argument("non-contiguous charging sequence labels");
  }
  left.target_id = right.target_id;
  left.travel_time_h += right.travel_time_h;
  left.service_time_h += right.service_time_h;
  if (right.station_legs.empty()) {
    left.tail_energy_wh += right.tail_energy_wh;
    return;
  }
  ChargingSequenceLeg first = right.station_legs.front();
  first.energy_from_previous_station_wh += left.tail_energy_wh;
  left.station_legs.push_back(first);
  left.station_legs.insert(left.station_legs.end(),
                           right.station_legs.begin() + 1,
                           right.station_legs.end());
  left.tail_energy_wh = right.tail_energy_wh;
}

ChargingSequenceLabel NonlinearSequenceEvaluator::concatenate(
    ChargingSequenceLabel left, const ChargingSequenceLabel& right) const {
  append(left, right);
  return left;
}

ChargingSequenceLabel NonlinearSequenceEvaluator::route_label(
    const Route& route) const {
  ChargingSequenceLabel result;
  int previous = instance_.depot().id;
  for (const int node_id : route.visits) {
    append(result, edge_label(previous, {}, node_id));
    previous = node_id;
  }
  append(result, edge_label(previous, {}, instance_.depot().id));
  return result;
}

RouteCostSummary NonlinearSequenceEvaluator::evaluate(
    const ChargingSequenceLabel& route) const {
  if (!route.initialized) {
    throw std::invalid_argument("cannot evaluate an empty sequence label");
  }
  return evaluate_data(route.travel_time_h, route.service_time_h,
                       route.station_legs, route.tail_energy_wh);
}

RouteCostSummary NonlinearSequenceEvaluator::evaluate(
    std::span<const ChargingSequenceLabel* const> components) const {
  if (components.empty()) {
    throw std::invalid_argument("cannot evaluate an empty sequence label");
  }
  std::size_t station_count = 0;
  double travel_time_h = 0.0;
  double service_time_h = 0.0;
  const ChargingSequenceLabel* previous = nullptr;
  for (const ChargingSequenceLabel* component : components) {
    if (component == nullptr || !component->initialized) {
      throw std::invalid_argument("sequence contains an empty component");
    }
    if (previous != nullptr &&
        previous->target_id != component->origin_id) {
      throw std::invalid_argument("non-contiguous charging sequence labels");
    }
    station_count += component->station_legs.size();
    travel_time_h += component->travel_time_h;
    service_time_h += component->service_time_h;
    previous = component;
  }

  constexpr std::size_t kInlineStations = 64;
  std::array<ChargingSequenceLeg, kInlineStations> inline_legs;
  std::vector<ChargingSequenceLeg> large_legs;
  if (station_count > kInlineStations) large_legs.resize(station_count);
  ChargingSequenceLeg* const station_legs =
      station_count <= kInlineStations ? inline_legs.data()
                                       : large_legs.data();
  std::size_t station_index = 0;
  double tail_energy_wh = 0.0;
  for (const ChargingSequenceLabel* component : components) {
    if (component->station_legs.empty()) {
      tail_energy_wh += component->tail_energy_wh;
      continue;
    }
    station_legs[station_index] = component->station_legs.front();
    station_legs[station_index].energy_from_previous_station_wh +=
        tail_energy_wh;
    ++station_index;
    for (std::size_t index = 1; index < component->station_legs.size();
         ++index) {
      station_legs[station_index++] = component->station_legs[index];
    }
    tail_energy_wh = component->tail_energy_wh;
  }
  return evaluate_data(
      travel_time_h, service_time_h,
      std::span<const ChargingSequenceLeg>{station_legs, station_count},
      tail_energy_wh);
}

RouteCostSummary NonlinearSequenceEvaluator::evaluate_data(
    double travel_time_h, double service_time_h,
    std::span<const ChargingSequenceLeg> station_legs,
    double tail_energy_wh) const {
  const double capacity = instance_.vehicle().battery_capacity_wh;
  const std::size_t station_count = station_legs.size();
  const double no_charge_duration_h = travel_time_h + service_time_h;
  if (station_count == 0) {
    const double deficit_wh = tail_energy_wh - capacity;
    const double shortfall_wh =
        deficit_wh > kEnergyToleranceWh ? deficit_wh : 0.0;
    const double route_duration_h = no_charge_duration_h;
    const double duration_excess_h = std::max(
        0.0, route_duration_h - instance_.vehicle().max_route_duration_h);
    return RouteCostSummary{
        .feasible = shortfall_wh <= kEnergyToleranceWh &&
                    duration_excess_h <= kTimeToleranceH,
        .raw_cost_h = travel_time_h,
        .travel_time_h = travel_time_h,
        .charging_time_h = 0.0,
        .service_time_h = service_time_h,
        .route_duration_h = route_duration_h,
        .energy_shortfall_wh = shortfall_wh,
        .duration_excess_h = duration_excess_h};
  }

  constexpr std::size_t kInlineStations = 64;
  const std::size_t count = station_count + 2;
  const bool inline_workspace = count <= kInlineStations;
  std::array<double, kInlineStations> inline_arrival;
  std::array<double, kInlineStations> inline_departure;
  std::array<double, kInlineStations> inline_capacity_slack;
  std::array<std::size_t, kInlineStations> inline_preceding;
  std::vector<double> large_arrival;
  std::vector<double> large_departure;
  std::vector<double> large_capacity_slack;
  std::vector<std::size_t> large_preceding;
  if (!inline_workspace) {
    large_arrival.resize(count);
    large_departure.resize(count);
    large_capacity_slack.resize(count);
    large_preceding.reserve(station_count);
  }
  double* const arrival = inline_workspace ? inline_arrival.data()
                                            : large_arrival.data();
  double* const departure = inline_workspace ? inline_departure.data()
                                              : large_departure.data();
  double* const capacity_slack =
      inline_workspace ? inline_capacity_slack.data()
                       : large_capacity_slack.data();
  std::fill_n(arrival, count, 0.0);
  std::fill_n(departure, count, 0.0);
  departure[0] = capacity;
  double remaining_shift_h = std::max(
      0.0,
      instance_.vehicle().max_route_duration_h - no_charge_duration_h);
  double shortfall_wh = 0.0;
  std::size_t preceding_count = 0;
  const auto preceding_station = [&](std::size_t index) {
    return inline_workspace ? inline_preceding[index]
                            : large_preceding[index];
  };
  const auto add_preceding_station = [&](std::size_t position) {
    if (inline_workspace) {
      inline_preceding[preceding_count] = position;
    } else {
      large_preceding.push_back(position);
    }
    ++preceding_count;
  };

  for (std::size_t position = 1; position < count; ++position) {
    const double energy_required =
        position <= station_count
            ? station_legs[position - 1].energy_from_previous_station_wh
            : tail_energy_wh;
    arrival[position] = departure[position - 1] - energy_required;
    departure[position] = arrival[position];
    while (arrival[position] < -kEnergyToleranceWh) {
      const double required_wh = -arrival[position];
      if (preceding_count == 0) {
        shortfall_wh += required_wh;
        arrival[position] += required_wh;
        departure[position] += required_wh;
        break;
      }

      double suffix_slack = std::numeric_limits<double>::infinity();
      const std::size_t earliest_station = preceding_station(0);
      for (std::size_t reverse = position; reverse-- > earliest_station;) {
        suffix_slack =
            std::min(suffix_slack, capacity - departure[reverse]);
        capacity_slack[reverse] = suffix_slack;
      }

      std::optional<std::size_t> best_station;
      double best_rate = -1.0;
      double best_available_wh = 0.0;
      for (std::size_t station_index = 0; station_index < preceding_count;
           ++station_index) {
        const std::size_t station_position =
            preceding_station(station_index);
        const int station_id =
            station_legs[station_position - 1].station_id;
        const ChargingCurve& curve = instance_.curve_for_station(station_id);
        const double current_energy = departure[station_position];
        if (current_energy < -kEnergyToleranceWh) continue;
        const double segment_slack =
            curve.next_breakpoint_energy(current_energy) - current_energy;
        double available_wh =
            std::min(capacity_slack[station_position], segment_slack);
        if (available_wh <= kEnergyToleranceWh) continue;
        const double rate = curve.marginal_rate_at(current_energy);
        if (!std::isinf(rate)) {
          available_wh = std::min(available_wh, remaining_shift_h * rate);
        }
        if (available_wh <= kEnergyToleranceWh) continue;
        if (rate > best_rate + 1e-9 ||
            (std::abs(rate - best_rate) <= 1e-9 &&
             (!best_station.has_value() ||
              station_position > *best_station))) {
          best_station = station_position;
          best_rate = rate;
          best_available_wh = available_wh;
        }
      }

      if (!best_station.has_value()) {
        shortfall_wh += required_wh;
        arrival[position] += required_wh;
        departure[position] += required_wh;
        break;
      }

      const double charged_wh = std::min(required_wh, best_available_wh);
      const int station_id =
          station_legs[*best_station - 1].station_id;
      const ChargingCurve& curve = instance_.curve_for_station(station_id);
      const double charge_duration_h = curve.charge_duration(
          departure[*best_station], departure[*best_station] + charged_wh);
      departure[*best_station] += charged_wh;
      for (std::size_t suffix = *best_station + 1; suffix <= position;
           ++suffix) {
        arrival[suffix] += charged_wh;
        departure[suffix] += charged_wh;
      }
      remaining_shift_h =
          std::max(0.0, remaining_shift_h - charge_duration_h);
    }

    if (position <= station_count) {
      add_preceding_station(position);
    }
  }

  double charging_time_h = 0.0;
  for (std::size_t position = 1; position <= station_count; ++position) {
    const int station_id = station_legs[position - 1].station_id;
    const double start = std::clamp(arrival[position], 0.0, capacity);
    const double finish = std::clamp(departure[position], start, capacity);
    charging_time_h += instance_.curve_for_station(station_id).charge_duration(
        start, finish);
  }
  const double raw_cost_h = travel_time_h + charging_time_h;
  const double route_duration_h = raw_cost_h + service_time_h;
  const double duration_excess_h = std::max(
      {no_charge_duration_h - instance_.vehicle().max_route_duration_h,
       route_duration_h - instance_.vehicle().max_route_duration_h, 0.0});
  return RouteCostSummary{
      .feasible = shortfall_wh <= kEnergyToleranceWh &&
                  duration_excess_h <= kTimeToleranceH,
      .raw_cost_h = raw_cost_h,
      .travel_time_h = travel_time_h,
      .charging_time_h = charging_time_h,
      .service_time_h = service_time_h,
      .route_duration_h = route_duration_h,
      .energy_shortfall_wh = shortfall_wh,
      .duration_excess_h = duration_excess_h};
}

double generalized_cost(const RouteCostSummary& summary,
                        const Instance& instance, double penalty_lambda) {
  if (!std::isfinite(penalty_lambda) || penalty_lambda < 0.0) {
    throw std::invalid_argument("penalty lambda must be finite and nonnegative");
  }
  return summary.raw_cost_h +
         penalty_lambda *
             (summary.energy_shortfall_wh /
                  instance.vehicle().consumption_wh_per_km +
              summary.duration_excess_h);
}

double NonlinearSequenceEvaluator::generalized_cost(
    const ChargingSequenceLabel& route, double penalty_lambda) const {
  return ils_sp::generalized_cost(evaluate(route), instance_, penalty_lambda);
}

double NonlinearSequenceEvaluator::generalized_cost(
    std::span<const ChargingSequenceLabel* const> components,
    double penalty_lambda) const {
  return ils_sp::generalized_cost(evaluate(components), instance_,
                                  penalty_lambda);
}

RouteEvaluator::RouteEvaluator(const Instance& instance,
                               std::size_t cache_capacity)
    : instance_(instance), cache_(cache_capacity) {}

void RouteEvaluator::validate_route(const Route& route) const {
  std::unordered_set<int> customers;
  for (const int node_id : route.visits) {
    const Node& current = instance_.node(node_id);
    if (current.kind == NodeKind::Depot) {
      throw std::invalid_argument("the depot is implicit in a Route");
    }
    if (current.kind == NodeKind::Customer &&
        !customers.insert(node_id).second) {
      throw std::invalid_argument("a route repeats a customer");
    }
  }
}

std::shared_ptr<const RouteEvaluation> RouteEvaluator::evaluate(
    const Route& route) {
  if (auto cached = cache_.get(route.visits); cached.has_value()) {
    return *cached;
  }
  validate_route(route);
  auto result = std::make_shared<const RouteEvaluation>(evaluate_uncached(route));
  cache_.put(route.visits, result);
  return result;
}

std::shared_ptr<const RouteEvaluation> RouteEvaluator::evaluate_schedule(
    const Route& route,
    const std::vector<std::optional<double>>& station_departure_wh_by_visit) {
  validate_route(route);
  if (station_departure_wh_by_visit.size() != route.visits.size()) {
    throw std::invalid_argument("charging schedule length differs from route");
  }
  for (std::size_t index = 0; index < route.visits.size(); ++index) {
    const bool station =
        instance_.node(route.visits[index]).kind == NodeKind::Station;
    if (station != station_departure_wh_by_visit[index].has_value()) {
      throw std::invalid_argument(
          "charging schedule positions differ from station visits");
    }
  }

  const int depot_id = instance_.depot().id;
  std::vector<int> path;
  path.reserve(route.visits.size() + 2);
  path.push_back(depot_id);
  path.insert(path.end(), route.visits.begin(), route.visits.end());
  path.push_back(depot_id);
  const std::size_t count = path.size();
  const double capacity = instance_.vehicle().battery_capacity_wh;
  std::vector<double> arrival(count, 0.0);
  std::vector<double> departure(count, 0.0);
  departure.front() = capacity;
  double travel_time_h = 0.0;
  double charging_time_h = 0.0;
  std::vector<ChargingEvent> events;

  for (std::size_t position = 1; position < count; ++position) {
    const int origin = path[position - 1];
    const int target = path[position];
    travel_time_h += instance_.travel_time_h(origin, target);
    const double energy =
        departure[position - 1] - instance_.energy_wh(origin, target);
    if (energy < -kEnergyToleranceWh) {
      throw std::invalid_argument("charging schedule contains an unreachable leg");
    }
    arrival[position] = std::max(0.0, energy);
    departure[position] = arrival[position];
    if (position == count - 1 ||
        instance_.node(target).kind != NodeKind::Station) {
      continue;
    }
    const std::size_t visit_index = position - 1;
    const double scheduled = *station_departure_wh_by_visit[visit_index];
    if (!std::isfinite(scheduled) ||
        scheduled < arrival[position] - kEnergyToleranceWh ||
        scheduled > capacity + kEnergyToleranceWh) {
      throw std::invalid_argument("invalid station departure energy");
    }
    departure[position] = std::clamp(scheduled, arrival[position], capacity);
    const ChargingCurve& curve = instance_.curve_for_station(target);
    const double duration =
        curve.charge_duration(arrival[position], departure[position]);
    charging_time_h += duration;
    if (departure[position] > arrival[position] + kEnergyToleranceWh) {
      events.push_back(ChargingEvent{.visit_index = visit_index,
                                     .station_id = target,
                                     .arrival_energy_wh = arrival[position],
                                     .departure_energy_wh = departure[position],
                                     .duration_h = duration});
    }
  }

  double service_time_h = 0.0;
  for (const int node_id : route.visits) {
    const Node& current = instance_.node(node_id);
    if (current.kind == NodeKind::Customer) {
      service_time_h += current.service_time_h;
    }
  }
  const double raw_cost_h = travel_time_h + charging_time_h;
  const double route_duration_h = raw_cost_h + service_time_h;
  const double excess = std::max(
      0.0, route_duration_h - instance_.vehicle().max_route_duration_h);
  auto result = std::make_shared<const RouteEvaluation>(RouteEvaluation{
      .route = route,
      .feasible = excess <= kTimeToleranceH,
      .raw_cost_h = raw_cost_h,
      .travel_time_h = travel_time_h,
      .charging_time_h = charging_time_h,
      .service_time_h = service_time_h,
      .route_duration_h = route_duration_h,
      .energy_shortfall_wh = 0.0,
      .duration_excess_h = excess,
      .arrival_energy_wh = std::move(arrival),
      .departure_energy_wh = std::move(departure),
      .charging_events = std::move(events)});
  return result;
}

RouteEvaluation RouteEvaluator::evaluate_uncached(const Route& route) const {
  const int depot_id = instance_.depot().id;
  std::vector<int> path;
  path.reserve(route.visits.size() + 2);
  path.push_back(depot_id);
  path.insert(path.end(), route.visits.begin(), route.visits.end());
  path.push_back(depot_id);
  const std::size_t count = path.size();
  const double capacity = instance_.vehicle().battery_capacity_wh;

  std::vector<double> arrival(count, 0.0);
  std::vector<double> departure(count, 0.0);
  departure.front() = capacity;
  double travel_time_h = 0.0;
  for (std::size_t position = 1; position < count; ++position) {
    travel_time_h += instance_.travel_time_h(path[position - 1], path[position]);
    arrival[position] = departure[position - 1] -
                        instance_.energy_wh(path[position - 1], path[position]);
    departure[position] = arrival[position];
  }

  double service_time_h = 0.0;
  for (const int node_id : route.visits) {
    const Node& current = instance_.node(node_id);
    if (current.kind == NodeKind::Customer) {
      service_time_h += current.service_time_h;
    }
  }
  const double no_charge_duration_h = travel_time_h + service_time_h;
  double duration_excess_h = std::max(
      0.0,
      no_charge_duration_h - instance_.vehicle().max_route_duration_h);
  double remaining_shift_h = std::max(
      0.0,
      instance_.vehicle().max_route_duration_h - no_charge_duration_h);

  double shortfall_wh = 0.0;
  std::vector<std::size_t> preceding_stations;
  for (std::size_t position = 1; position < count; ++position) {
    while (arrival[position] < -kEnergyToleranceWh) {
      const double required_wh = -arrival[position];
      if (preceding_stations.empty()) {
        shortfall_wh += required_wh;
        for (std::size_t suffix = position; suffix < count; ++suffix) {
          arrival[suffix] += required_wh;
          departure[suffix] += required_wh;
        }
        break;
      }

      std::vector<double> capacity_slack(position, 0.0);
      double suffix_slack = std::numeric_limits<double>::infinity();
      const std::size_t earliest_station = preceding_stations.front();
      for (std::size_t reverse = position; reverse-- > earliest_station;) {
        suffix_slack =
            std::min(suffix_slack, capacity - departure[reverse]);
        capacity_slack[reverse] = suffix_slack;
      }

      std::optional<std::size_t> best_station;
      double best_rate = -1.0;
      double best_available_wh = 0.0;
      for (const std::size_t station_position : preceding_stations) {
        const int station_id = path[station_position];
        const ChargingCurve& curve = instance_.curve_for_station(station_id);
        const double current_energy = departure[station_position];
        if (current_energy < -kEnergyToleranceWh) {
          continue;
        }
        const double segment_slack =
            curve.next_breakpoint_energy(current_energy) - current_energy;
        double available_wh =
            std::min(capacity_slack[station_position], segment_slack);
        if (available_wh <= kEnergyToleranceWh) {
          continue;
        }
        const double rate = curve.marginal_rate_at(current_energy);
        if (!std::isinf(rate)) {
          available_wh = std::min(available_wh, remaining_shift_h * rate);
        }
        if (available_wh <= kEnergyToleranceWh) {
          continue;
        }
        if (rate > best_rate + 1e-9 ||
            (std::abs(rate - best_rate) <= 1e-9 &&
             (!best_station.has_value() ||
              station_position > *best_station))) {
          best_station = station_position;
          best_rate = rate;
          best_available_wh = available_wh;
        }
      }

      if (!best_station.has_value()) {
        shortfall_wh += required_wh;
        for (std::size_t suffix = position; suffix < count; ++suffix) {
          arrival[suffix] += required_wh;
          departure[suffix] += required_wh;
        }
        break;
      }

      const double charged_wh = std::min(required_wh, best_available_wh);
      const ChargingCurve& curve =
          instance_.curve_for_station(path[*best_station]);
      const double charge_duration_h = curve.charge_duration(
          departure[*best_station], departure[*best_station] + charged_wh);
      departure[*best_station] += charged_wh;
      for (std::size_t suffix = *best_station + 1; suffix < count; ++suffix) {
        arrival[suffix] += charged_wh;
        departure[suffix] += charged_wh;
      }
      remaining_shift_h =
          std::max(0.0, remaining_shift_h - charge_duration_h);
    }

    if (position + 1 < count &&
        instance_.node(path[position]).kind == NodeKind::Station) {
      preceding_stations.push_back(position);
    }
  }

  std::vector<ChargingEvent> events;
  double charging_time_h = 0.0;
  for (std::size_t position = 1; position + 1 < count; ++position) {
    const int node_id = path[position];
    if (instance_.node(node_id).kind != NodeKind::Station) {
      continue;
    }
    const double start = std::clamp(arrival[position], 0.0, capacity);
    const double end = std::clamp(departure[position], start, capacity);
    const double duration_h =
        instance_.curve_for_station(node_id).charge_duration(start, end);
    charging_time_h += duration_h;
    if (end > start + kEnergyToleranceWh) {
      events.push_back(ChargingEvent{.visit_index = position - 1,
                                     .station_id = node_id,
                                     .arrival_energy_wh = start,
                                     .departure_energy_wh = end,
                                     .duration_h = duration_h});
    }
  }
  const double raw_cost_h = travel_time_h + charging_time_h;
  const double route_duration_h = raw_cost_h + service_time_h;
  duration_excess_h =
      std::max({duration_excess_h,
                route_duration_h - instance_.vehicle().max_route_duration_h,
                0.0});
  return RouteEvaluation{
      .route = route,
      .feasible = shortfall_wh <= kEnergyToleranceWh &&
                  duration_excess_h <= kTimeToleranceH,
      .raw_cost_h = raw_cost_h,
      .travel_time_h = travel_time_h,
      .charging_time_h = charging_time_h,
      .service_time_h = service_time_h,
      .route_duration_h = route_duration_h,
      .energy_shortfall_wh = shortfall_wh,
      .duration_excess_h = duration_excess_h,
      .arrival_energy_wh = std::move(arrival),
      .departure_energy_wh = std::move(departure),
      .charging_events = std::move(events)};
}

std::size_t RouteEvaluator::cache_size() const noexcept { return cache_.size(); }
std::size_t RouteEvaluator::cache_hits() const noexcept { return cache_.hits(); }
std::size_t RouteEvaluator::cache_misses() const noexcept {
  return cache_.misses();
}
std::size_t RouteEvaluator::cache_evictions() const noexcept {
  return cache_.evictions();
}

PlanFactory::PlanFactory(const Instance& instance, RouteEvaluator& evaluator,
                         std::size_t exact_cache_capacity)
    : instance_(instance),
      evaluator_(evaluator),
      exact_cache_(exact_cache_capacity) {}

Plan PlanFactory::make_plan(std::vector<int> customer_ids,
                            std::vector<std::vector<int>> gaps) {
  if (customer_ids.empty()) {
    throw std::invalid_argument("a route plan needs at least one customer");
  }
  std::unordered_set<int> seen;
  for (const int customer_id : customer_ids) {
    if (instance_.node(customer_id).kind != NodeKind::Customer ||
        !seen.insert(customer_id).second) {
      throw std::invalid_argument("invalid customer sequence");
    }
  }
  for (const auto& gap : gaps) {
    std::unordered_set<int> gap_stations;
    for (const int station_id : gap) {
      if (instance_.node(station_id).kind != NodeKind::Station ||
          !gap_stations.insert(station_id).second) {
        throw std::invalid_argument("a gap must be a simple station path");
      }
    }
  }
  Route route = expand_route(customer_ids, gaps);
  return Plan{.customer_ids = std::move(customer_ids),
              .gaps = std::move(gaps),
              .route = route,
              .evaluation = evaluator_.evaluate(route),
              .exact_charging = false};
}

std::optional<Plan> PlanFactory::exact_plan(
    const std::vector<int>& customer_ids) {
  return exact_cache_.get(customer_ids);
}

PlanFactory::GapMap PlanFactory::gap_map(const Solution& reference) const {
  GapMap result;
  const int depot_id = instance_.depot().id;
  // A gap is inherited only in the direction in which it exists in the
  // incumbent.  Reversing an anchor pair creates a new direct gap.
  for (const Plan& plan : reference.plans) {
    std::vector<int> anchors;
    anchors.reserve(plan.customer_ids.size() + 2);
    anchors.push_back(depot_id);
    anchors.insert(anchors.end(), plan.customer_ids.begin(),
                   plan.customer_ids.end());
    anchors.push_back(depot_id);
    for (std::size_t gap = 0; gap < plan.gaps.size(); ++gap) {
      result.try_emplace(AnchorPair{anchors[gap], anchors[gap + 1]},
                         plan.gaps[gap]);
    }
  }
  return result;
}

Plan PlanFactory::plan_sequence(const std::vector<int>& customer_ids,
                                const Solution* reference) {
  if (reference != nullptr) {
    return plan_sequence(customer_ids, gap_map(*reference));
  }
  std::vector<std::vector<int>> gaps = direct_gaps(customer_ids.size());
  return make_plan(customer_ids, std::move(gaps));
}

Plan PlanFactory::plan_sequence(const std::vector<int>& customer_ids,
                                const GapMap& reference_gaps) {
  std::vector<std::vector<int>> gaps = direct_gaps(customer_ids.size());
  std::vector<int> anchors;
  anchors.reserve(customer_ids.size() + 2);
  anchors.push_back(instance_.depot().id);
  anchors.insert(anchors.end(), customer_ids.begin(), customer_ids.end());
  anchors.push_back(instance_.depot().id);
  for (std::size_t gap = 0; gap + 1 < anchors.size(); ++gap) {
    const auto found =
        reference_gaps.find(AnchorPair{anchors[gap], anchors[gap + 1]});
    if (found != reference_gaps.end()) {
      gaps[gap] = found->second;
    }
  }
  return make_plan(customer_ids, std::move(gaps));
}

void PlanFactory::publish_exact(Plan plan) {
  if (!plan.exact_charging || plan.evaluation == nullptr ||
      !plan.evaluation->feasible) {
    throw std::invalid_argument("exact cache only accepts feasible exact plans");
  }
  const std::vector<int> customer_key = plan.customer_ids;
  exact_cache_.put(customer_key, std::move(plan));
}

struct RouteSequenceScorer::Implementation {
  static constexpr std::size_t kMaximumLabelComponents = 512;

  struct CustomerLocation {
    std::size_t route{};
    std::size_t position{};
  };

  struct LabelTree {
    std::size_t base{1};
    std::size_t edge_count{};
    std::vector<ChargingSequenceLabel> nodes;

    LabelTree() = default;

    LabelTree(std::vector<ChargingSequenceLabel> leaves,
              const NonlinearSequenceEvaluator& evaluator)
        : edge_count(leaves.size()) {
      while (base < edge_count) base *= 2;
      nodes.resize(base * 2);
      for (std::size_t index = 0; index < leaves.size(); ++index) {
        nodes[base + index] = std::move(leaves[index]);
      }
      for (std::size_t index = base; index-- > 1;) {
        nodes[index] = evaluator.concatenate(nodes[index * 2],
                                             nodes[index * 2 + 1]);
      }
    }

    void append_range(std::size_t begin, std::size_t end,
                      ChargingSequenceLabel& result,
                      const NonlinearSequenceEvaluator& evaluator) const {
      if (begin > end || end > edge_count) {
        throw std::out_of_range("sequence-label range is invalid");
      }
      if (begin == end) return;
      std::array<std::size_t, 128> left_nodes{};
      std::array<std::size_t, 128> right_nodes{};
      std::size_t left_count = 0;
      std::size_t right_count = 0;
      std::size_t left = begin + base;
      std::size_t right = end + base;
      while (left < right) {
        if ((left & 1U) != 0U) left_nodes[left_count++] = left++;
        if ((right & 1U) != 0U) right_nodes[right_count++] = --right;
        left /= 2;
        right /= 2;
      }
      for (std::size_t index = 0; index < left_count; ++index) {
        evaluator.append(result, nodes[left_nodes[index]]);
      }
      while (right_count > 0) {
        --right_count;
        evaluator.append(result, nodes[right_nodes[right_count]]);
      }
    }

    void append_range_components(
        std::size_t begin, std::size_t end,
        std::span<const ChargingSequenceLabel*> components,
        std::size_t& component_count) const {
      if (begin > end || end > edge_count) {
        throw std::out_of_range("sequence-label range is invalid");
      }
      if (begin == end) return;
      std::array<std::size_t, 128> left_nodes;
      std::array<std::size_t, 128> right_nodes;
      std::size_t left_count = 0;
      std::size_t right_count = 0;
      std::size_t left = begin + base;
      std::size_t right = end + base;
      while (left < right) {
        if ((left & 1U) != 0U) left_nodes[left_count++] = left++;
        if ((right & 1U) != 0U) right_nodes[right_count++] = --right;
        left /= 2;
        right /= 2;
      }
      const auto append_node = [&](std::size_t node) {
        if (component_count >= components.size()) {
          throw std::length_error("route label exceeds component bound");
        }
        components[component_count++] = &nodes[node];
      };
      for (std::size_t index = 0; index < left_count; ++index) {
        append_node(left_nodes[index]);
      }
      while (right_count > 0) append_node(right_nodes[--right_count]);
    }
  };

  struct RouteProfile {
    std::vector<int> physical_visits;
    std::vector<int> customers;
    LabelTree forward;
    LabelTree reverse;

    RouteProfile(const Plan& plan,
                 const NonlinearSequenceEvaluator& evaluator)
        : physical_visits(plan.route.visits), customers(plan.customer_ids) {
      std::vector<ChargingSequenceLabel> forward_edges;
      std::vector<ChargingSequenceLabel> reverse_edges;
      if (customers.size() > 1) {
        forward_edges.reserve(customers.size() - 1);
        reverse_edges.reserve(customers.size() - 1);
      }
      for (std::size_t index = 0; index + 1 < customers.size(); ++index) {
        forward_edges.push_back(evaluator.edge_label(
            customers[index], plan.gaps[index + 1], customers[index + 1]));
      }
      for (std::size_t index = customers.size(); index > 1; --index) {
        reverse_edges.push_back(evaluator.edge_label(
            customers[index - 1], {}, customers[index - 2]));
      }
      forward = LabelTree(std::move(forward_edges), evaluator);
      reverse = LabelTree(std::move(reverse_edges), evaluator);
    }

    void append_slice(std::size_t begin, std::size_t end, bool reversed,
                      ChargingSequenceLabel& result,
                      const NonlinearSequenceEvaluator& evaluator) const {
      if (begin >= end || end > customers.size()) {
        throw std::out_of_range("customer subsequence is invalid");
      }
      if (end == begin + 1) return;
      if (!reversed) {
        forward.append_range(begin, end - 1, result, evaluator);
        return;
      }
      reverse.append_range(customers.size() - end,
                           customers.size() - begin - 1, result, evaluator);
    }

    void append_slice_components(
        std::size_t begin, std::size_t end, bool reversed,
        std::span<const ChargingSequenceLabel*> components,
        std::size_t& component_count) const {
      if (begin >= end || end > customers.size()) {
        throw std::out_of_range("customer subsequence is invalid");
      }
      if (end == begin + 1) return;
      if (!reversed) {
        forward.append_range_components(begin, end - 1, components,
                                        component_count);
        return;
      }
      reverse.append_range_components(customers.size() - end,
                                      customers.size() - begin - 1,
                                      components, component_count);
    }
  };

  const Instance& instance;
  NonlinearSequenceEvaluator evaluator;
  std::vector<RouteProfile> profiles;
  std::unordered_map<int, CustomerLocation> locations;
  std::vector<std::size_t> anchor_index_by_id;
  std::size_t anchor_count{};
  std::vector<ChargingSequenceLabel> direct_edges;
  std::vector<std::size_t> reference_edge_index;
  std::vector<std::size_t> active_reference_edges;
  std::vector<ChargingSequenceLabel> reference_edges;

  Implementation(const Instance& instance_value, const Solution& reference,
                 const PlanFactory::GapMap& gaps)
      : instance(instance_value),
        evaluator(instance_value) {
    std::vector<int> anchors;
    anchors.reserve(instance.customer_ids().size() + 1);
    anchors.push_back(instance.depot().id);
    anchors.insert(anchors.end(), instance.customer_ids().begin(),
                   instance.customer_ids().end());
    if (std::any_of(anchors.begin(), anchors.end(),
                    [](int id) { return id < 0; })) {
      throw std::invalid_argument(
          "Montoya anchor ids must be nonnegative for dense scoring");
    }
    const int maximum_anchor =
        *std::max_element(anchors.begin(), anchors.end());
    const std::size_t missing = std::numeric_limits<std::size_t>::max();
    anchor_index_by_id.assign(static_cast<std::size_t>(maximum_anchor) + 1,
                              missing);
    for (std::size_t index = 0; index < anchors.size(); ++index) {
      anchor_index_by_id[static_cast<std::size_t>(anchors[index])] = index;
    }
    anchor_count = anchors.size();
    direct_edges.reserve(anchor_count * anchor_count);
    for (const int origin : anchors) {
      for (const int target : anchors) {
        direct_edges.push_back(evaluator.edge_label(origin, {}, target));
      }
    }
    reference_edge_index.assign(direct_edges.size(), missing);
    reset(reference, gaps);
  }

  std::size_t dense_edge_index(int origin_id, int target_id) const {
    const std::size_t missing = std::numeric_limits<std::size_t>::max();
    const auto anchor_index = [&](int id) {
      if (id < 0 || static_cast<std::size_t>(id) >= anchor_index_by_id.size() ||
          anchor_index_by_id[static_cast<std::size_t>(id)] == missing) {
        throw std::out_of_range("route scorer edge has a non-anchor endpoint");
      }
      return anchor_index_by_id[static_cast<std::size_t>(id)];
    };
    return anchor_index(origin_id) * anchor_count + anchor_index(target_id);
  }

  void reset(const Solution& reference, const PlanFactory::GapMap& gaps) {
    std::vector<bool> reused(profiles.size(), false);
    std::vector<RouteProfile> updated_profiles;
    updated_profiles.reserve(reference.plans.size());
    for (const Plan& plan : reference.plans) {
      std::optional<std::size_t> reusable;
      for (std::size_t index = 0; index < profiles.size(); ++index) {
        if (!reused[index] &&
            profiles[index].physical_visits == plan.route.visits &&
            profiles[index].customers == plan.customer_ids) {
          reusable = index;
          break;
        }
      }
      if (reusable.has_value()) {
        reused[*reusable] = true;
        updated_profiles.push_back(std::move(profiles[*reusable]));
      } else {
        updated_profiles.emplace_back(plan, evaluator);
      }
    }
    profiles = std::move(updated_profiles);
    const std::size_t missing = std::numeric_limits<std::size_t>::max();
    for (const std::size_t edge_index : active_reference_edges) {
      reference_edge_index[edge_index] = missing;
    }
    active_reference_edges.clear();
    reference_edges.clear();
    reference_edges.reserve(gaps.size());
    active_reference_edges.reserve(gaps.size());
    for (const auto& [pair, path] : gaps) {
      if (path.empty()) continue;
      const std::size_t edge_index = dense_edge_index(pair.origin, pair.target);
      reference_edge_index[edge_index] = reference_edges.size();
      active_reference_edges.push_back(edge_index);
      reference_edges.push_back(
          evaluator.edge_label(pair.origin, path, pair.target));
    }
    locations.clear();
    locations.reserve(instance.customer_ids().size());
    for (std::size_t route = 0; route < reference.plans.size(); ++route) {
      const Plan& plan = reference.plans[route];
      for (std::size_t position = 0; position < plan.customer_ids.size();
           ++position) {
        locations.emplace(plan.customer_ids[position],
                          CustomerLocation{route, position});
      }
    }
  }

  const ChargingSequenceLabel& edge(int origin_id, int target_id) const {
    const std::size_t edge_index = dense_edge_index(origin_id, target_id);
    const std::size_t override_index = reference_edge_index[edge_index];
    return override_index == std::numeric_limits<std::size_t>::max()
               ? direct_edges[edge_index]
               : reference_edges[override_index];
  }

  ChargingSequenceLabel label(
      const std::vector<int>& customer_sequence) const {
    if (customer_sequence.empty()) {
      throw std::invalid_argument("a scored route needs at least one customer");
    }
    ChargingSequenceLabel result;
    int previous = instance.depot().id;
    std::size_t sequence_position = 0;
    while (sequence_position < customer_sequence.size()) {
      const int first_customer = customer_sequence[sequence_position];
      const auto first_location = locations.find(first_customer);
      if (first_location == locations.end()) {
        evaluator.append(result, edge(previous, first_customer));
        previous = first_customer;
        ++sequence_position;
        continue;
      }

      const CustomerLocation first = first_location->second;
      int direction = 0;
      if (sequence_position + 1 < customer_sequence.size()) {
        const auto next_location =
            locations.find(customer_sequence[sequence_position + 1]);
        if (next_location != locations.end() &&
            next_location->second.route == first.route) {
          if (next_location->second.position == first.position + 1) {
            direction = 1;
          } else if (first.position == next_location->second.position + 1) {
            direction = -1;
          }
        }
      }

      std::size_t end_position = sequence_position + 1;
      std::size_t last_original_position = first.position;
      if (direction != 0) {
        while (end_position < customer_sequence.size()) {
          const auto next_location =
              locations.find(customer_sequence[end_position]);
          if (next_location == locations.end() ||
              next_location->second.route != first.route) {
            break;
          }
          const std::ptrdiff_t expected =
              static_cast<std::ptrdiff_t>(first.position) +
              static_cast<std::ptrdiff_t>(direction) *
                  static_cast<std::ptrdiff_t>(end_position -
                                              sequence_position);
          if (expected < 0 ||
              next_location->second.position !=
                  static_cast<std::size_t>(expected)) {
            break;
          }
          last_original_position = next_location->second.position;
          ++end_position;
        }
      }

      evaluator.append(result, edge(previous, first_customer));
      if (direction >= 0) {
        profiles[first.route].append_slice(
            first.position, last_original_position + 1, false, result,
            evaluator);
      } else {
        profiles[first.route].append_slice(
            last_original_position, first.position + 1, true, result,
            evaluator);
      }
      previous = customer_sequence[end_position - 1];
      sequence_position = end_position;
    }
    evaluator.append(result, edge(previous, instance.depot().id));
    return result;
  }

  ChargingSequenceLabel insertion_label(std::size_t route_index,
                                        std::size_t position,
                                        int customer_id) const {
    if (route_index >= profiles.size()) {
      throw std::out_of_range("insertion scorer has invalid route");
    }
    if (locations.contains(customer_id)) {
      throw std::invalid_argument(
          "insertion scorer customer is already in the incumbent");
    }
    const RouteProfile& profile = profiles[route_index];
    if (position > profile.customers.size()) {
      throw std::out_of_range("insertion scorer has invalid position");
    }

    ChargingSequenceLabel result;
    int previous = instance.depot().id;
    if (position > 0) {
      evaluator.append(result,
                       edge(instance.depot().id, profile.customers.front()));
      profile.append_slice(0, position, false, result, evaluator);
      previous = profile.customers[position - 1];
    }
    evaluator.append(result, edge(previous, customer_id));
    previous = customer_id;
    if (position < profile.customers.size()) {
      evaluator.append(result,
                       edge(customer_id, profile.customers[position]));
      profile.append_slice(position, profile.customers.size(), false, result,
                           evaluator);
      previous = profile.customers.back();
    }
    evaluator.append(result, edge(previous, instance.depot().id));
    return result;
  }

  ChargingSequenceLabel label(
      std::span<const IncumbentRouteSlice> slices) const {
    if (slices.empty()) {
      throw std::invalid_argument("a scored route needs at least one slice");
    }
    ChargingSequenceLabel result;
    int previous = instance.depot().id;
    bool has_customer = false;
    for (const IncumbentRouteSlice& slice : slices) {
      if (slice.route_index >= profiles.size()) {
        throw std::out_of_range("incumbent route slice has invalid route");
      }
      const RouteProfile& profile = profiles[slice.route_index];
      if (slice.begin >= slice.end || slice.end > profile.customers.size()) {
        throw std::out_of_range("incumbent route slice has invalid bounds");
      }
      const int first_customer =
          slice.reversed ? profile.customers[slice.end - 1]
                         : profile.customers[slice.begin];
      evaluator.append(result, edge(previous, first_customer));
      profile.append_slice(slice.begin, slice.end, slice.reversed, result,
                           evaluator);
      previous = slice.reversed ? profile.customers[slice.begin]
                                : profile.customers[slice.end - 1];
      has_customer = true;
    }
    if (!has_customer) {
      throw std::invalid_argument("a scored route needs at least one customer");
    }
    evaluator.append(result, edge(previous, instance.depot().id));
    return result;
  }

  std::size_t collect_components(
      std::span<const IncumbentRouteSlice> slices,
      std::span<const ChargingSequenceLabel*> components) const {
    if (slices.empty()) {
      throw std::invalid_argument("a scored route needs at least one slice");
    }
    std::size_t component_count = 0;
    const auto append_component = [&](const ChargingSequenceLabel& label) {
      if (component_count >= components.size()) {
        throw std::length_error("route label exceeds component bound");
      }
      components[component_count++] = &label;
    };
    int previous = instance.depot().id;
    for (const IncumbentRouteSlice& slice : slices) {
      if (slice.route_index >= profiles.size()) {
        throw std::out_of_range("incumbent route slice has invalid route");
      }
      const RouteProfile& profile = profiles[slice.route_index];
      if (slice.begin >= slice.end || slice.end > profile.customers.size()) {
        throw std::out_of_range("incumbent route slice has invalid bounds");
      }
      const int first_customer =
          slice.reversed ? profile.customers[slice.end - 1]
                         : profile.customers[slice.begin];
      append_component(edge(previous, first_customer));
      profile.append_slice_components(slice.begin, slice.end, slice.reversed,
                                      components, component_count);
      previous = slice.reversed ? profile.customers[slice.begin]
                                : profile.customers[slice.end - 1];
    }
    append_component(edge(previous, instance.depot().id));
    return component_count;
  }

  RouteCostSummary evaluate_slices(
      std::span<const IncumbentRouteSlice> slices) const {
    std::array<const ChargingSequenceLabel*, kMaximumLabelComponents>
        components;
    const std::size_t count = collect_components(slices, components);
    return evaluator.evaluate(
        std::span<const ChargingSequenceLabel* const>{components.data(),
                                                      count});
  }

  double generalized_cost_slices(
      std::span<const IncumbentRouteSlice> slices,
      double penalty_lambda) const {
    std::array<const ChargingSequenceLabel*, kMaximumLabelComponents>
        components;
    const std::size_t count = collect_components(slices, components);
    return evaluator.generalized_cost(
        std::span<const ChargingSequenceLabel* const>{components.data(),
                                                      count},
        penalty_lambda);
  }
};

RouteSequenceScorer::RouteSequenceScorer(
    const Instance& instance, const Solution& reference,
    const PlanFactory::GapMap& reference_gaps)
    : implementation_(std::make_unique<Implementation>(instance, reference,
                                                        reference_gaps)) {}

RouteSequenceScorer::~RouteSequenceScorer() = default;
RouteSequenceScorer::RouteSequenceScorer(RouteSequenceScorer&&) noexcept =
    default;
RouteSequenceScorer& RouteSequenceScorer::operator=(
    RouteSequenceScorer&&) noexcept = default;

void RouteSequenceScorer::reset(
    const Solution& reference,
    const PlanFactory::GapMap& reference_gaps) {
  implementation_->reset(reference, reference_gaps);
}

RouteCostSummary RouteSequenceScorer::evaluate(
    const std::vector<int>& customer_sequence) const {
  return implementation_->evaluator.evaluate(
      implementation_->label(customer_sequence));
}

RouteCostSummary RouteSequenceScorer::evaluate_insertion(
    std::size_t route_index, std::size_t position, int customer_id) const {
  return implementation_->evaluator.evaluate(
      implementation_->insertion_label(route_index, position, customer_id));
}

double RouteSequenceScorer::generalized_cost(
    const std::vector<int>& customer_sequence, double penalty_lambda) const {
  return implementation_->evaluator.generalized_cost(
      implementation_->label(customer_sequence), penalty_lambda);
}

RouteCostSummary RouteSequenceScorer::evaluate(
    std::span<const IncumbentRouteSlice> slices) const {
  return implementation_->evaluate_slices(slices);
}

double RouteSequenceScorer::generalized_cost(
    std::span<const IncumbentRouteSlice> slices,
    double penalty_lambda) const {
  return implementation_->generalized_cost_slices(slices, penalty_lambda);
}

}  // namespace ils_sp
