#include "ils_sp/copt_models.hpp"
#include "ils_sp/detail/copt_result_validation.hpp"
#include "ils_sp/detail/input_contract.hpp"
#include "callbackbase.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdio>
#include <deque>
#include <exception>
#include <limits>
#include <iostream>
#include <map>
#include <mutex>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

namespace ils_sp {
namespace {

using ModelClock = std::chrono::steady_clock;
constexpr double kDeadlineGuardSeconds = 0.01;

[[nodiscard]] double wall_seconds(ModelClock::time_point start) {
  return std::chrono::duration<double>(ModelClock::now() - start).count();
}

template <typename Result>
[[nodiscard]] Result time_limit_result(ModelClock::time_point start) {
  Result result;
  result.status = COPT_STATUS_TIMEOUT;
  result.status_name = copt_status_name(COPT_STATUS_TIMEOUT);
  result.wall_runtime_s = wall_seconds(start);
  return result;
}

struct ModelTimeBudget {
  bool expired{};
  std::optional<double> applied_time_limit_s;
};

[[nodiscard]] bool deadline_exhausted(CoptDeadline deadline) {
  return deadline.has_value() &&
         std::chrono::duration<double>(*deadline - ModelClock::now()).count() <=
             kDeadlineGuardSeconds;
}

[[nodiscard]] CoptDeadline bounded_deadline(
    ModelClock::time_point start, CoptDeadline deadline,
    std::optional<double> maximum_wall_time_s) {
  if (!maximum_wall_time_s.has_value()) return deadline;
  const auto local_deadline =
      detail::checked_deadline_after<ModelClock>(start, *maximum_wall_time_s);
  if (!deadline.has_value() || local_deadline < *deadline) {
    return local_deadline;
  }
  return deadline;
}

[[nodiscard]] ModelTimeBudget model_time_budget(
    const CoptConfig& config, CoptDeadline deadline) {
  std::optional<double> limit = config.time_limit_s;
  if (deadline.has_value()) {
    const double remaining =
        std::chrono::duration<double>(*deadline - ModelClock::now()).count() -
        kDeadlineGuardSeconds;
    if (remaining <= 0.0) {
      return ModelTimeBudget{.expired = true,
                             .applied_time_limit_s = std::nullopt};
    }
    limit = limit.has_value() ? std::min(*limit, remaining) : remaining;
  }
  return ModelTimeBudget{.expired = false,
                         .applied_time_limit_s = limit};
}

void apply_time_budget(Model& model, const ModelTimeBudget& budget) {
  if (budget.applied_time_limit_s.has_value()) {
    model.SetDblParam(COPT_DBLPARAM_TIMELIMIT, *budget.applied_time_limit_s);
  }
}

[[nodiscard]] CoptModelSize model_size(Model& model) {
  return CoptModelSize{
      .variables = static_cast<std::size_t>(
          model.GetIntAttr(COPT_INTATTR_COLS)),
      .binary_variables = static_cast<std::size_t>(
          model.GetIntAttr(COPT_INTATTR_BINS)),
      .linear_constraints = static_cast<std::size_t>(
          model.GetIntAttr(COPT_INTATTR_ROWS)),
      .sos_constraints = static_cast<std::size_t>(
          model.GetIntAttr(COPT_INTATTR_SOSS)),
      .linear_nonzeros = static_cast<std::size_t>(
          model.GetIntAttr(COPT_INTATTR_ELEMS))};
}

using detail::finite_copt_value;

[[nodiscard]] std::optional<double> finite_optional(double value) {
  return finite_copt_value(value) ? std::optional<double>(value)
                                    : std::nullopt;
}

[[nodiscard]] std::optional<double> relative_gap(
    std::optional<double> incumbent, std::optional<double> bound) {
  if (!incumbent.has_value() || !bound.has_value()) return std::nullopt;
  return std::abs(*incumbent - *bound) /
         (1e-10 + std::abs(*incumbent));
}

[[nodiscard]] bool rounded_target_reached(
    double objective_h, std::optional<double> rounded_target_objective_h) {
  return rounded_target_objective_h.has_value() &&
         detail::roundable_objective(objective_h) &&
         std::llround(objective_h * 100.0) <=
             std::llround(*rounded_target_objective_h * 100.0);
}

struct PiecewiseCurve {
  std::vector<double> knots;
  std::vector<Var> weights;

  void set_start(Model& model, double input, bool active = true) const {
    if (!active) {
      for (const Var& weight : weights) model.SetMipStart(weight, 0.0);
      return;
    }
    std::vector<double> values(knots.size(), 0.0);
    input = std::clamp(input, knots.front(), knots.back());
    auto upper = std::upper_bound(knots.begin(), knots.end(), input);
    if (upper == knots.end()) {
      values.back() = 1.0;
    } else {
      const auto right = static_cast<std::size_t>(upper - knots.begin());
      const auto left = right - 1;
      values[right] = (input - knots[left]) / (knots[right] - knots[left]);
      values[left] = 1.0 - values[right];
    }
    for (std::size_t i = 0; i < weights.size(); ++i) {
      model.SetMipStart(weights[i], values[i]);
    }
  }
};

[[nodiscard]] PiecewiseCurve add_piecewise_curve(
    Model& model, const Var& input, const Var& output,
    const Var& selected, int count,
    const double* x, const double* y, const std::string& name) {
  PiecewiseCurve curve;
  curve.knots.assign(x, x + count);
  VarArray sos;
  Expr sum = 0.0, input_sum = 0.0, output_sum = 0.0;
  for (int i = 0; i < count; ++i) {
    Var weight = model.AddVar(0.0, 1.0, 0.0, COPT_CONTINUOUS,
                             (name + "_weight_" + std::to_string(i)).c_str());
    curve.weights.push_back(weight);
    sos.PushBack(weight);
    sum += weight;
    input_sum += x[i] * weight;
    output_sum += y[i] * weight;
  }
  // Froger (2019) FRVCP (97)-(98), (104)-(105): both PWL graphs
  // are scaled by path selection. Inactive paths have *all* weights zero,
  // not a unit weight at the zero-energy knot. Keep the per-path graph.
  model.AddConstr(sum == selected, (name + "_sum").c_str());
  model.AddConstr(input == input_sum, (name + "_input").c_str());
  model.AddConstr(output == output_sum, (name + "_output").c_str());
  model.AddSos(sos, x, COPT_SOS_TYPE2);
  return curve;
}

struct FrvcpRouteLowerBound {
  double objective_h{};
  double service_time_h{};
};

struct FrvcpStationSpan {
  std::size_t first{};
  std::size_t last{};
};

struct FrvcpPathSpec {
  std::vector<int> stations;
  std::vector<double> leg_energy_wh;
  double travel_time_h{};
  double first_leg_h{};
  double last_leg_h{};
  std::size_t fastest_curve_rank{};
  std::size_t slowest_curve_rank{};
  bool required_for_start{};
  // Model-only spans; external path occurrences and catalog identity stay intact.
  std::vector<FrvcpStationSpan> model_station_spans;
};

struct FrvcpGapPathCatalog {
  int origin{};
  int target{};
  std::vector<FrvcpPathSpec> paths;
};

struct FrvcpPathCatalog {
  bool complete{true};
  std::vector<FrvcpGapPathCatalog> gaps;
  FrvcpNetworkSize size;
};

void saturating_add(std::size_t& value, std::size_t increment) {
  const std::size_t maximum = std::numeric_limits<std::size_t>::max();
  value = increment > maximum - value ? maximum : value + increment;
}

[[nodiscard]] std::unordered_map<int, std::size_t> station_curve_ranks(
    const Instance& instance) {
  std::vector<std::string> curve_types;
  curve_types.reserve(instance.charging_curves().size());
  for (const auto& [type, curve] : instance.charging_curves()) {
    (void)curve;
    curve_types.push_back(type);
  }
  std::sort(curve_types.begin(), curve_types.end());
  for (std::size_t left = 0; left < curve_types.size(); ++left) {
    for (std::size_t right = left + 1; right < curve_types.size(); ++right) {
      const ChargingCurve& left_curve =
          instance.charging_curves().at(curve_types[left]);
      const ChargingCurve& right_curve =
          instance.charging_curves().at(curve_types[right]);
      if (!charging_curve_no_slower(left_curve, right_curve) &&
          !charging_curve_no_slower(right_curve, left_curve)) {
        throw std::invalid_argument(
            "exact TRC path dominance requires comparable charging curves");
      }
    }
  }
  // Pairwise epsilon equivalence is not transitive and cannot be used by
  // std::sort. Compare exact finite signatures on one common energy grid.
  // Tolerance remains only in the separate dominance validation below.
  std::vector<double> energy_breakpoints;
  for (const auto& type : curve_types) {
    for (const ChargingPoint& point : instance.charging_curves().at(type).points()) {
      energy_breakpoints.push_back(point.energy_wh);
    }
  }
  std::sort(energy_breakpoints.begin(), energy_breakpoints.end());
  energy_breakpoints.erase(
      std::unique(energy_breakpoints.begin(), energy_breakpoints.end()),
      energy_breakpoints.end());
  std::map<std::string, std::vector<double>> curve_signatures;
  for (const auto& type : curve_types) {
    const ChargingCurve& curve = instance.charging_curves().at(type);
    auto& signature = curve_signatures[type];
    for (std::size_t index = 1; index < energy_breakpoints.size(); ++index) {
      const double duration = curve.charge_duration(
          energy_breakpoints[index - 1], energy_breakpoints[index]);
      if (!std::isfinite(duration)) {
        throw std::invalid_argument("charging curve ordering requires finite durations");
      }
      signature.push_back(duration);
    }
  }
  std::sort(curve_types.begin(), curve_types.end(),
            [&](const std::string& left, const std::string& right) {
              const auto& first = curve_signatures.at(left);
              const auto& second = curve_signatures.at(right);
              if (first != second) return first < second;
              return left < right;
            });
  // A rank is consumed as a total dominance chain. Validate every pair:
  // adjacent epsilon comparisons alone would not establish transitivity.
  for (std::size_t left = 0; left < curve_types.size(); ++left) {
    for (std::size_t right = left + 1; right < curve_types.size(); ++right) {
      if (!charging_curve_no_slower(
              instance.charging_curves().at(curve_types[left]),
              instance.charging_curves().at(curve_types[right]))) {
        throw std::invalid_argument(
            "exact TRC path dominance requires a consistent charging curve order");
      }
    }
  }
  std::unordered_map<std::string, std::size_t> rank_by_type;
  for (std::size_t rank = 0; rank < curve_types.size(); ++rank) {
    rank_by_type.emplace(curve_types[rank], rank);
  }
  std::unordered_map<int, std::size_t> result;
  result.reserve(instance.station_ids().size());
  for (const int station_id : instance.station_ids()) {
    result.emplace(station_id,
                   rank_by_type.at(instance.node(station_id).station_type));
  }
  return result;
}

[[nodiscard]] bool frvcp_path_dominates(const FrvcpPathSpec& left,
                                        const FrvcpPathSpec& right) {
  if (left.stations.empty() || right.stations.empty()) return false;
  // TRC Section 4.3: even the slowest charger on the dominating RP must
  // be no slower than the fastest charger on the dominated RP.  Montoya
  // has no station waiting time, so the waiting-function condition vanishes.
  return left.slowest_curve_rank <= right.fastest_curve_rank &&
         left.first_leg_h <= right.first_leg_h + kTimeToleranceH &&
         left.last_leg_h <= right.last_leg_h + kTimeToleranceH &&
         left.travel_time_h <= right.travel_time_h + kTimeToleranceH;
}

// Optimistic marginal charging rate across the whole Montoya instance.
// An instantaneous positive-energy segment gives an infinite upper rate.
[[nodiscard]] double maximum_charging_rate(const Instance& instance) {
  double maximum = 0.0;
  for (const int station_id : instance.station_ids()) {
    const auto& points = instance.curve_for_station(station_id).points();
    for (std::size_t index = 1; index < points.size(); ++index) {
      const double energy = points[index].energy_wh - points[index - 1].energy_wh;
      if (energy <= 0.0) continue;
      const double time = points[index].time_h - points[index - 1].time_h;
      if (time <= 1e-15) return std::numeric_limits<double>::infinity();
      maximum = std::max(maximum, energy / time);
    }
  }
  return maximum;
}

[[nodiscard]] FrvcpPathCatalog enumerate_frvcp_paths(
    const Instance& instance, const std::vector<int>& anchors,
    const std::vector<std::vector<int>>* incumbent_gaps,
    double route_objective_upper_bound_h, CoptDeadline deadline) {
  FrvcpPathCatalog result;
  result.size.customer_count = anchors.size() - 2;
  result.size.gap_count = anchors.size() - 1;
  result.gaps.reserve(anchors.size() - 1);
  const auto curve_rank = station_curve_ranks(instance);
  const auto& stations = instance.station_ids();
  const std::size_t station_count = stations.size();
  const double capacity = instance.vehicle().battery_capacity_wh;
  double direct_route_travel_h = 0.0;
  for (std::size_t gap = 0; gap + 1 < anchors.size(); ++gap) {
    direct_route_travel_h +=
        instance.travel_time_h(anchors[gap], anchors[gap + 1]);
  }

  for (std::size_t gap_index = 0; gap_index + 1 < anchors.size();
       ++gap_index) {
    if (deadline_exhausted(deadline)) {
      result.complete = false;
      break;
    }
    FrvcpGapPathCatalog gap{.origin = anchors[gap_index],
                            .target = anchors[gap_index + 1],
                            .paths = {}};
    const double direct_gap_h = instance.travel_time_h(gap.origin, gap.target);
    const double path_travel_upper_h =
        route_objective_upper_bound_h -
        (direct_route_travel_h - direct_gap_h) + kTimeToleranceH;

    const auto make_path = [&](const std::vector<int>& path,
                               bool required_for_start) {
      FrvcpPathSpec spec;
      spec.stations = path;
      spec.required_for_start = required_for_start;
      spec.leg_energy_wh.reserve(path.size() + 1);
      int previous = gap.origin;
      for (std::size_t index = 0; index < path.size(); ++index) {
        const int station_id = path[index];
        const double energy = instance.energy_wh(previous, station_id);
        spec.leg_energy_wh.push_back(energy);
        spec.travel_time_h += instance.travel_time_h(previous, station_id);
        const std::size_t rank = curve_rank.at(station_id);
        if (index == 0) {
          spec.fastest_curve_rank = rank;
          spec.slowest_curve_rank = rank;
          spec.first_leg_h = instance.travel_time_h(previous, station_id);
        } else {
          spec.fastest_curve_rank = std::min(spec.fastest_curve_rank, rank);
          spec.slowest_curve_rank = std::max(spec.slowest_curve_rank, rank);
        }
        previous = station_id;
      }
      const double last_energy = instance.energy_wh(previous, gap.target);
      spec.leg_energy_wh.push_back(last_energy);
      spec.last_leg_h = instance.travel_time_h(previous, gap.target);
      spec.travel_time_h += spec.last_leg_h;
      if (path.empty()) {
        spec.first_leg_h = spec.last_leg_h;
      }
      return spec;
    };

    auto add_path = [&](FrvcpPathSpec candidate) {
      saturating_add(result.size.enumerated_paths, 1);
      saturating_add(result.size.enumerated_station_visits,
                     candidate.stations.size());
      auto duplicate = std::find_if(
          gap.paths.begin(), gap.paths.end(), [&](const FrvcpPathSpec& current) {
            return current.stations == candidate.stations;
          });
      if (duplicate != gap.paths.end()) {
        duplicate->required_for_start =
            duplicate->required_for_start || candidate.required_for_start;
        return;
      }
      const bool dominated = std::any_of(
          gap.paths.begin(), gap.paths.end(), [&](const FrvcpPathSpec& current) {
            return frvcp_path_dominates(current, candidate);
          });
      if (dominated && !candidate.required_for_start) return;
      std::erase_if(gap.paths, [&](const FrvcpPathSpec& current) {
        return !current.required_for_start &&
               frvcp_path_dominates(candidate, current);
      });
      gap.paths.push_back(std::move(candidate));
    };

    if (instance.energy_wh(gap.origin, gap.target) <=
            capacity + kEnergyToleranceWh &&
        direct_gap_h <= path_travel_upper_h) {
      add_path(make_path({}, incumbent_gaps != nullptr &&
                                 (*incumbent_gaps)[gap_index].empty()));
    }
    if (incumbent_gaps != nullptr &&
        !(*incumbent_gaps)[gap_index].empty()) {
      FrvcpPathSpec incumbent =
          make_path((*incumbent_gaps)[gap_index], true);
      const bool individually_reachable =
          std::all_of(incumbent.leg_energy_wh.begin(),
                      incumbent.leg_energy_wh.end(), [&](double energy) {
                        return energy <= capacity + kEnergyToleranceWh;
                      });
      if (!individually_reachable) {
        throw std::invalid_argument(
            "FRVCP incumbent contains a battery-infeasible station leg");
      }
      add_path(std::move(incumbent));
    }

    // A reverse shortest-path relaxation supplies an admissible completion
    // bound.  CSPs are elementary as in Froger et al. (2019); every complete
    // elementary path under the incumbent/duration cutoff is enumerated before
    // the paper's dominance rule is applied.
    const double infinity = std::numeric_limits<double>::infinity();
    std::vector<double> shortest_to_target(station_count, infinity);
    using DistanceNode = std::pair<double, std::size_t>;
    std::priority_queue<DistanceNode, std::vector<DistanceNode>,
                        std::greater<DistanceNode>>
        queue;
    for (std::size_t index = 0; index < station_count; ++index) {
      if (instance.energy_wh(stations[index], gap.target) <=
          capacity + kEnergyToleranceWh) {
        shortest_to_target[index] =
            instance.travel_time_h(stations[index], gap.target);
        queue.emplace(shortest_to_target[index], index);
      }
    }
    while (!queue.empty()) {
      const auto [distance, next] = queue.top();
      queue.pop();
      if (distance != shortest_to_target[next]) continue;
      for (std::size_t previous = 0; previous < station_count; ++previous) {
        if (previous == next ||
            instance.energy_wh(stations[previous], stations[next]) >
                capacity + kEnergyToleranceWh) {
          continue;
        }
        const double candidate =
            instance.travel_time_h(stations[previous], stations[next]) +
            distance;
        if (candidate + kTimeToleranceH < shortest_to_target[previous]) {
          shortest_to_target[previous] = candidate;
          queue.emplace(candidate, previous);
        }
      }
    }

    std::vector<std::vector<std::size_t>> successors(station_count + 1);
    const std::size_t origin_index = station_count;
    for (std::size_t previous = 0; previous <= station_count; ++previous) {
      const int previous_id =
          previous == origin_index ? gap.origin : stations[previous];
      auto& next_nodes = successors[previous];
      for (std::size_t next = 0; next < station_count; ++next) {
        if (previous == next || !std::isfinite(shortest_to_target[next]) ||
            instance.energy_wh(previous_id, stations[next]) >
                capacity + kEnergyToleranceWh) {
          continue;
        }
        next_nodes.push_back(next);
      }
      std::sort(next_nodes.begin(), next_nodes.end(), [&](std::size_t left,
                                                          std::size_t right) {
        const double left_bound =
            instance.travel_time_h(previous_id, stations[left]) +
            shortest_to_target[left];
        const double right_bound =
            instance.travel_time_h(previous_id, stations[right]) +
            shortest_to_target[right];
        if (left_bound != right_bound) return left_bound < right_bound;
        return stations[left] < stations[right];
      });
    }

    std::vector<bool> used(station_count, false);
    std::vector<int> path;
    std::size_t extensions = 0;
    const auto enumerate = [&](auto&& self, std::size_t previous,
                               double prefix_travel_h) -> void {
      if (!result.complete) return;
      if ((extensions++ & 1023U) == 0U && deadline_exhausted(deadline)) {
        result.complete = false;
        return;
      }
      const int previous_id =
          previous == origin_index ? gap.origin : stations[previous];
      for (const std::size_t next : successors[previous]) {
        if (used[next]) continue;
        const double next_prefix =
            prefix_travel_h +
            instance.travel_time_h(previous_id, stations[next]);
        if (next_prefix + shortest_to_target[next] >
            path_travel_upper_h + kTimeToleranceH) {
          continue;
        }
        used[next] = true;
        path.push_back(stations[next]);
        if (instance.energy_wh(stations[next], gap.target) <=
            capacity + kEnergyToleranceWh) {
          const double total = next_prefix +
                               instance.travel_time_h(stations[next],
                                                      gap.target);
          if (total <= path_travel_upper_h + kTimeToleranceH) {
            add_path(make_path(path, false));
          }
        }
        self(self, next, next_prefix);
        path.pop_back();
        used[next] = false;
        if (!result.complete) return;
      }
    };
    enumerate(enumerate, origin_index, 0.0);
    if (!result.complete) break;

    std::sort(gap.paths.begin(), gap.paths.end(),
              [](const FrvcpPathSpec& left, const FrvcpPathSpec& right) {
                return std::tuple{left.travel_time_h, left.stations.size(),
                                  left.stations} <
                       std::tuple{right.travel_time_h, right.stations.size(),
                                  right.stations};
              });
    result.gaps.push_back(std::move(gap));
  }

  // TRC 4.2 / Froger 5.1: a noncharging customer couples the legs before
  // and after it. Apply necessary feasibility conditions AFTER the existing
  // dominance pass; do not revive paths discarded by the paper's rule.
  if (result.complete) {
    const int depot = instance.depot().id;
    for (FrvcpGapPathCatalog& gap : result.gaps) {
      double into_origin = instance.energy_wh(depot, gap.origin);
      double out_of_target = instance.energy_wh(gap.target, depot);
      for (const int station : stations) {
        into_origin = std::min(into_origin, instance.energy_wh(station, gap.origin));
        out_of_target = std::min(out_of_target, instance.energy_wh(gap.target, station));
      }
      std::erase_if(gap.paths, [&](const FrvcpPathSpec& path_spec) {
        if (path_spec.required_for_start) return false;
        if (path_spec.stations.empty()) {
          return into_origin + path_spec.leg_energy_wh.front() + out_of_target >
                 capacity + kEnergyToleranceWh;
        }
        return into_origin + path_spec.leg_energy_wh.front() > capacity + kEnergyToleranceWh ||
               path_spec.leg_energy_wh.back() + out_of_target > capacity + kEnergyToleranceWh;
      });
    }

    // A path must also fit into a full route, even if every other gap takes
    // its shortest surviving alternative and all charging uses the fastest
    // possible rate. This is a hard Tmax bound, not an incumbent-quality cap.
    const bool nonempty = std::all_of(result.gaps.begin(), result.gaps.end(),
        [](const auto& gap) { return !gap.paths.empty(); });
    if (nonempty) {
      std::vector<double> minimum_travel, minimum_energy;
      double total_travel = 0.0, total_energy = 0.0, service = 0.0;
      for (std::size_t i = 1; i + 1 < anchors.size(); ++i) {
        service += instance.node(anchors[i]).service_time_h;
      }
      for (const auto& gap : result.gaps) {
        double travel = std::numeric_limits<double>::infinity();
        double energy = std::numeric_limits<double>::infinity();
        for (const auto& path_spec : gap.paths) {
          travel = std::min(travel, path_spec.travel_time_h);
          energy = std::min(energy, std::accumulate(path_spec.leg_energy_wh.begin(),
                                                   path_spec.leg_energy_wh.end(), 0.0));
        }
        minimum_travel.push_back(travel);
        minimum_energy.push_back(energy);
        total_travel += travel;
        total_energy += energy;
      }
      const double maximum_rate = maximum_charging_rate(instance);
      for (std::size_t i = 0; i < result.gaps.size(); ++i) {
        std::erase_if(result.gaps[i].paths, [&](const FrvcpPathSpec& path_spec) {
          if (path_spec.required_for_start) return false;
          const double travel = total_travel - minimum_travel[i] + path_spec.travel_time_h;
          const double energy = total_energy - minimum_energy[i] +
              std::accumulate(path_spec.leg_energy_wh.begin(), path_spec.leg_energy_wh.end(), 0.0);
          const double required = std::max(0.0, energy - capacity);
          const double charging = required <= kEnergyToleranceWh ? 0.0 :
              maximum_rate > 0.0 ? required / maximum_rate :
              std::numeric_limits<double>::infinity();
          return travel + charging + service >
                 instance.vehicle().max_route_duration_h + kTimeToleranceH;
        });
      }
    }
  }
  // Enumerated counts retain the original work; active counts describe only
  // the final model. Stable erasure preserves surviving path/tie order.
  for (const auto& gap : result.gaps) {
    saturating_add(result.size.active_paths, gap.paths.size());
    for (const auto& path_spec : gap.paths) {
      saturating_add(result.size.active_station_visits, path_spec.stations.size());
    }
  }
  result.size.candidate_station_copies =
      result.size.enumerated_station_visits;
  result.size.active_station_copies = result.size.active_station_visits;
  return result;
}

// For a fixed customer order, Euclidean direct anchor legs minimize both
// travel time and energy.  Any feasible station detour therefore consumes at
// least the direct energy.  Energy consumed beyond the initial battery must be
// recharged, and even the fastest marginal segment bounds that charging time.
[[nodiscard]] FrvcpRouteLowerBound frvcp_route_lower_bound(
    const Instance& instance, const std::vector<int>& customer_ids) {
  std::vector<int> anchors;
  anchors.reserve(customer_ids.size() + 2);
  anchors.push_back(instance.depot().id);
  anchors.insert(anchors.end(), customer_ids.begin(), customer_ids.end());
  anchors.push_back(instance.depot().id);

  double direct_travel_h = 0.0;
  double direct_energy_wh = 0.0;
  for (std::size_t index = 0; index + 1 < anchors.size(); ++index) {
    direct_travel_h += instance.travel_time_h(anchors[index],
                                               anchors[index + 1]);
    direct_energy_wh += instance.energy_wh(anchors[index],
                                            anchors[index + 1]);
  }
  double service_time_h = 0.0;
  for (const int customer_id : customer_ids) {
    service_time_h += instance.node(customer_id).service_time_h;
  }

  const double required_charge_wh = std::max(
      0.0, direct_energy_wh - instance.vehicle().battery_capacity_wh);
  if (required_charge_wh <= kEnergyToleranceWh) {
    return FrvcpRouteLowerBound{.objective_h = direct_travel_h,
                                .service_time_h = service_time_h};
  }

  const double maximum_rate_wh_per_h = maximum_charging_rate(instance);
  const double charging_lower_bound_h = maximum_rate_wh_per_h > 0.0
      ? required_charge_wh / maximum_rate_wh_per_h
      : std::numeric_limits<double>::infinity();
  return FrvcpRouteLowerBound{
      .objective_h = direct_travel_h + charging_lower_bound_h,
      .service_time_h = service_time_h};
}

class SetPartitioningCallback final : public CallbackBase {
 public:
  SetPartitioningCallback(SetPartitioningProgressCallback progress,
                          std::optional<double> target,
                          std::optional<double> stall_limit,
                          std::optional<double> initial_incumbent,
                          double interval, CoptModelSize size, double build_s,
                          bool serialize_callbacks)
      : progress_(std::move(progress)), target_(target), stall_limit_(stall_limit),
        best_(initial_incumbent), last_improvement_(initial_incumbent.has_value()
            ? std::optional<double>(0.0) : std::nullopt), interval_(interval),
        next_progress_(interval), size_(size), build_s_(build_s),
        serialize_callbacks_(serialize_callbacks) {}

  void start_clock() { started_ = ModelClock::now(); }
  void rethrow_callback_error() const {
    if (error_) std::rethrow_exception(error_);
  }
  [[nodiscard]] bool incumbent_stall_limit_reached() const noexcept {
    return stalled_;
  }

 private:
  void callback() override {
    // Protect our own mutable observer state when COPT uses multiple threads.
    // Callback context validity remains the responsibility of the public API.
    std::unique_lock<std::mutex> lock(callback_mutex_, std::defer_lock);
    if (serialize_callbacks_) lock.lock();
    if (error_) return;
    try {
      const double runtime = wall_seconds(started_);
      // With no stopping policy, callback values only feed periodic reporting.
      // Avoid querying COPT at every node/relaxation between reports. Target and
      // stall policies still inspect every callback, including improvements.
      if (!target_ && !stall_limit_ &&
          (!progress_ || runtime + 1e-9 < next_progress_)) {
        return;
      }
      const bool has_incumbent = GetIntInfo(COPT_CBINFO_HASINCUMBENT) != 0;
      const auto incumbent = has_incumbent
          ? finite_optional(GetDblInfo(COPT_CBINFO_BESTOBJ)) : std::nullopt;
      if (incumbent && (!best_ || *incumbent < *best_ - kCostTolerance)) {
        best_ = incumbent;
        last_improvement_ = runtime;
      }
      const bool reached = incumbent && rounded_target_reached(*incumbent, target_);
      const bool stalled = stall_limit_ && last_improvement_ &&
          runtime - *last_improvement_ + 1e-9 >= *stall_limit_;
      if (progress_ && (reached || stalled || runtime + 1e-9 >= next_progress_)) {
        const auto bound = finite_optional(GetDblInfo(COPT_CBINFO_BESTBND));
        progress_(SetPartitioningProgress{
            .event = SetPartitioningProgressEvent::Search,
            .runtime_s = runtime, .model_build_runtime_s = build_s_,
            .incumbent_objective_h = incumbent, .best_bound_h = bound,
            .relative_gap = relative_gap(incumbent, bound),
            // COPT does not expose node counts through callback information.
            .explored_nodes = std::nullopt,
            .solution_count = has_incumbent ? 1 : 0,
            .rounded_target_reached = reached,
            .incumbent_stall_limit_reached = stalled, .model_size = size_});
        next_progress_ = (std::floor(runtime / interval_) + 1.0) * interval_;
      }
      if (reached || stalled) {
        stalled_ = stalled;
        Interrupt();
      }
    } catch (...) {
      error_ = std::current_exception();
      Interrupt();
    }
  }
  SetPartitioningProgressCallback progress_;
  std::optional<double> target_, stall_limit_, best_, last_improvement_;
  double interval_{}, next_progress_{};
  CoptModelSize size_;
  double build_s_{};
  bool serialize_callbacks_{};
  std::mutex callback_mutex_;
  ModelClock::time_point started_{};
  bool stalled_{};
  std::exception_ptr error_;
};

}  // namespace

std::string copt_status_name(int status) {
  switch (status) {
    case COPT_STATUS_UNSTARTED: return "UNSTARTED";
    case COPT_STATUS_OPTIMAL: return "OPTIMAL";
    case COPT_STATUS_INFEASIBLE: return "INFEASIBLE";
    case COPT_STATUS_UNBOUNDED: return "UNBOUNDED";
    case COPT_STATUS_INF_OR_UNB: return "INF_OR_UNBD";
    case COPT_STATUS_NUMERICAL: return "NUMERICAL";
    case COPT_STATUS_NODELIMIT: return "NODE_LIMIT";
    case COPT_STATUS_IMPRECISE: return "IMPRECISE";
    case COPT_STATUS_TIMEOUT: return "TIME_LIMIT";
    case COPT_STATUS_UNFINISHED: return "UNFINISHED";
    case COPT_STATUS_INTERRUPTED: return "INTERRUPTED";
    case COPT_STATUS_ITERLIMIT: return "ITERATION_LIMIT";
    case COPT_STATUS_MEMLIMIT: return "MEMORY_LIMIT";
    default: return "STATUS_" + std::to_string(status);
  }
}

CoptSetPartitioning::CoptSetPartitioning(const Instance& instance,
                                             Envr& environment,
                                             CoptConfig config)
    : instance_(instance), environment_(environment), config_(config) {
  if (config_.threads < 1 ||
      (config_.sp_integer_threads &&
       (*config_.sp_integer_threads < 1 || *config_.sp_integer_threads > 128)) ||
      !std::isfinite(config_.mip_gap) ||
      config_.mip_gap < 0.0 ||
      (config_.time_limit_s.has_value() &&
       (!std::isfinite(*config_.time_limit_s) ||
        *config_.time_limit_s <= 0.0))) {
    throw std::invalid_argument("invalid Copt configuration");
  }
}

void CoptSetPartitioning::configure(Model& model, bool integer) const {
  model.SetIntParam(COPT_INTPARAM_LOGGING, config_.output_flag);
  model.SetIntParam(COPT_INTPARAM_THREADS,
                    integer ? config_.sp_integer_threads.value_or(config_.threads)
                            : config_.threads);
  if (integer) model.SetDblParam(COPT_DBLPARAM_RELGAP, config_.mip_gap);
}

std::vector<std::vector<std::size_t>>
CoptSetPartitioning::customer_incidence(
    const std::vector<RouteColumn>& columns) const {
  if (columns.empty()) {
    throw std::invalid_argument("set partitioning needs at least one column");
  }
  const std::vector<int>& customers = instance_.customer_ids();
  std::unordered_map<int, std::size_t> row_by_customer;
  row_by_customer.reserve(customers.size());
  for (std::size_t row = 0; row < customers.size(); ++row) {
    row_by_customer.emplace(customers[row], row);
  }

  std::vector<std::vector<std::size_t>> incidence(customers.size());
  for (std::size_t column_index = 0; column_index < columns.size();
       ++column_index) {
    if (!finite_copt_value(columns[column_index].cost_h)) {
      throw std::invalid_argument("route column cost must be finite");
    }
    for (const int customer_id : columns[column_index].customer_set) {
      const auto row = row_by_customer.find(customer_id);
      if (row == row_by_customer.end()) {
        throw std::invalid_argument("route column contains unknown customer " +
                                    std::to_string(customer_id));
      }
      incidence[row->second].push_back(column_index);
    }
  }
  for (std::size_t row = 0; row < customers.size(); ++row) {
    if (incidence[row].empty()) {
      throw std::invalid_argument("route pool does not cover customer " +
                                  std::to_string(customers[row]));
    }
  }
  return incidence;
}

SetPartitioningResult CoptSetPartitioning::solve_integer(
    const std::vector<RouteColumn>& columns,
    SetPartitioningSolveOptions options) {
  const auto call_started = ModelClock::now();
  if (options.rounded_target_objective_h) {
    detail::validate_rounded_target(*options.rounded_target_objective_h);
  }
  if (!std::isfinite(options.progress_interval_s) ||
      options.progress_interval_s <= 0.0 ||
      (options.rounded_target_objective_h.has_value() &&
       (!std::isfinite(*options.rounded_target_objective_h) ||
        *options.rounded_target_objective_h < 0.0)) ||
      (options.incumbent_stall_limit_s.has_value() &&
       (!std::isfinite(*options.incumbent_stall_limit_s) ||
        *options.incumbent_stall_limit_s <= 0.0)) ||
      (options.maximum_wall_time_s.has_value() &&
       (!std::isfinite(*options.maximum_wall_time_s) ||
        *options.maximum_wall_time_s <= 0.0))) {
    throw std::invalid_argument("invalid set-partitioning solve options");
  }
  const CoptDeadline deadline = bounded_deadline(
      call_started, options.deadline, options.maximum_wall_time_s);
  const std::vector<std::vector<std::size_t>> incidence =
      customer_incidence(columns);
  std::vector<bool> selected_start(columns.size(), false);
  std::optional<double> mip_start_objective_h;
  if (!options.mip_start_indices.empty()) {
    for (const std::size_t index : options.mip_start_indices) {
      if (index >= columns.size() || selected_start[index]) {
        throw std::invalid_argument("invalid duplicate SP MIP-start column");
      }
      selected_start[index] = true;
    }
    for (const std::vector<std::size_t>& covering_columns : incidence) {
      const std::size_t count = static_cast<std::size_t>(std::count_if(
          covering_columns.begin(), covering_columns.end(),
          [&](std::size_t index) { return selected_start[index]; }));
      if (count != 1) {
        throw std::invalid_argument(
            "SP MIP start is not an exact customer partition");
      }
    }
    double objective_h = 0.0;
    Solution start_solution;
    for (std::size_t index = 0; index < columns.size(); ++index) {
      if (selected_start[index]) {
        const Plan& plan = columns[index].plan;
        if (!plan.evaluation || !plan.evaluation->feasible) {
          throw std::invalid_argument("SP MIP start contains an infeasible plan");
        }
        start_solution.plans.push_back(plan);
        objective_h += columns[index].cost_h;
      }
    }
    start_solution.validate_partition(instance_);
    if (!finite_copt_value(objective_h)) {
      throw std::invalid_argument("SP MIP-start objective must be finite");
    }
    mip_start_objective_h = objective_h;
  }

  // Preserve every budget/early-stop entry's original solve order. Paper-default
  // calls may solve the already-required LP first, and feed it back only later.
  std::optional<LinearRelaxationResult> reused_relaxation;
  detail::SpDualCertificate certificate;
  bool lp_prepass_attempted = false;
  bool lp_prepass_exception = false;
  double lp_prepass_wall_s = 0.0;
  if (options.allow_lp_certificate && instance_.customer_ids().size() > 10 &&
      mip_start_objective_h &&
      !deadline && !config_.time_limit_s &&
      !options.rounded_target_objective_h && !options.incumbent_stall_limit_s &&
      !options.maximum_wall_time_s) {
    lp_prepass_attempted = true;
    const auto lp_started = ModelClock::now();
    try {
      reused_relaxation = solve_relaxation(columns, std::nullopt, true);
    } catch (const CoptException&) {
      // Do not let a speculative LP backend failure prevent the original MIP.
      // The ordinary LP feedback position retries using the original path.
      lp_prepass_exception = true;
    }
    lp_prepass_wall_s = wall_seconds(lp_started);
    if (reused_relaxation && reused_relaxation->status == COPT_STATUS_OPTIMAL) {
      std::vector<double> costs;
      costs.reserve(columns.size());
      for (const auto& column : columns) costs.push_back(column.cost_h);
      std::vector<double> duals;
      duals.reserve(instance_.customer_ids().size());
      for (int customer : instance_.customer_ids()) {
        const auto it = reused_relaxation->duals.find(customer);
        duals.push_back(it == reused_relaxation->duals.end()
            ? std::numeric_limits<double>::quiet_NaN() : it->second);
      }
      certificate = detail::certify_sp_start(costs, incidence,
          options.mip_start_indices, duals,
          reused_relaxation->effective_abs_gap, config_.mip_gap, true);
    }
  }

  const auto retained_start = [&](int status) {
    SetPartitioningResult result;
    result.status = status;
    result.status_name = copt_status_name(status);
    result.wall_runtime_s = wall_seconds(call_started);
    result.mip_start_routes = options.mip_start_indices.size();
    result.mip_start_objective_h = mip_start_objective_h;
    result.lp_prepass_attempted = lp_prepass_attempted;
    result.lp_prepass_exception = lp_prepass_exception;
    result.lp_prepass_wall_s = lp_prepass_wall_s;
    result.lp_certificate = certificate;
    result.reused_relaxation = reused_relaxation;
    // This is the validated application incumbent, not a backend solution.
    // Keep solution_count=0 when COPT has not supplied a usable solution.
    if (mip_start_objective_h) {
      result.objective_h = mip_start_objective_h;
      for (const std::size_t index : options.mip_start_indices) {
        result.selected_columns.push_back(columns[index]);
      }
      result.rounded_target_reached = rounded_target_reached(
          *mip_start_objective_h, options.rounded_target_objective_h);
    }
    return result;
  };
  if (certificate.certified) {
    SetPartitioningResult result = retained_start(COPT_STATUS_OPTIMAL);
    result.status_name = "LP_CERTIFIED_START";
    result.lp_certified_start = true;
    result.best_bound_h = certificate.lower_h;
    if (certificate.lower_h >= 0.0 && certificate.upper_h > 0.0)
      result.relative_gap = certificate.relative_gap_upper;
    // Match the MIP decoding's column order, not the supplied start's order.
    result.selected_columns.clear();
    for (std::size_t j = 0; j < columns.size(); ++j)
      if (selected_start[j]) result.selected_columns.push_back(columns[j]);
    return result;
  }
  if (model_time_budget(config_, deadline).expired) {
    return retained_start(COPT_STATUS_TIMEOUT);
  }
  Model model = environment_.CreateModel(("ils_sp_partition_" + instance_.name()).c_str());
  configure(model, true);
  std::vector<Var> theta;
  theta.reserve(columns.size());
  std::size_t next_fixed_zero = 0;
  for (std::size_t index = 0; index < columns.size(); ++index) {
    if ((index & 255U) == 0U && deadline_exhausted(deadline)) {
      return retained_start(COPT_STATUS_TIMEOUT);
    }
    const bool fixed_zero = next_fixed_zero < certificate.fixed_zero_columns.size() &&
        certificate.fixed_zero_columns[next_fixed_zero] == index;
    next_fixed_zero += fixed_zero;
    theta.push_back(model.AddVar(0.0, fixed_zero ? 0.0 : 1.0, columns[index].cost_h, COPT_BINARY, ("theta_" + std::to_string(index)).c_str()));
  }
  if (deadline_exhausted(deadline)) {
    return retained_start(COPT_STATUS_TIMEOUT);
  }
  const std::vector<int>& customers = instance_.customer_ids();
  for (std::size_t row = 0; row < customers.size(); ++row) {
    if (deadline_exhausted(deadline)) {
      return retained_start(COPT_STATUS_TIMEOUT);
    }
    Expr cover = 0.0;
    for (const std::size_t column_index : incidence[row]) {
      cover += theta[column_index];
    }
    model.AddConstr(cover == 1.0, ("cover_" + std::to_string(customers[row])).c_str());
  }

  if (mip_start_objective_h) {
    for (std::size_t index = 0; index < columns.size(); ++index) {
      model.SetMipStart(theta[index], selected_start[index] ? 1.0 : 0.0);
    }
  }

  model.SetObjSense(COPT_MINIMIZE);
  const CoptModelSize size = model_size(model);
  const double model_build_runtime_s =
      std::max(0.0, wall_seconds(call_started) - lp_prepass_wall_s);
  if (options.progress) {
    options.progress(SetPartitioningProgress{
        .event = SetPartitioningProgressEvent::ModelBuilt,
        .runtime_s = 0.0,
        .model_build_runtime_s = model_build_runtime_s,
        .incumbent_objective_h = std::nullopt,
        .best_bound_h = std::nullopt,
        .relative_gap = std::nullopt,
        .explored_nodes = 0.0,
        .solution_count = 0,
        .rounded_target_reached = false,
        .incumbent_stall_limit_reached = false,
        .model_size = size});
  }
  ModelTimeBudget budget = model_time_budget(config_, deadline);
  if (budget.expired) {
    SetPartitioningResult result =
        retained_start(COPT_STATUS_TIMEOUT);
    result.model_build_runtime_s = model_build_runtime_s;
    result.mip_start_routes = options.mip_start_indices.size();
    result.mip_start_objective_h = mip_start_objective_h;
    result.model_size = size;
    return result;
  }
  apply_time_budget(model, budget);
  SetPartitioningCallback callback(
      options.progress, options.rounded_target_objective_h,
      options.incumbent_stall_limit_s, mip_start_objective_h,
      options.progress_interval_s, size, model_build_runtime_s,
      config_.sp_integer_threads.value_or(config_.threads) > 1);
  if (options.progress || options.incumbent_stall_limit_s.has_value() ||
      options.rounded_target_objective_h.has_value()) {
    model.SetCallback(&callback, COPT_CBCONTEXT_INCUMBENT |
                      COPT_CBCONTEXT_MIPNODE | COPT_CBCONTEXT_MIPRELAX);
  }
  if (mip_start_objective_h) model.LoadMipStart();
  budget = model_time_budget(config_, deadline);
  if (budget.expired) {
    SetPartitioningResult result = retained_start(COPT_STATUS_TIMEOUT);
    result.model_build_runtime_s = model_build_runtime_s;
    result.model_size = size;
    return result;
  }
  apply_time_budget(model, budget);
  callback.start_clock();
  model.Solve();
  callback.rethrow_callback_error();
  const int status = model.GetIntAttr(COPT_INTATTR_STATUS);
  const int solution_count = model.GetIntAttr(COPT_INTATTR_HASSOL);
  const std::optional<double> objective_h =
      solution_count > 0
          ? std::optional<double>(model.GetDblAttr(COPT_DBLATTR_OBJVAL))
          : std::nullopt;
  const bool invalid_objective = objective_h && !finite_copt_value(*objective_h);
  const bool stopped_without_solution = solution_count == 0 &&
      (status == COPT_STATUS_TIMEOUT || status == COPT_STATUS_INTERRUPTED ||
       status == COPT_STATUS_IMPRECISE || status == COPT_STATUS_NUMERICAL);
  if (invalid_objective || stopped_without_solution) {
    SetPartitioningResult result = retained_start(
        invalid_objective ? COPT_STATUS_NUMERICAL : status);
    result.runtime_s = model.GetDblAttr(COPT_DBLATTR_SOLVINGTIME);
    result.model_build_runtime_s = model_build_runtime_s;
    result.applied_time_limit_s = budget.applied_time_limit_s;
    result.model_size = size;
    return result;
  }
  const std::optional<double> best_bound_h =
      finite_optional(model.GetDblAttr(COPT_DBLATTR_BESTBND));
  std::vector<RouteColumn> selected;
  if (solution_count > 0) {
    for (std::size_t index = 0; index < columns.size(); ++index) {
      if (theta[index].Get(COPT_DBLINFO_VALUE) > 0.5) {
        selected.push_back(columns[index]);
      }
    }
    Solution assembled;
    for (const RouteColumn& column : selected) {
      assembled.plans.push_back(column.plan);
    }
    assembled.validate_partition(instance_);
  }
  return SetPartitioningResult{
      .status = status,
      .status_name = copt_status_name(status),
      .selected_columns = std::move(selected),
      .objective_h = objective_h,
      .runtime_s = model.GetDblAttr(COPT_DBLATTR_SOLVINGTIME),
      .wall_runtime_s = wall_seconds(call_started),
      .model_build_runtime_s = model_build_runtime_s,
      .applied_time_limit_s = budget.applied_time_limit_s,
      .best_bound_h = best_bound_h,
      .relative_gap = relative_gap(objective_h, best_bound_h),
      .explored_nodes = static_cast<double>(model.GetIntAttr(COPT_INTATTR_NODECNT)),
      .solution_count = solution_count,
      .mip_start_routes = options.mip_start_indices.size(),
      .mip_start_objective_h = mip_start_objective_h,
      .rounded_target_reached =
          objective_h.has_value() &&
          rounded_target_reached(*objective_h,
                                 options.rounded_target_objective_h),
      .incumbent_stall_limit_reached =
          callback.incumbent_stall_limit_reached(),
      .model_size = size,
      .lp_prepass_attempted = lp_prepass_attempted,
      .lp_prepass_exception = lp_prepass_exception,
      .lp_prepass_wall_s = lp_prepass_wall_s,
      .lp_certificate = certificate,
      .reused_relaxation = std::move(reused_relaxation)};
}

LinearRelaxationResult CoptSetPartitioning::solve_relaxation(
    const std::vector<RouteColumn>& columns, CoptDeadline deadline,
    bool read_certificate_abs_gap) {
  const auto call_started = ModelClock::now();
  const std::vector<std::vector<std::size_t>> incidence =
      customer_incidence(columns);
  if (model_time_budget(config_, deadline).expired) {
    return time_limit_result<LinearRelaxationResult>(call_started);
  }
  Model model = environment_.CreateModel(("ils_sp_partition_lp_" + instance_.name()).c_str());
  configure(model, false);
  std::vector<Var> theta;
  theta.reserve(columns.size());
  for (std::size_t index = 0; index < columns.size(); ++index) {
    if ((index & 255U) == 0U && deadline_exhausted(deadline)) {
      return time_limit_result<LinearRelaxationResult>(call_started);
    }
    theta.push_back(model.AddVar(0.0, 1.0, columns[index].cost_h, COPT_CONTINUOUS, ("theta_" + std::to_string(index)).c_str()));
  }
  if (deadline_exhausted(deadline)) {
    return time_limit_result<LinearRelaxationResult>(call_started);
  }
  std::unordered_map<int, Constraint> cover_constraints;
  const std::vector<int>& customers = instance_.customer_ids();
  for (std::size_t row = 0; row < customers.size(); ++row) {
    if (deadline_exhausted(deadline)) {
      return time_limit_result<LinearRelaxationResult>(call_started);
    }
    Expr cover = 0.0;
    for (const std::size_t column_index : incidence[row]) {
      cover += theta[column_index];
    }
    cover_constraints.emplace(
        customers[row],
        model.AddConstr(cover == 1.0, ("cover_" + std::to_string(customers[row])).c_str()));
  }
  model.SetObjSense(COPT_MINIMIZE);
  const CoptModelSize size = model_size(model);
  const ModelTimeBudget budget = model_time_budget(config_, deadline);
  if (budget.expired) {
    LinearRelaxationResult result =
        time_limit_result<LinearRelaxationResult>(call_started);
    result.model_size = size;
    return result;
  }
  apply_time_budget(model, budget);
  model.Solve();
  int status = model.GetIntAttr(COPT_INTATTR_STATUS);
  std::unordered_map<int, double> duals;
  std::optional<double> objective_h;
  if (status == COPT_STATUS_OPTIMAL) {
    objective_h = finite_optional(model.GetDblAttr(COPT_DBLATTR_OBJVAL));
    for (auto& [customer_id, constraint] : cover_constraints) {
      const double dual = constraint.Get(COPT_DBLINFO_DUAL);
      if (!finite_copt_value(dual)) {
        status = COPT_STATUS_NUMERICAL;
        break;
      }
      duals.emplace(customer_id, dual);
    }
    if (!objective_h || status != COPT_STATUS_OPTIMAL) {
      status = COPT_STATUS_NUMERICAL;
      objective_h.reset();
      duals.clear();
    }
  }
  return LinearRelaxationResult{
      .status = status,
      .status_name = copt_status_name(status),
      .objective_h = objective_h,
      .duals = std::move(duals),
      .runtime_s = model.GetDblAttr(COPT_DBLATTR_SOLVINGTIME),
      .wall_runtime_s = wall_seconds(call_started),
      .applied_time_limit_s = budget.applied_time_limit_s,
      .model_size = size,
      .effective_abs_gap = read_certificate_abs_gap
          ? model.GetDblParam(COPT_DBLPARAM_ABSGAP) : 0.0};
}

CoptFrvcp::CoptFrvcp(const Instance& instance, RouteEvaluator& evaluator,
                         PlanFactory& factory, Envr& environment,
                         CoptConfig config,
                         std::size_t infeasible_cache_capacity)
    : instance_(instance),
      evaluator_(evaluator),
      factory_(factory),
      environment_(environment),
      config_(config),
      model_cache_(infeasible_cache_capacity) {
  if (config_.threads < 1 || !std::isfinite(config_.mip_gap) ||
      config_.mip_gap < 0.0 ||
      (config_.time_limit_s.has_value() &&
       (!std::isfinite(*config_.time_limit_s) ||
        *config_.time_limit_s <= 0.0))) {
    throw std::invalid_argument("invalid FRVCP Copt configuration");
  }
}

void CoptFrvcp::configure(Model& model) const {
  model.SetIntParam(COPT_INTPARAM_LOGGING, config_.output_flag);
  model.SetIntParam(COPT_INTPARAM_THREADS, config_.threads);
  model.SetDblParam(COPT_DBLPARAM_RELGAP, config_.mip_gap);
}

FrvcpResult CoptFrvcp::optimize(
    const std::vector<int>& customer_ids, CoptDeadline deadline,
    FrvcpProgressCallback progress) {
  return optimize_impl(customer_ids, nullptr, deadline, progress);
}

FrvcpResult CoptFrvcp::optimize(const Plan& incumbent_plan,
                                  CoptDeadline deadline,
                                  FrvcpProgressCallback progress) {
  return optimize_impl(incumbent_plan.customer_ids, &incumbent_plan,
                       deadline, progress);
}

FrvcpResult CoptFrvcp::optimize_impl(
    const std::vector<int>& customer_ids, const Plan* incumbent_plan,
    CoptDeadline deadline, const FrvcpProgressCallback& progress) {
  const auto call_started = ModelClock::now();
  if (customer_ids.empty() ||
      sorted_customer_set(customer_ids).size() != customer_ids.size()) {
    throw std::invalid_argument("FRVCP requires a nonempty simple customer order");
  }
  if (incumbent_plan != nullptr &&
      (incumbent_plan->customer_ids != customer_ids ||
       incumbent_plan->gaps.size() != customer_ids.size() + 1 ||
       incumbent_plan->route !=
           expand_route(customer_ids, incumbent_plan->gaps) ||
       factory_.has_adjacent_station_visits(incumbent_plan->route.visits) ||
       incumbent_plan->evaluation == nullptr ||
       !incumbent_plan->evaluation->feasible ||
       !finite_copt_value(incumbent_plan->evaluation->raw_cost_h))) {
    throw std::invalid_argument("FRVCP MIP start must be a feasible plan for "
                                "the same customer order");
  }
  const FrvcpRouteLowerBound lower_bound =
      frvcp_route_lower_bound(instance_, customer_ids);
  const std::optional<double> finite_lower_bound =
      std::isfinite(lower_bound.objective_h)
          ? std::optional<double>(lower_bound.objective_h)
          : std::nullopt;
  const auto exact_result = [&](Plan plan) {
    const double objective = plan.evaluation->raw_cost_h;
    return FrvcpResult{.status = COPT_STATUS_OPTIMAL,
                       .status_name = "OPTIMAL",
                       .plan = std::move(plan),
                       .model_objective_h = objective,
                       .runtime_s = 0.0,
                       .wall_runtime_s = wall_seconds(call_started),
                       .cache_hit = true,
                       .applied_time_limit_s = std::nullopt,
                       .objective_lower_bound_h = finite_lower_bound,
                       .mip_start_objective_h = std::nullopt,
                       .model_size = {},
                       .network_size = FrvcpNetworkSize{
                           .customer_count = customer_ids.size(),
                           .gap_count = customer_ids.size() + 1}};
  };
  if (auto exact = factory_.exact_plan(customer_ids); exact.has_value()) {
    ++cache_hits_;
    return exact_result(std::move(*exact));
  }
  // This label now requires an independent fixed-customer-order certificate,
  // rather than OPTIMAL on the incumbent-dependent TRC-pruned path set.
  if (incumbent_plan != nullptr && incumbent_plan->exact_charging) {
    factory_.publish_exact(*incumbent_plan);
    ++cache_hits_;
    return exact_result(*incumbent_plan);
  }
  const ModelRequest request{
      .customers = customer_ids,
      .incumbent_visits = incumbent_plan != nullptr
          ? std::optional<std::vector<int>>(incumbent_plan->route.visits)
          : std::nullopt,
      .objective_upper_bound_bits = std::bit_cast<std::uint64_t>(
          incumbent_plan != nullptr ? incumbent_plan->evaluation->raw_cost_h
              : instance_.vehicle().max_route_duration_h - lower_bound.service_time_h)};
  if (auto cached = model_cache_.get(request); cached.has_value()) {
    ++cache_hits_;
    cached->runtime_s = 0.0;
    cached->wall_runtime_s = wall_seconds(call_started);
    cached->path_generation_runtime_s = 0.0;
    cached->model_build_runtime_s = 0.0;
    cached->cache_hit = true;
    cached->applied_time_limit_s = std::nullopt;
    cached->mip_start_objective_h = std::nullopt;
    cached->mip_start_supplied = false;
    cached->model_size = {};
    cached->network_size = FrvcpNetworkSize{
        .customer_count = customer_ids.size(), .gap_count = customer_ids.size() + 1};
    return *cached;
  }

  // With Euclidean travel, each direct anchor leg is a shortest path.  If the
  // all-direct fixed-order route is feasible on the initial battery, every
  // station detour has no smaller travel time and adds nonnegative charging
  // time.  This is therefore an exact certificate, not a heuristic shortcut.
  Plan direct = factory_.make_plan(customer_ids,
                                   direct_gaps(customer_ids.size()));
  if (direct.evaluation->feasible) {
    direct.exact_charging = true;
    const double objective = direct.evaluation->raw_cost_h;
    factory_.publish_exact(direct);
    return FrvcpResult{.status = COPT_STATUS_OPTIMAL,
                       .status_name = "OPTIMAL",
                       .plan = std::move(direct),
                       .model_objective_h = objective,
                       .runtime_s = 0.0,
                       .wall_runtime_s = wall_seconds(call_started),
                       .cache_hit = false,
                       .applied_time_limit_s = std::nullopt,
                       .objective_lower_bound_h = finite_lower_bound,
                       .mip_start_objective_h = std::nullopt,
                       .mip_start_supplied = false,
                       .analytic_optimum = true,
                       .lower_bound_proved_infeasible = false,
                       .model_size = {},
                       .network_size = FrvcpNetworkSize{
                           .customer_count = customer_ids.size(),
                           .gap_count = customer_ids.size() + 1}};
  }

  if (!std::isfinite(lower_bound.objective_h) ||
      lower_bound.objective_h + lower_bound.service_time_h >
          instance_.vehicle().max_route_duration_h + kTimeToleranceH) {
    FrvcpResult certificate{.status = COPT_STATUS_INFEASIBLE,
                       .status_name = "INFEASIBLE",
                       .plan = std::nullopt,
                       .model_objective_h = std::nullopt,
                       .runtime_s = 0.0,
                       .wall_runtime_s = wall_seconds(call_started),
                       .cache_hit = false,
                       .applied_time_limit_s = std::nullopt,
                       .objective_lower_bound_h = finite_lower_bound,
                       .mip_start_objective_h = std::nullopt,
                       .mip_start_supplied = false,
                       .analytic_optimum = false,
                       .lower_bound_proved_infeasible = true,
                       .model_size = {},
                       .network_size = FrvcpNetworkSize{
                           .customer_count = customer_ids.size(),
                           .gap_count = customer_ids.size() + 1}};
    model_cache_.put(request, certificate);
    return certificate;
  }
  ++solve_calls_;
  FrvcpResult result = solve(customer_ids, incumbent_plan,
                             lower_bound.objective_h, deadline, progress);
  if (result.status == COPT_STATUS_OPTIMAL || result.status == COPT_STATUS_INFEASIBLE) {
    model_cache_.put(request, result);
  }
  return result;
}

FrvcpResult CoptFrvcp::solve(const std::vector<int>& customer_ids,
                               const Plan* incumbent_plan,
                               double objective_lower_bound_h,
                               CoptDeadline deadline,
                               const FrvcpProgressCallback& progress) {
  const auto call_started = ModelClock::now();
  const CoptDeadline operation_deadline = deadline;
  std::vector<int> anchors;
  anchors.reserve(customer_ids.size() + 2);
  anchors.push_back(instance_.depot().id);
  anchors.insert(anchors.end(), customer_ids.begin(), customer_ids.end());
  anchors.push_back(instance_.depot().id);

  double service_time_h = 0.0;
  for (const int customer_id : customer_ids) {
    service_time_h += instance_.node(customer_id).service_time_h;
  }
  const double route_objective_upper_bound_h =
      incumbent_plan != nullptr
          ? incumbent_plan->evaluation->raw_cost_h
          : instance_.vehicle().max_route_duration_h - service_time_h;
  const auto path_generation_started = ModelClock::now();
  FrvcpPathCatalog catalog = enumerate_frvcp_paths(
      instance_, anchors,
      incumbent_plan != nullptr ? &incumbent_plan->gaps : nullptr,
      route_objective_upper_bound_h, operation_deadline);
  const double path_generation_runtime_s =
      wall_seconds(path_generation_started);

  const auto terminal_result = [&](int status) {
    const bool retain_incumbent =
        incumbent_plan != nullptr &&
        (status == COPT_STATUS_TIMEOUT || status == COPT_STATUS_INTERRUPTED ||
         status == COPT_STATUS_IMPRECISE || status == COPT_STATUS_NUMERICAL);
    return FrvcpResult{
        .status = status,
        .status_name = copt_status_name(status),
        .plan = retain_incumbent ? std::optional<Plan>(*incumbent_plan)
                                 : std::nullopt,
        .model_objective_h =
            retain_incumbent
                ? std::optional<double>(
                      incumbent_plan->evaluation->raw_cost_h)
                : std::nullopt,
        .runtime_s = 0.0,
        .wall_runtime_s = wall_seconds(call_started),
        .path_generation_runtime_s = path_generation_runtime_s,
        .model_build_runtime_s = 0.0,
        .cache_hit = false,
        .applied_time_limit_s = std::nullopt,
        .objective_lower_bound_h = objective_lower_bound_h,
        .mip_start_objective_h =
            retain_incumbent
                ? std::optional<double>(
                      incumbent_plan->evaluation->raw_cost_h)
                : std::nullopt,
        .mip_start_supplied = false,
        .analytic_optimum = false,
        .lower_bound_proved_infeasible = false,
        .model_size = {},
        .network_size = catalog.size};
  };
  if (!catalog.complete) return terminal_result(COPT_STATUS_TIMEOUT);
  if (catalog.gaps.size() != anchors.size() - 1 ||
      std::any_of(catalog.gaps.begin(), catalog.gaps.end(),
                  [](const FrvcpGapPathCatalog& gap) {
                    return gap.paths.empty();
                  })) {
    return terminal_result(COPT_STATUS_INFEASIBLE);
  }

  // Compress only variables for consecutive occurrences of the same physical
  // station. Single-valued F(E) and exactly zero self legs make charge times
  // telescope. No path, cache key, fallback Plan or external visit is removed.
  std::size_t compressed_runs = 0, merged_occurrences = 0;
  std::size_t skipped_non_strict_energy = 0, skipped_nonzero_self_leg = 0;
  catalog.size.active_station_copies = 0;
  for (auto& gap : catalog.gaps) {
    for (auto& spec : gap.paths) {
      for (std::size_t first = 0; first < spec.stations.size();) {
        std::size_t last = first;
        while (last + 1 < spec.stations.size() &&
               spec.stations[last + 1] == spec.stations[first]) ++last;
        bool merge = last > first;
        if (merge) {
          const int id = spec.stations[first];
          const auto& points = instance_.curve_for_station(id).points();
          const bool strict_energy = std::adjacent_find(
              points.begin(), points.end(), [](const ChargingPoint& a,
                                               const ChargingPoint& b) {
                return !(a.energy_wh < b.energy_wh);
              }) == points.end();
          const bool zero_self_leg = instance_.travel_time_h(id, id) == 0.0 &&
              std::all_of(spec.leg_energy_wh.begin() + first + 1,
                          spec.leg_energy_wh.begin() + last + 1,
                          [](double energy) { return energy == 0.0; });
          if (!strict_energy) ++skipped_non_strict_energy;
          if (!zero_self_leg) ++skipped_nonzero_self_leg;
          merge = strict_energy && zero_self_leg;
        }
        if (merge) {
          spec.model_station_spans.push_back({first, last});
          ++compressed_runs;
          saturating_add(merged_occurrences, last - first);
        } else {
          for (std::size_t index = first; index <= last; ++index) {
            spec.model_station_spans.push_back({index, index});
          }
        }
        first = last + 1;
      }
      saturating_add(catalog.size.active_station_copies,
                     spec.model_station_spans.size());
    }
  }
  // Rare model-boundary diagnostics only, never a route-evaluator hot path.
  if (compressed_runs || skipped_non_strict_energy || skipped_nonzero_self_leg) {
    std::clog << "frvcp_station_model_runs original_occurrences="
              << catalog.size.active_station_visits
              << " model_copies=" << catalog.size.active_station_copies
              << " compressed_runs=" << compressed_runs
              << " merged_occurrences=" << merged_occurrences
              << " skipped_non_strict_energy=" << skipped_non_strict_energy
              << " skipped_nonzero_self_leg=" << skipped_nonzero_self_leg << '\n';
  }

  struct PathStationVariables {
    int station_id{};
    std::size_t first_occurrence{};
    std::size_t last_occurrence{};
    Var arrival;
    Var departure;
    Var curve_arrival;
    Var curve_departure;
    Var charge_duration;
    PiecewiseCurve arrival_curve;
    PiecewiseCurve departure_curve;
  };
  struct PathVariables {
    const FrvcpPathSpec* spec{};
    Var selected;
    Var origin_energy;
    Var target_energy;
    std::vector<PathStationVariables> stations;
  };
  struct GapVariables {
    std::vector<PathVariables> paths;
  };
  struct CurvePoints {
    std::vector<double> energy;
    std::vector<double> time;
  };

  const auto model_build_started = ModelClock::now();
  Model model = environment_.CreateModel(("frvcp_path_" + instance_.name()).c_str());
  configure(model);
  const double capacity = instance_.vehicle().battery_capacity_wh;
  std::unordered_map<std::string, CurvePoints> curve_points;
  for (const auto& [type, curve] : instance_.charging_curves()) {
    CurvePoints points;
    points.energy.reserve(curve.points().size());
    points.time.reserve(curve.points().size());
    for (const ChargingPoint& point : curve.points()) {
      points.energy.push_back(point.energy_wh);
      points.time.push_back(point.time_h);
    }
    curve_points.emplace(type, std::move(points));
  }

  std::vector<Var> anchor_energy;
  anchor_energy.reserve(anchors.size());
  for (std::size_t index = 0; index < anchors.size(); ++index) {
    const double lower = index == 0 ? capacity : 0.0;
    anchor_energy.push_back(model.AddVar(lower, capacity, 0.0, COPT_CONTINUOUS, ("anchor_energy_" + std::to_string(index)).c_str()));
  }

  std::vector<GapVariables> gaps;
  gaps.reserve(catalog.gaps.size());
  Expr objective = 0.0;
  for (std::size_t gap_index = 0; gap_index < catalog.gaps.size();
       ++gap_index) {
    if (deadline_exhausted(operation_deadline)) {
      return terminal_result(COPT_STATUS_TIMEOUT);
    }
    const FrvcpGapPathCatalog& path_catalog = catalog.gaps[gap_index];
    GapVariables gap;
    gap.paths.reserve(path_catalog.paths.size());
    Expr selection = 0.0;
    Expr origin_flow = 0.0;
    Expr target_flow = 0.0;
    for (std::size_t path_index = 0;
         path_index < path_catalog.paths.size(); ++path_index) {
      const FrvcpPathSpec& spec = path_catalog.paths[path_index];
      const std::string suffix = std::to_string(gap_index) + "_" +
                                 std::to_string(path_index);
      Var selected = model.AddVar(0.0, 1.0, 0.0, COPT_BINARY, ("path_selected_" + suffix).c_str());
      Var origin_energy = model.AddVar(0.0, capacity, 0.0, COPT_CONTINUOUS, ("path_origin_energy_" + suffix).c_str());
      Var target_energy = model.AddVar(0.0, capacity, 0.0, COPT_CONTINUOUS, ("path_target_energy_" + suffix).c_str());
      model.AddConstr(origin_energy <= capacity * selected, ("path_origin_active_" + suffix).c_str());
      model.AddConstr(target_energy <= capacity * selected, ("path_target_active_" + suffix).c_str());
      objective += spec.travel_time_h * selected;
      selection += selected;
      origin_flow += origin_energy;
      target_flow += target_energy;

      PathVariables variables{.spec = &spec,
                              .selected = selected,
                              .origin_energy = origin_energy,
                              .target_energy = target_energy,
                              .stations = {}};
      variables.stations.reserve(spec.model_station_spans.size());
      Expr previous_departure = origin_energy;
      for (const FrvcpStationSpan span : spec.model_station_spans) {
        const std::size_t station_index = span.first;
        const int station_id = spec.stations[station_index];
        const std::string station_suffix =
            suffix + "_" + std::to_string(station_index);
        const ChargingCurve& curve = instance_.curve_for_station(station_id);
        const double curve_upper = curve.time_at_energy(capacity);
        Var arrival = model.AddVar(0.0, capacity, 0.0, COPT_CONTINUOUS, ("path_station_arrival_" + station_suffix).c_str());
        Var departure = model.AddVar(0.0, capacity, 0.0, COPT_CONTINUOUS, ("path_station_departure_" + station_suffix).c_str());
        Var curve_arrival = model.AddVar(0.0, curve_upper, 0.0, COPT_CONTINUOUS, ("path_curve_arrival_" + station_suffix).c_str());
        Var curve_departure = model.AddVar(0.0, curve_upper, 0.0, COPT_CONTINUOUS, ("path_curve_departure_" + station_suffix).c_str());
        Var charge_duration = model.AddVar(0.0, curve_upper, 0.0, COPT_CONTINUOUS, ("path_charge_duration_" + station_suffix).c_str());
        model.AddConstr(arrival <= capacity * selected, ("path_station_arrival_active_" + station_suffix).c_str());
        model.AddConstr(departure <= capacity * selected, ("path_station_departure_active_" + station_suffix).c_str());
        model.AddConstr(departure >= arrival, ("path_station_charge_" + station_suffix).c_str());
        model.AddConstr(arrival == previous_departure -
                                       spec.leg_energy_wh[station_index] *
                                           selected, ("path_station_energy_" + station_suffix).c_str());
        const CurvePoints& points =
            curve_points.at(instance_.node(station_id).station_type);
        PiecewiseCurve arrival_pwl = add_piecewise_curve(model,
            arrival, curve_arrival, selected, static_cast<int>(points.energy.size()),
            points.energy.data(), points.time.data(),
            "path_arrival_curve_" + station_suffix);
        PiecewiseCurve departure_pwl = add_piecewise_curve(model,
            departure, curve_departure, selected,
            static_cast<int>(points.energy.size()), points.energy.data(),
            points.time.data(), "path_departure_curve_" + station_suffix);
        model.AddConstr(charge_duration == curve_departure - curve_arrival, ("path_charge_time_" + station_suffix).c_str());
        objective += charge_duration;
        variables.stations.push_back(PathStationVariables{
            .station_id = station_id,
            .first_occurrence = span.first,
            .last_occurrence = span.last,
            .arrival = arrival,
            .departure = departure,
            .curve_arrival = curve_arrival,
            .curve_departure = curve_departure,
            .charge_duration = charge_duration,
            .arrival_curve = std::move(arrival_pwl),
            .departure_curve = std::move(departure_pwl)});
        previous_departure = departure;
      }
      model.AddConstr(target_energy ==
              previous_departure - spec.leg_energy_wh.back() * selected, ("path_target_energy_balance_" + suffix).c_str());
      gap.paths.push_back(std::move(variables));
    }
    model.AddConstr(selection == 1.0, ("select_path_" + std::to_string(gap_index)).c_str());
    model.AddConstr(origin_flow == anchor_energy[gap_index], ("path_origin_flow_" + std::to_string(gap_index)).c_str());
    model.AddConstr(target_flow == anchor_energy[gap_index + 1], ("path_target_flow_" + std::to_string(gap_index)).c_str());
    gaps.push_back(std::move(gap));
  }
  model.AddConstr(objective + service_time_h <=
                      instance_.vehicle().max_route_duration_h, "route_duration");
  model.SetObjective(objective, COPT_MINIMIZE);

  bool mip_start_supplied = false;
  std::optional<double> mip_start_objective_h;
  if (incumbent_plan != nullptr) {
    const RouteEvaluation& evaluation = *incumbent_plan->evaluation;
    const std::size_t expected_energy_count =
        incumbent_plan->route.visits.size() + 2;
    bool representable =
        evaluation.arrival_energy_wh.size() == expected_energy_count &&
        evaluation.departure_energy_wh.size() == expected_energy_count;
    std::vector<std::size_t> anchor_positions(anchors.size(), 0);
    std::vector<std::vector<std::size_t>> station_positions(gaps.size());
    std::vector<std::size_t> selected_path_indices(gaps.size(), 0);
    std::size_t visit_index = 0;
    std::size_t route_position = 0;
    for (std::size_t gap_index = 0;
         representable && gap_index < gaps.size(); ++gap_index) {
      const auto selected_path = std::find_if(
          gaps[gap_index].paths.begin(), gaps[gap_index].paths.end(),
          [&](const PathVariables& candidate) {
            return candidate.spec->stations ==
                   incumbent_plan->gaps[gap_index];
          });
      if (selected_path == gaps[gap_index].paths.end()) {
        representable = false;
        break;
      }
      selected_path_indices[gap_index] = static_cast<std::size_t>(
          std::distance(gaps[gap_index].paths.begin(), selected_path));
      for (const int station_id : incumbent_plan->gaps[gap_index]) {
        ++route_position;
        if (visit_index >= incumbent_plan->route.visits.size() ||
            incumbent_plan->route.visits[visit_index++] != station_id) {
          representable = false;
          break;
        }
        station_positions[gap_index].push_back(route_position);
      }
      if (!representable) break;
      ++route_position;
      if (gap_index < customer_ids.size() &&
          (visit_index >= incumbent_plan->route.visits.size() ||
           incumbent_plan->route.visits[visit_index++] !=
               customer_ids[gap_index])) {
        representable = false;
        break;
      }
      anchor_positions[gap_index + 1] = route_position;
    }
    representable = representable &&
                    visit_index == incumbent_plan->route.visits.size() &&
                    route_position + 1 == expected_energy_count;

    if (representable) {
      for (std::size_t index = 0; index < anchor_energy.size(); ++index) {
        model.SetMipStart(anchor_energy[index], std::clamp(evaluation.departure_energy_wh[anchor_positions[index]],
                       0.0, capacity));
      }
      for (std::size_t gap_index = 0; gap_index < gaps.size(); ++gap_index) {
        for (PathVariables& path : gaps[gap_index].paths) {
          model.SetMipStart(path.selected, 0.0);
          model.SetMipStart(path.origin_energy, 0.0);
          model.SetMipStart(path.target_energy, 0.0);
          for (PathStationVariables& station : path.stations) {
            model.SetMipStart(station.arrival, 0.0);
            model.SetMipStart(station.departure, 0.0);
            model.SetMipStart(station.curve_arrival, 0.0);
            model.SetMipStart(station.curve_departure, 0.0);
            model.SetMipStart(station.charge_duration, 0.0);
            station.arrival_curve.set_start(model, 0.0, false);
            station.departure_curve.set_start(model, 0.0, false);
          }
        }
        PathVariables& selected =
            gaps[gap_index].paths[selected_path_indices[gap_index]];
        model.SetMipStart(selected.selected, 1.0);
        model.SetMipStart(selected.origin_energy, evaluation.departure_energy_wh[anchor_positions[gap_index]]);
        model.SetMipStart(selected.target_energy, evaluation.departure_energy_wh[anchor_positions[gap_index + 1]]);
        for (std::size_t station_index = 0;
             station_index < selected.stations.size(); ++station_index) {
          PathStationVariables& station = selected.stations[station_index];
          const std::size_t first_position =
              station_positions[gap_index][station.first_occurrence];
          const std::size_t last_position =
              station_positions[gap_index][station.last_occurrence];
          const double arrival = evaluation.arrival_energy_wh[first_position];
          const double departure = evaluation.departure_energy_wh[last_position];
          const ChargingCurve& curve =
              instance_.curve_for_station(station.station_id);
          model.SetMipStart(station.arrival, arrival);
          model.SetMipStart(station.departure, departure);
          station.arrival_curve.set_start(model, arrival);
          station.departure_curve.set_start(model, departure);
          model.SetMipStart(station.curve_arrival, curve.time_at_energy(arrival));
          model.SetMipStart(station.curve_departure, curve.time_at_energy(departure));
          model.SetMipStart(station.charge_duration, curve.charge_duration(arrival, departure));
        }
      }
      mip_start_supplied = true;
      mip_start_objective_h = incumbent_plan->evaluation->raw_cost_h;
    }
  }

  const CoptModelSize size = model_size(model);
  const double model_build_runtime_s = wall_seconds(model_build_started);
  const auto emit = [&](FrvcpProgressEvent event) {
    if (!progress) return;
    progress(FrvcpProgress{
        .event = event,
        .wall_runtime_s = wall_seconds(call_started),
        .path_generation_runtime_s = path_generation_runtime_s,
        .model_build_runtime_s = model_build_runtime_s,
        .mip_start_supplied = mip_start_supplied,
        .model_size = size,
        .network_size = catalog.size});
  };
  emit(FrvcpProgressEvent::ModelBuilt);
  ModelTimeBudget budget = model_time_budget(config_, operation_deadline);
  if (budget.expired) {
    FrvcpResult result = terminal_result(COPT_STATUS_TIMEOUT);
    result.model_build_runtime_s = model_build_runtime_s;
    result.model_size = size;
    result.mip_start_supplied = mip_start_supplied;
    result.mip_start_objective_h = mip_start_objective_h;
    return result;
  }
  apply_time_budget(model, budget);
  if (mip_start_supplied) {
    model.LoadMipStart();
    emit(FrvcpProgressEvent::MipStartLoaded);
  }
  emit(FrvcpProgressEvent::SolveStarted);
  // Loading the start and observation callbacks also consume the wall budget.
  budget = model_time_budget(config_, operation_deadline);
  if (budget.expired) {
    FrvcpResult result = terminal_result(COPT_STATUS_TIMEOUT);
    result.model_build_runtime_s = model_build_runtime_s;
    result.model_size = size;
    result.mip_start_supplied = mip_start_supplied;
    result.mip_start_objective_h = mip_start_objective_h;
    return result;
  }
  apply_time_budget(model, budget);
  model.Solve();
  const int status = model.GetIntAttr(COPT_INTATTR_STATUS);
  const int solution_count = model.GetIntAttr(COPT_INTATTR_HASSOL);
  const double runtime = model.GetDblAttr(COPT_DBLATTR_SOLVINGTIME);
  const auto unusable_solution = [&](int result_status) {
    const bool retain_incumbent =
        incumbent_plan != nullptr &&
        (result_status == COPT_STATUS_TIMEOUT || result_status == COPT_STATUS_INTERRUPTED ||
         result_status == COPT_STATUS_IMPRECISE || result_status == COPT_STATUS_NUMERICAL);
    return FrvcpResult{
        .status = result_status,
        .status_name = copt_status_name(result_status),
        .plan = retain_incumbent ? std::optional<Plan>(*incumbent_plan)
                                 : std::nullopt,
        .model_objective_h =
            retain_incumbent
                ? std::optional<double>(
                      incumbent_plan->evaluation->raw_cost_h)
                : std::nullopt,
        .runtime_s = runtime,
        .wall_runtime_s = wall_seconds(call_started),
        .path_generation_runtime_s = path_generation_runtime_s,
        .model_build_runtime_s = model_build_runtime_s,
        .cache_hit = false,
        .applied_time_limit_s = budget.applied_time_limit_s,
        .objective_lower_bound_h = objective_lower_bound_h,
        .mip_start_objective_h = mip_start_objective_h,
        .mip_start_supplied = mip_start_supplied,
        .analytic_optimum = false,
        .lower_bound_proved_infeasible = false,
        .model_size = size,
        .network_size = catalog.size};
  };
  if (solution_count == 0) return unusable_solution(status);
  const double model_objective = model.GetDblAttr(COPT_DBLATTR_OBJVAL);
  if (!finite_copt_value(model_objective)) {
    return unusable_solution(COPT_STATUS_NUMERICAL);
  }

  std::vector<std::vector<int>> station_paths;
  std::vector<std::vector<double>> station_departures;
  station_paths.reserve(gaps.size());
  station_departures.reserve(gaps.size());
  for (GapVariables& gap : gaps) {
    PathVariables* selected = nullptr;
    for (PathVariables& candidate : gap.paths) {
      if (candidate.selected.Get(COPT_DBLINFO_VALUE) > 0.5) {
        if (selected != nullptr) {
          throw std::runtime_error("FRVCP selected multiple paths in one gap");
        }
        selected = &candidate;
      }
    }
    if (selected == nullptr) {
      throw std::runtime_error("FRVCP did not select a path in one gap");
    }
    station_paths.push_back(selected->spec->stations);
    std::vector<double> departures;
    departures.reserve(selected->spec->stations.size());
    for (PathStationVariables& station : selected->stations) {
      // Canonical lift: charge at the first original occurrence; every later
      // identical occurrence has zero charge. Preserve all external indices.
      const double departure = station.departure.Get(COPT_DBLINFO_VALUE);
      for (std::size_t occurrence = station.first_occurrence;
           occurrence <= station.last_occurrence; ++occurrence) {
        departures.push_back(departure);
      }
    }
    station_departures.push_back(std::move(departures));
  }

  Route route = expand_route(customer_ids, station_paths);
  if(factory_.has_adjacent_station_visits(route.visits)) {
    std::fprintf(stderr,"frvcp_plan_rejected reason=adjacent_identical_station retain_incumbent=%d\n",incumbent_plan!=nullptr);
    return unusable_solution(COPT_STATUS_NUMERICAL);
  }
  std::vector<std::optional<double>> schedule(route.visits.size());
  std::size_t visit_index = 0;
  for (std::size_t gap_index = 0; gap_index < station_paths.size();
       ++gap_index) {
    for (std::size_t station_index = 0;
         station_index < station_paths[gap_index].size(); ++station_index) {
      schedule[visit_index++] = station_departures[gap_index][station_index];
    }
    if (gap_index < customer_ids.size()) ++visit_index;
  }
  std::shared_ptr<const RouteEvaluation> evaluation;
  try {
    evaluation = evaluator_.evaluate_schedule(route, schedule);
  } catch (const std::invalid_argument&) {
    // Backend feasibility tolerances can admit a physically invalid schedule.
    return unusable_solution(COPT_STATUS_NUMERICAL);
  }
  if (!detail::valid_frvcp_replay(*evaluation, model_objective)) {
    return unusable_solution(COPT_STATUS_NUMERICAL);
  }
  // OPTIMAL certifies this pruned model only. A feasible objective meeting an
  // independent lower bound additionally certifies the full customer order.
  const bool fixed_order_certificate = status == COPT_STATUS_OPTIMAL &&
      std::isfinite(objective_lower_bound_h) &&
      evaluation->raw_cost_h <= objective_lower_bound_h + kCostTolerance;
  Plan plan{.customer_ids = customer_ids,
            .gaps = std::move(station_paths),
            .route = std::move(route),
            .evaluation = std::move(evaluation),
            .exact_charging = fixed_order_certificate};
  if (fixed_order_certificate) factory_.publish_exact(plan);
  return FrvcpResult{
      .status = status,
      .status_name = copt_status_name(status),
      .plan = std::move(plan),
      .model_objective_h = model_objective,
      .runtime_s = runtime,
      .wall_runtime_s = wall_seconds(call_started),
      .path_generation_runtime_s = path_generation_runtime_s,
      .model_build_runtime_s = model_build_runtime_s,
      .cache_hit = false,
      .applied_time_limit_s = budget.applied_time_limit_s,
      .objective_lower_bound_h = objective_lower_bound_h,
      .mip_start_objective_h = mip_start_objective_h,
      .mip_start_supplied = mip_start_supplied,
      .analytic_optimum = false,
      .lower_bound_proved_infeasible = false,
      .model_size = size,
      .network_size = catalog.size};
}


}  // namespace ils_sp
