#include "ils_sp/gurobi_models.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <exception>
#include <limits>
#include <map>
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
  result.status = GRB_TIME_LIMIT;
  result.status_name = gurobi_status_name(GRB_TIME_LIMIT);
  result.wall_runtime_s = wall_seconds(start);
  return result;
}

struct ModelTimeBudget {
  bool expired{};
  std::optional<double> applied_time_limit_s;
};

[[nodiscard]] bool deadline_exhausted(GurobiDeadline deadline) {
  return deadline.has_value() &&
         std::chrono::duration<double>(*deadline - ModelClock::now()).count() <=
             kDeadlineGuardSeconds;
}

[[nodiscard]] GurobiDeadline bounded_deadline(
    ModelClock::time_point start, GurobiDeadline deadline,
    std::optional<double> maximum_wall_time_s) {
  if (!maximum_wall_time_s.has_value()) return deadline;
  const auto local_deadline =
      start + std::chrono::duration_cast<ModelClock::duration>(
                  std::chrono::duration<double>(*maximum_wall_time_s));
  if (!deadline.has_value() || local_deadline < *deadline) {
    return local_deadline;
  }
  return deadline;
}

[[nodiscard]] ModelTimeBudget model_time_budget(
    const GurobiConfig& config, GurobiDeadline deadline) {
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

void apply_time_budget(GRBModel& model, const ModelTimeBudget& budget) {
  if (budget.applied_time_limit_s.has_value()) {
    model.set(GRB_DoubleParam_TimeLimit, *budget.applied_time_limit_s);
  }
}

[[nodiscard]] GurobiModelSize model_size(GRBModel& model) {
  model.update();
  return GurobiModelSize{
      .variables = static_cast<std::size_t>(
          model.get(GRB_IntAttr_NumVars)),
      .binary_variables = static_cast<std::size_t>(
          model.get(GRB_IntAttr_NumBinVars)),
      .linear_constraints = static_cast<std::size_t>(
          model.get(GRB_IntAttr_NumConstrs)),
      .general_constraints = static_cast<std::size_t>(
          model.get(GRB_IntAttr_NumGenConstrs)),
      .linear_nonzeros = static_cast<std::size_t>(
          model.get(GRB_IntAttr_NumNZs))};
}

[[nodiscard]] bool finite_gurobi_value(double value) {
  return std::isfinite(value) && std::abs(value) < GRB_INFINITY;
}

[[nodiscard]] std::optional<double> finite_optional(double value) {
  return finite_gurobi_value(value) ? std::optional<double>(value)
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
         std::llround(objective_h * 100.0) <=
             std::llround(*rounded_target_objective_h * 100.0);
}

[[nodiscard]] double rounded_target_stop_value(double target_h) {
  const double first_failing_value =
      (static_cast<double>(std::llround(target_h * 100.0)) + 0.5) / 100.0;
  return std::nextafter(first_failing_value,
                        -std::numeric_limits<double>::infinity());
}

struct FrvcpRouteLowerBound {
  double objective_h{};
  double service_time_h{};
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
  std::sort(curve_types.begin(), curve_types.end(),
            [&](const std::string& left, const std::string& right) {
              const bool left_no_slower = charging_curve_no_slower(
                  instance.charging_curves().at(left),
                  instance.charging_curves().at(right));
              const bool right_no_slower = charging_curve_no_slower(
                  instance.charging_curves().at(right),
                  instance.charging_curves().at(left));
              if (left_no_slower != right_no_slower) return left_no_slower;
              return left < right;
            });
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

[[nodiscard]] FrvcpPathCatalog enumerate_frvcp_paths(
    const Instance& instance, const std::vector<int>& anchors,
    const std::vector<std::vector<int>>* incumbent_gaps,
    double route_objective_upper_bound_h, GurobiDeadline deadline) {
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
    saturating_add(result.size.active_paths, gap.paths.size());
    for (const FrvcpPathSpec& path_spec : gap.paths) {
      saturating_add(result.size.active_station_visits,
                     path_spec.stations.size());
    }
    result.gaps.push_back(std::move(gap));
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

  double maximum_rate_wh_per_h = 0.0;
  bool zero_time_charge_segment = false;
  for (const int station_id : instance.station_ids()) {
    const std::vector<ChargingPoint>& points =
        instance.curve_for_station(station_id).points();
    for (std::size_t index = 1; index < points.size(); ++index) {
      const double energy_delta =
          points[index].energy_wh - points[index - 1].energy_wh;
      if (energy_delta <= kEnergyToleranceWh) continue;
      const double time_delta =
          points[index].time_h - points[index - 1].time_h;
      if (time_delta <= 1e-15) {
        zero_time_charge_segment = true;
        break;
      }
      maximum_rate_wh_per_h =
          std::max(maximum_rate_wh_per_h, energy_delta / time_delta);
    }
    if (zero_time_charge_segment) break;
  }
  const double charging_lower_bound_h =
      zero_time_charge_segment
          ? 0.0
          : (maximum_rate_wh_per_h > 0.0
                 ? required_charge_wh / maximum_rate_wh_per_h
                 : std::numeric_limits<double>::infinity());
  return FrvcpRouteLowerBound{
      .objective_h = direct_travel_h + charging_lower_bound_h,
      .service_time_h = service_time_h};
}

class SetPartitioningCallback final : public GRBCallback {
 public:
  SetPartitioningCallback(SetPartitioningProgressCallback progress,
                          std::optional<double> rounded_target_objective_h,
                          std::optional<double> incumbent_stall_limit_s,
                          std::optional<double> initial_incumbent_h,
                          double progress_interval_s,
                          GurobiModelSize size,
                          double model_build_runtime_s)
      : progress_(std::move(progress)),
        rounded_target_objective_h_(rounded_target_objective_h),
        incumbent_stall_limit_s_(incumbent_stall_limit_s),
        best_incumbent_h_(initial_incumbent_h),
        last_improvement_runtime_s_(initial_incumbent_h.has_value()
                                        ? std::optional<double>(0.0)
                                        : std::nullopt),
        progress_interval_s_(progress_interval_s),
        next_progress_s_(progress_interval_s),
        size_(size),
        model_build_runtime_s_(model_build_runtime_s) {}

  void rethrow_callback_error() const {
    if (callback_error_ != nullptr) std::rethrow_exception(callback_error_);
  }
  [[nodiscard]] bool incumbent_stall_limit_reached() const noexcept {
    return incumbent_stall_limit_reached_;
  }

 private:
  void callback() override {
    if (callback_error_ != nullptr ||
        (where != GRB_CB_MIP && where != GRB_CB_MIPSOL)) {
      return;
    }
    try {
      const double runtime_s = getDoubleInfo(GRB_CB_RUNTIME);
      const bool mip_solution = where == GRB_CB_MIPSOL;
      const double incumbent_value = getDoubleInfo(
          mip_solution ? GRB_CB_MIPSOL_OBJBST : GRB_CB_MIP_OBJBST);
      const double bound_value = getDoubleInfo(
          mip_solution ? GRB_CB_MIPSOL_OBJBND : GRB_CB_MIP_OBJBND);
      const double candidate_value =
          mip_solution ? getDoubleInfo(GRB_CB_MIPSOL_OBJ) : incumbent_value;
      const std::optional<double> observed_incumbent =
          finite_optional(candidate_value);
      if (observed_incumbent.has_value() &&
          (!best_incumbent_h_.has_value() ||
           *observed_incumbent < *best_incumbent_h_ - kCostTolerance)) {
        best_incumbent_h_ = observed_incumbent;
        last_improvement_runtime_s_ = runtime_s;
      }
      const bool target_reached =
          mip_solution && finite_gurobi_value(candidate_value) &&
          rounded_target_reached(candidate_value,
                                 rounded_target_objective_h_);
      const bool stalled =
          incumbent_stall_limit_s_.has_value() &&
          last_improvement_runtime_s_.has_value() &&
          runtime_s - *last_improvement_runtime_s_ + 1e-9 >=
              *incumbent_stall_limit_s_;

      const std::optional<double> incumbent =
          best_incumbent_h_.has_value() ? best_incumbent_h_
                                        : finite_optional(incumbent_value);
      const std::optional<double> bound = finite_optional(bound_value);
      if (progress_ &&
          (stalled || runtime_s + 1e-9 >= next_progress_s_)) {
        progress_(SetPartitioningProgress{
            .event = SetPartitioningProgressEvent::Search,
            .runtime_s = runtime_s,
            .model_build_runtime_s = model_build_runtime_s_,
            .incumbent_objective_h = incumbent,
            .best_bound_h = bound,
            .relative_gap = relative_gap(incumbent, bound),
            .explored_nodes = getDoubleInfo(
                mip_solution ? GRB_CB_MIPSOL_NODCNT : GRB_CB_MIP_NODCNT),
            .solution_count = getIntInfo(
                mip_solution ? GRB_CB_MIPSOL_SOLCNT : GRB_CB_MIP_SOLCNT),
            .rounded_target_reached = target_reached,
            .incumbent_stall_limit_reached = stalled,
            .model_size = size_});
        next_progress_s_ =
            (std::floor(runtime_s / progress_interval_s_) + 1.0) *
            progress_interval_s_;
      }
      if (stalled) {
        incumbent_stall_limit_reached_ = true;
        abort();
      }
    } catch (...) {
      callback_error_ = std::current_exception();
      abort();
    }
  }

  SetPartitioningProgressCallback progress_;
  std::optional<double> rounded_target_objective_h_;
  std::optional<double> incumbent_stall_limit_s_;
  std::optional<double> best_incumbent_h_;
  std::optional<double> last_improvement_runtime_s_;
  double progress_interval_s_{};
  double next_progress_s_{};
  GurobiModelSize size_;
  double model_build_runtime_s_{};
  bool incumbent_stall_limit_reached_{};
  std::exception_ptr callback_error_;
};

}  // namespace

std::string gurobi_status_name(int status) {
  switch (status) {
    case GRB_LOADED:
      return "LOADED";
    case GRB_OPTIMAL:
      return "OPTIMAL";
    case GRB_INFEASIBLE:
      return "INFEASIBLE";
    case GRB_INF_OR_UNBD:
      return "INF_OR_UNBD";
    case GRB_UNBOUNDED:
      return "UNBOUNDED";
    case GRB_TIME_LIMIT:
      return "TIME_LIMIT";
    case GRB_INTERRUPTED:
      return "INTERRUPTED";
    case GRB_NUMERIC:
      return "NUMERIC";
    case GRB_SUBOPTIMAL:
      return "SUBOPTIMAL";
    case GRB_USER_OBJ_LIMIT:
      return "USER_OBJ_LIMIT";
    default:
      return "STATUS_" + std::to_string(status);
  }
}

GurobiSetPartitioning::GurobiSetPartitioning(const Instance& instance,
                                             GRBEnv& environment,
                                             GurobiConfig config)
    : instance_(instance), environment_(environment), config_(config) {
  if (config_.threads < 1 || !std::isfinite(config_.mip_gap) ||
      config_.mip_gap < 0.0 ||
      (config_.time_limit_s.has_value() &&
       (!std::isfinite(*config_.time_limit_s) ||
        *config_.time_limit_s <= 0.0))) {
    throw std::invalid_argument("invalid Gurobi configuration");
  }
}

void GurobiSetPartitioning::configure(GRBModel& model, bool integer) const {
  model.set(GRB_IntParam_OutputFlag, config_.output_flag);
  model.set(GRB_IntParam_Threads, config_.threads);
  model.set(GRB_IntParam_Seed, config_.seed);
  if (integer) model.set(GRB_DoubleParam_MIPGap, config_.mip_gap);
}

std::vector<std::vector<std::size_t>>
GurobiSetPartitioning::customer_incidence(
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

SetPartitioningResult GurobiSetPartitioning::solve_integer(
    const std::vector<RouteColumn>& columns,
    SetPartitioningSolveOptions options) {
  const auto call_started = ModelClock::now();
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
  const GurobiDeadline deadline = bounded_deadline(
      call_started, options.deadline, options.maximum_wall_time_s);
  const std::vector<std::vector<std::size_t>> incidence =
      customer_incidence(columns);
  if (model_time_budget(config_, deadline).expired) {
    return time_limit_result<SetPartitioningResult>(call_started);
  }
  GRBModel model(environment_);
  model.set(GRB_StringAttr_ModelName, "ils_sp_partition_" + instance_.name());
  configure(model, true);
  std::vector<GRBVar> theta;
  theta.reserve(columns.size());
  for (std::size_t index = 0; index < columns.size(); ++index) {
    if ((index & 255U) == 0U && deadline_exhausted(deadline)) {
      return time_limit_result<SetPartitioningResult>(call_started);
    }
    theta.push_back(model.addVar(0.0, 1.0, columns[index].cost_h, GRB_BINARY,
                                 "theta_" + std::to_string(index)));
  }
  if (deadline_exhausted(deadline)) {
    return time_limit_result<SetPartitioningResult>(call_started);
  }
  const std::vector<int>& customers = instance_.customer_ids();
  for (std::size_t row = 0; row < customers.size(); ++row) {
    if (deadline_exhausted(deadline)) {
      return time_limit_result<SetPartitioningResult>(call_started);
    }
    GRBLinExpr cover = 0.0;
    for (const std::size_t column_index : incidence[row]) {
      cover += theta[column_index];
    }
    model.addConstr(cover == 1.0,
                    "cover_" + std::to_string(customers[row]));
  }

  std::optional<double> mip_start_objective_h;
  if (!options.mip_start_indices.empty()) {
    std::vector<bool> selected(columns.size(), false);
    for (const std::size_t index : options.mip_start_indices) {
      if (index >= columns.size() || selected[index]) {
        throw std::invalid_argument("invalid duplicate SP MIP-start column");
      }
      selected[index] = true;
    }
    for (const std::vector<std::size_t>& covering_columns : incidence) {
      const std::size_t count = static_cast<std::size_t>(std::count_if(
          covering_columns.begin(), covering_columns.end(),
          [&](std::size_t index) { return selected[index]; }));
      if (count != 1) {
        throw std::invalid_argument(
            "SP MIP start is not an exact customer partition");
      }
    }
    double objective_h = 0.0;
    for (std::size_t index = 0; index < columns.size(); ++index) {
      theta[index].set(GRB_DoubleAttr_Start, selected[index] ? 1.0 : 0.0);
      if (selected[index]) objective_h += columns[index].cost_h;
    }
    mip_start_objective_h = objective_h;
  }

  model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
  const GurobiModelSize size = model_size(model);
  const double model_build_runtime_s = wall_seconds(call_started);
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
  const ModelTimeBudget budget = model_time_budget(config_, deadline);
  if (budget.expired) {
    SetPartitioningResult result =
        time_limit_result<SetPartitioningResult>(call_started);
    result.model_build_runtime_s = model_build_runtime_s;
    result.mip_start_routes = options.mip_start_indices.size();
    result.mip_start_objective_h = mip_start_objective_h;
    result.model_size = size;
    return result;
  }
  apply_time_budget(model, budget);
  if (options.rounded_target_objective_h.has_value()) {
    model.set(GRB_DoubleParam_BestObjStop,
              rounded_target_stop_value(
                  *options.rounded_target_objective_h));
  }
  SetPartitioningCallback callback(
      options.progress, options.rounded_target_objective_h,
      options.incumbent_stall_limit_s, mip_start_objective_h,
      options.progress_interval_s, size, model_build_runtime_s);
  if (options.progress || options.incumbent_stall_limit_s.has_value()) {
    model.setCallback(&callback);
  }
  model.optimize();
  callback.rethrow_callback_error();
  const int status = model.get(GRB_IntAttr_Status);
  const int solution_count = model.get(GRB_IntAttr_SolCount);
  const std::optional<double> objective_h =
      solution_count > 0
          ? std::optional<double>(model.get(GRB_DoubleAttr_ObjVal))
          : std::nullopt;
  const std::optional<double> best_bound_h =
      finite_optional(model.get(GRB_DoubleAttr_ObjBound));
  std::vector<RouteColumn> selected;
  if (solution_count > 0) {
    for (std::size_t index = 0; index < columns.size(); ++index) {
      if (theta[index].get(GRB_DoubleAttr_X) > 0.5) {
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
      .status_name = gurobi_status_name(status),
      .selected_columns = std::move(selected),
      .objective_h = objective_h,
      .runtime_s = model.get(GRB_DoubleAttr_Runtime),
      .wall_runtime_s = wall_seconds(call_started),
      .model_build_runtime_s = model_build_runtime_s,
      .applied_time_limit_s = budget.applied_time_limit_s,
      .best_bound_h = best_bound_h,
      .relative_gap = relative_gap(objective_h, best_bound_h),
      .explored_nodes = model.get(GRB_DoubleAttr_NodeCount),
      .solution_count = solution_count,
      .mip_start_routes = options.mip_start_indices.size(),
      .mip_start_objective_h = mip_start_objective_h,
      .rounded_target_reached =
          objective_h.has_value() &&
          rounded_target_reached(*objective_h,
                                 options.rounded_target_objective_h),
      .incumbent_stall_limit_reached =
          callback.incumbent_stall_limit_reached(),
      .model_size = size};
}

LinearRelaxationResult GurobiSetPartitioning::solve_relaxation(
    const std::vector<RouteColumn>& columns, GurobiDeadline deadline) {
  const auto call_started = ModelClock::now();
  const std::vector<std::vector<std::size_t>> incidence =
      customer_incidence(columns);
  if (model_time_budget(config_, deadline).expired) {
    return time_limit_result<LinearRelaxationResult>(call_started);
  }
  GRBModel model(environment_);
  model.set(GRB_StringAttr_ModelName,
            "ils_sp_partition_lp_" + instance_.name());
  configure(model, false);
  std::vector<GRBVar> theta;
  theta.reserve(columns.size());
  for (std::size_t index = 0; index < columns.size(); ++index) {
    if ((index & 255U) == 0U && deadline_exhausted(deadline)) {
      return time_limit_result<LinearRelaxationResult>(call_started);
    }
    theta.push_back(model.addVar(0.0, 1.0, columns[index].cost_h,
                                 GRB_CONTINUOUS,
                                 "theta_" + std::to_string(index)));
  }
  if (deadline_exhausted(deadline)) {
    return time_limit_result<LinearRelaxationResult>(call_started);
  }
  std::unordered_map<int, GRBConstr> cover_constraints;
  const std::vector<int>& customers = instance_.customer_ids();
  for (std::size_t row = 0; row < customers.size(); ++row) {
    if (deadline_exhausted(deadline)) {
      return time_limit_result<LinearRelaxationResult>(call_started);
    }
    GRBLinExpr cover = 0.0;
    for (const std::size_t column_index : incidence[row]) {
      cover += theta[column_index];
    }
    cover_constraints.emplace(
        customers[row],
        model.addConstr(cover == 1.0,
                        "cover_" + std::to_string(customers[row])));
  }
  model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
  const GurobiModelSize size = model_size(model);
  const ModelTimeBudget budget = model_time_budget(config_, deadline);
  if (budget.expired) {
    LinearRelaxationResult result =
        time_limit_result<LinearRelaxationResult>(call_started);
    result.model_size = size;
    return result;
  }
  apply_time_budget(model, budget);
  model.optimize();
  const int status = model.get(GRB_IntAttr_Status);
  std::unordered_map<int, double> duals;
  if (status == GRB_OPTIMAL) {
    for (auto& [customer_id, constraint] : cover_constraints) {
      duals.emplace(customer_id, constraint.get(GRB_DoubleAttr_Pi));
    }
  }
  return LinearRelaxationResult{
      .status = status,
      .status_name = gurobi_status_name(status),
      .objective_h = status == GRB_OPTIMAL
                         ? std::optional<double>(
                               model.get(GRB_DoubleAttr_ObjVal))
                         : std::nullopt,
      .duals = std::move(duals),
      .runtime_s = model.get(GRB_DoubleAttr_Runtime),
      .wall_runtime_s = wall_seconds(call_started),
      .applied_time_limit_s = budget.applied_time_limit_s,
      .model_size = size};
}

GurobiFrvcp::GurobiFrvcp(const Instance& instance, RouteEvaluator& evaluator,
                         PlanFactory& factory, GRBEnv& environment,
                         GurobiConfig config,
                         std::size_t infeasible_cache_capacity)
    : instance_(instance),
      evaluator_(evaluator),
      factory_(factory),
      environment_(environment),
      config_(config),
      infeasible_cache_(infeasible_cache_capacity) {
  if (config_.threads < 1 || !std::isfinite(config_.mip_gap) ||
      config_.mip_gap < 0.0 ||
      (config_.time_limit_s.has_value() &&
       (!std::isfinite(*config_.time_limit_s) ||
        *config_.time_limit_s <= 0.0))) {
    throw std::invalid_argument("invalid FRVCP Gurobi configuration");
  }
}

void GurobiFrvcp::configure(GRBModel& model) const {
  model.set(GRB_IntParam_OutputFlag, config_.output_flag);
  model.set(GRB_IntParam_Threads, config_.threads);
  model.set(GRB_IntParam_Seed, config_.seed);
  model.set(GRB_DoubleParam_MIPGap, config_.mip_gap);
}

FrvcpResult GurobiFrvcp::optimize(
    const std::vector<int>& customer_ids, GurobiDeadline deadline) {
  return optimize_impl(customer_ids, nullptr, deadline);
}

FrvcpResult GurobiFrvcp::optimize(const Plan& incumbent_plan,
                                  GurobiDeadline deadline) {
  return optimize_impl(incumbent_plan.customer_ids, &incumbent_plan,
                       deadline);
}

FrvcpResult GurobiFrvcp::optimize_impl(
    const std::vector<int>& customer_ids, const Plan* incumbent_plan,
    GurobiDeadline deadline) {
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
       incumbent_plan->evaluation == nullptr ||
       !incumbent_plan->evaluation->feasible)) {
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
    return FrvcpResult{.status = GRB_OPTIMAL,
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
  // TRC marks an optimized route so that it is never optimized again.  The
  // label remains a certificate even if the bounded ordered-route LRU evicts
  // its entry; restore that entry without rebuilding a model.
  if (incumbent_plan != nullptr && incumbent_plan->exact_charging) {
    factory_.publish_exact(*incumbent_plan);
    ++cache_hits_;
    return exact_result(*incumbent_plan);
  }
  if (auto infeasible = infeasible_cache_.get(customer_ids);
      infeasible.has_value()) {
    ++cache_hits_;
    return FrvcpResult{.status = GRB_INFEASIBLE,
                       .status_name = "INFEASIBLE",
                       .plan = std::nullopt,
                       .model_objective_h = std::nullopt,
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
    return FrvcpResult{.status = GRB_OPTIMAL,
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
    infeasible_cache_.put(customer_ids, true);
    return FrvcpResult{.status = GRB_INFEASIBLE,
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
  }
  ++solve_calls_;
  FrvcpResult result = solve(customer_ids, incumbent_plan,
                             lower_bound.objective_h, deadline);
  if (result.status == GRB_INFEASIBLE) {
    infeasible_cache_.put(customer_ids, true);
  }
  return result;
}

FrvcpResult GurobiFrvcp::solve(const std::vector<int>& customer_ids,
                               const Plan* incumbent_plan,
                               double objective_lower_bound_h,
                               GurobiDeadline deadline) {
  const auto call_started = ModelClock::now();
  const GurobiDeadline operation_deadline = deadline;
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
        (status == GRB_TIME_LIMIT || status == GRB_INTERRUPTED ||
         status == GRB_SUBOPTIMAL);
    return FrvcpResult{
        .status = status,
        .status_name = gurobi_status_name(status),
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
  if (!catalog.complete) return terminal_result(GRB_TIME_LIMIT);
  if (catalog.gaps.size() != anchors.size() - 1 ||
      std::any_of(catalog.gaps.begin(), catalog.gaps.end(),
                  [](const FrvcpGapPathCatalog& gap) {
                    return gap.paths.empty();
                  })) {
    return terminal_result(GRB_INFEASIBLE);
  }

  struct PathStationVariables {
    int station_id{};
    GRBVar arrival;
    GRBVar departure;
    GRBVar curve_arrival;
    GRBVar curve_departure;
    GRBVar charge_duration;
  };
  struct PathVariables {
    const FrvcpPathSpec* spec{};
    GRBVar selected;
    GRBVar origin_energy;
    GRBVar target_energy;
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
  GRBModel model(environment_);
  model.set(GRB_StringAttr_ModelName, "frvcp_path_" + instance_.name());
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

  std::vector<GRBVar> anchor_energy;
  anchor_energy.reserve(anchors.size());
  for (std::size_t index = 0; index < anchors.size(); ++index) {
    const double lower = index == 0 ? capacity : 0.0;
    anchor_energy.push_back(model.addVar(
        lower, capacity, 0.0, GRB_CONTINUOUS,
        "anchor_energy_" + std::to_string(index)));
  }

  std::vector<GapVariables> gaps;
  gaps.reserve(catalog.gaps.size());
  GRBLinExpr objective = 0.0;
  for (std::size_t gap_index = 0; gap_index < catalog.gaps.size();
       ++gap_index) {
    if (deadline_exhausted(operation_deadline)) {
      return terminal_result(GRB_TIME_LIMIT);
    }
    const FrvcpGapPathCatalog& path_catalog = catalog.gaps[gap_index];
    GapVariables gap;
    gap.paths.reserve(path_catalog.paths.size());
    GRBLinExpr selection = 0.0;
    GRBLinExpr origin_flow = 0.0;
    GRBLinExpr target_flow = 0.0;
    for (std::size_t path_index = 0;
         path_index < path_catalog.paths.size(); ++path_index) {
      const FrvcpPathSpec& spec = path_catalog.paths[path_index];
      const std::string suffix = std::to_string(gap_index) + "_" +
                                 std::to_string(path_index);
      GRBVar selected = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                                     "path_selected_" + suffix);
      GRBVar origin_energy = model.addVar(
          0.0, capacity, 0.0, GRB_CONTINUOUS,
          "path_origin_energy_" + suffix);
      GRBVar target_energy = model.addVar(
          0.0, capacity, 0.0, GRB_CONTINUOUS,
          "path_target_energy_" + suffix);
      model.addConstr(origin_energy <= capacity * selected,
                      "path_origin_active_" + suffix);
      model.addConstr(target_energy <= capacity * selected,
                      "path_target_active_" + suffix);
      objective += spec.travel_time_h * selected;
      selection += selected;
      origin_flow += origin_energy;
      target_flow += target_energy;

      PathVariables variables{.spec = &spec,
                              .selected = selected,
                              .origin_energy = origin_energy,
                              .target_energy = target_energy,
                              .stations = {}};
      variables.stations.reserve(spec.stations.size());
      GRBLinExpr previous_departure = origin_energy;
      for (std::size_t station_index = 0;
           station_index < spec.stations.size(); ++station_index) {
        const int station_id = spec.stations[station_index];
        const std::string station_suffix =
            suffix + "_" + std::to_string(station_index);
        const ChargingCurve& curve = instance_.curve_for_station(station_id);
        const double curve_upper = curve.time_at_energy(capacity);
        GRBVar arrival = model.addVar(
            0.0, capacity, 0.0, GRB_CONTINUOUS,
            "path_station_arrival_" + station_suffix);
        GRBVar departure = model.addVar(
            0.0, capacity, 0.0, GRB_CONTINUOUS,
            "path_station_departure_" + station_suffix);
        GRBVar curve_arrival = model.addVar(
            0.0, curve_upper, 0.0, GRB_CONTINUOUS,
            "path_curve_arrival_" + station_suffix);
        GRBVar curve_departure = model.addVar(
            0.0, curve_upper, 0.0, GRB_CONTINUOUS,
            "path_curve_departure_" + station_suffix);
        GRBVar charge_duration = model.addVar(
            0.0, curve_upper, 0.0, GRB_CONTINUOUS,
            "path_charge_duration_" + station_suffix);
        model.addConstr(arrival <= capacity * selected,
                        "path_station_arrival_active_" + station_suffix);
        model.addConstr(departure <= capacity * selected,
                        "path_station_departure_active_" + station_suffix);
        model.addConstr(departure >= arrival,
                        "path_station_charge_" + station_suffix);
        model.addConstr(arrival == previous_departure -
                                       spec.leg_energy_wh[station_index] *
                                           selected,
                        "path_station_energy_" + station_suffix);
        const CurvePoints& points =
            curve_points.at(instance_.node(station_id).station_type);
        model.addGenConstrPWL(
            arrival, curve_arrival, static_cast<int>(points.energy.size()),
            points.energy.data(), points.time.data(),
            "path_arrival_curve_" + station_suffix);
        model.addGenConstrPWL(
            departure, curve_departure,
            static_cast<int>(points.energy.size()), points.energy.data(),
            points.time.data(), "path_departure_curve_" + station_suffix);
        model.addConstr(charge_duration == curve_departure - curve_arrival,
                        "path_charge_time_" + station_suffix);
        objective += charge_duration;
        variables.stations.push_back(PathStationVariables{
            .station_id = station_id,
            .arrival = arrival,
            .departure = departure,
            .curve_arrival = curve_arrival,
            .curve_departure = curve_departure,
            .charge_duration = charge_duration});
        previous_departure = departure;
      }
      model.addConstr(
          target_energy ==
              previous_departure - spec.leg_energy_wh.back() * selected,
          "path_target_energy_balance_" + suffix);
      gap.paths.push_back(std::move(variables));
    }
    model.addConstr(selection == 1.0,
                    "select_path_" + std::to_string(gap_index));
    model.addConstr(origin_flow == anchor_energy[gap_index],
                    "path_origin_flow_" + std::to_string(gap_index));
    model.addConstr(target_flow == anchor_energy[gap_index + 1],
                    "path_target_flow_" + std::to_string(gap_index));
    gaps.push_back(std::move(gap));
  }
  model.addConstr(objective + service_time_h <=
                      instance_.vehicle().max_route_duration_h,
                  "route_duration");
  model.setObjective(objective, GRB_MINIMIZE);

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
        anchor_energy[index].set(
            GRB_DoubleAttr_Start,
            std::clamp(evaluation.departure_energy_wh[anchor_positions[index]],
                       0.0, capacity));
      }
      for (std::size_t gap_index = 0; gap_index < gaps.size(); ++gap_index) {
        for (PathVariables& path : gaps[gap_index].paths) {
          path.selected.set(GRB_DoubleAttr_Start, 0.0);
          path.origin_energy.set(GRB_DoubleAttr_Start, 0.0);
          path.target_energy.set(GRB_DoubleAttr_Start, 0.0);
          for (PathStationVariables& station : path.stations) {
            station.arrival.set(GRB_DoubleAttr_Start, 0.0);
            station.departure.set(GRB_DoubleAttr_Start, 0.0);
            station.curve_arrival.set(GRB_DoubleAttr_Start, 0.0);
            station.curve_departure.set(GRB_DoubleAttr_Start, 0.0);
            station.charge_duration.set(GRB_DoubleAttr_Start, 0.0);
          }
        }
        PathVariables& selected =
            gaps[gap_index].paths[selected_path_indices[gap_index]];
        selected.selected.set(GRB_DoubleAttr_Start, 1.0);
        selected.origin_energy.set(
            GRB_DoubleAttr_Start,
            evaluation.departure_energy_wh[anchor_positions[gap_index]]);
        selected.target_energy.set(
            GRB_DoubleAttr_Start,
            evaluation.departure_energy_wh[anchor_positions[gap_index + 1]]);
        for (std::size_t station_index = 0;
             station_index < selected.stations.size(); ++station_index) {
          PathStationVariables& station = selected.stations[station_index];
          const std::size_t position =
              station_positions[gap_index][station_index];
          const double arrival = evaluation.arrival_energy_wh[position];
          const double departure = evaluation.departure_energy_wh[position];
          const ChargingCurve& curve =
              instance_.curve_for_station(station.station_id);
          station.arrival.set(GRB_DoubleAttr_Start, arrival);
          station.departure.set(GRB_DoubleAttr_Start, departure);
          station.curve_arrival.set(GRB_DoubleAttr_Start,
                                    curve.time_at_energy(arrival));
          station.curve_departure.set(GRB_DoubleAttr_Start,
                                      curve.time_at_energy(departure));
          station.charge_duration.set(
              GRB_DoubleAttr_Start,
              curve.charge_duration(arrival, departure));
        }
      }
      mip_start_supplied = true;
      mip_start_objective_h = incumbent_plan->evaluation->raw_cost_h;
    }
  }

  const GurobiModelSize size = model_size(model);
  const double model_build_runtime_s = wall_seconds(model_build_started);
  const ModelTimeBudget budget = model_time_budget(config_, operation_deadline);
  if (budget.expired) {
    FrvcpResult result = terminal_result(GRB_TIME_LIMIT);
    result.model_build_runtime_s = model_build_runtime_s;
    result.model_size = size;
    result.mip_start_supplied = mip_start_supplied;
    result.mip_start_objective_h = mip_start_objective_h;
    return result;
  }
  apply_time_budget(model, budget);
  model.optimize();
  const int status = model.get(GRB_IntAttr_Status);
  const int solution_count = model.get(GRB_IntAttr_SolCount);
  const double runtime = model.get(GRB_DoubleAttr_Runtime);
  if (solution_count == 0) {
    const bool retain_incumbent =
        incumbent_plan != nullptr &&
        (status == GRB_TIME_LIMIT || status == GRB_INTERRUPTED ||
         status == GRB_SUBOPTIMAL);
    return FrvcpResult{
        .status = status,
        .status_name = gurobi_status_name(status),
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
  }

  std::vector<std::vector<int>> station_paths;
  std::vector<std::vector<double>> station_departures;
  station_paths.reserve(gaps.size());
  station_departures.reserve(gaps.size());
  for (const GapVariables& gap : gaps) {
    const PathVariables* selected = nullptr;
    for (const PathVariables& candidate : gap.paths) {
      if (candidate.selected.get(GRB_DoubleAttr_X) > 0.5) {
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
    departures.reserve(selected->stations.size());
    for (const PathStationVariables& station : selected->stations) {
      departures.push_back(station.departure.get(GRB_DoubleAttr_X));
    }
    station_departures.push_back(std::move(departures));
  }

  Route route = expand_route(customer_ids, station_paths);
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
  auto evaluation = evaluator_.evaluate_schedule(route, schedule);
  const double model_objective = model.get(GRB_DoubleAttr_ObjVal);
  if (std::abs(evaluation->raw_cost_h - model_objective) > 1e-5) {
    throw std::runtime_error("FRVCP path model and replay objectives disagree");
  }
  Plan plan{.customer_ids = customer_ids,
            .gaps = std::move(station_paths),
            .route = std::move(route),
            .evaluation = std::move(evaluation),
            .exact_charging = status == GRB_OPTIMAL};
  if (status == GRB_OPTIMAL) factory_.publish_exact(plan);
  return FrvcpResult{
      .status = status,
      .status_name = gurobi_status_name(status),
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
