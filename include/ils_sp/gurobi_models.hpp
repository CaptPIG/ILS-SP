#pragma once

#include <chrono>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "gurobi_c++.h"

#include "ils_sp/evaluator.hpp"
#include "ils_sp/pool.hpp"

namespace ils_sp {

using GurobiDeadline =
    std::optional<std::chrono::steady_clock::time_point>;

struct GurobiConfig {
  int threads{1};
  int seed{0};
  int output_flag{0};
  double mip_gap{1e-4};
  std::optional<double> time_limit_s;
};

struct GurobiModelSize {
  std::size_t variables{};
  std::size_t binary_variables{};
  std::size_t linear_constraints{};
  std::size_t general_constraints{};
  std::size_t linear_nonzeros{};
};

enum class SetPartitioningProgressEvent { ModelBuilt, Search };

struct SetPartitioningProgress {
  SetPartitioningProgressEvent event{SetPartitioningProgressEvent::ModelBuilt};
  double runtime_s{};
  double model_build_runtime_s{};
  std::optional<double> incumbent_objective_h;
  std::optional<double> best_bound_h;
  std::optional<double> relative_gap;
  double explored_nodes{};
  int solution_count{};
  bool rounded_target_reached{};
  bool incumbent_stall_limit_reached{};
  GurobiModelSize model_size;
};

using SetPartitioningProgressCallback =
    std::function<void(const SetPartitioningProgress&)>;

struct SetPartitioningSolveOptions {
  GurobiDeadline deadline;
  std::vector<std::size_t> mip_start_indices;
  std::optional<double> rounded_target_objective_h;
  std::optional<double> incumbent_stall_limit_s;
  std::optional<double> maximum_wall_time_s;
  SetPartitioningProgressCallback progress;
  double progress_interval_s{30.0};
};

struct FrvcpNetworkSize {
  std::size_t customer_count{};
  std::size_t gap_count{};
  std::size_t enumerated_paths{};
  std::size_t active_paths{};
  std::size_t enumerated_station_visits{};
  std::size_t active_station_visits{};
  std::size_t candidate_station_copies{};
  std::size_t active_station_copies{};
  std::size_t candidate_arcs{};
  std::size_t active_arcs{};
};

struct SetPartitioningResult {
  int status{};
  std::string status_name;
  std::vector<RouteColumn> selected_columns;
  std::optional<double> objective_h;
  double runtime_s{};
  double wall_runtime_s{};
  double model_build_runtime_s{};
  std::optional<double> applied_time_limit_s;
  std::optional<double> best_bound_h;
  std::optional<double> relative_gap;
  double explored_nodes{};
  int solution_count{};
  std::size_t mip_start_routes{};
  std::optional<double> mip_start_objective_h;
  bool rounded_target_reached{};
  bool incumbent_stall_limit_reached{};
  GurobiModelSize model_size;
};

struct LinearRelaxationResult {
  int status{};
  std::string status_name;
  std::optional<double> objective_h;
  std::unordered_map<int, double> duals;
  double runtime_s{};
  double wall_runtime_s{};
  std::optional<double> applied_time_limit_s;
  GurobiModelSize model_size;
};

class GurobiSetPartitioning {
 public:
  GurobiSetPartitioning(const Instance& instance, GRBEnv& environment,
                        GurobiConfig config = {});

  [[nodiscard]] SetPartitioningResult solve_integer(
      const std::vector<RouteColumn>& columns,
      SetPartitioningSolveOptions options = {});
  [[nodiscard]] LinearRelaxationResult solve_relaxation(
      const std::vector<RouteColumn>& columns,
      GurobiDeadline deadline = std::nullopt);

 private:
  [[nodiscard]] std::vector<std::vector<std::size_t>> customer_incidence(
      const std::vector<RouteColumn>& columns) const;
  void configure(GRBModel& model, bool integer) const;

  const Instance& instance_;
  GRBEnv& environment_;
  GurobiConfig config_;
};

struct FrvcpResult {
  int status{};
  std::string status_name;
  std::optional<Plan> plan;
  std::optional<double> model_objective_h;
  double runtime_s{};
  double wall_runtime_s{};
  double path_generation_runtime_s{};
  double model_build_runtime_s{};
  bool cache_hit{};
  std::optional<double> applied_time_limit_s;
  std::optional<double> objective_lower_bound_h;
  std::optional<double> mip_start_objective_h;
  bool mip_start_supplied{};
  bool analytic_optimum{};
  bool lower_bound_proved_infeasible{};
  GurobiModelSize model_size;
  FrvcpNetworkSize network_size;
};

class GurobiFrvcp {
 public:
  GurobiFrvcp(const Instance& instance, RouteEvaluator& evaluator,
              PlanFactory& factory, GRBEnv& environment,
              GurobiConfig config = {},
              std::size_t infeasible_cache_capacity = 20'000);

  [[nodiscard]] FrvcpResult optimize(
      const std::vector<int>& customer_ids,
      GurobiDeadline deadline = std::nullopt);
  [[nodiscard]] FrvcpResult optimize(
      const Plan& incumbent_plan,
      GurobiDeadline deadline = std::nullopt);
  [[nodiscard]] std::size_t solve_calls() const noexcept { return solve_calls_; }
  [[nodiscard]] std::size_t cache_hits() const noexcept { return cache_hits_; }

 private:
  [[nodiscard]] FrvcpResult optimize_impl(
      const std::vector<int>& customer_ids, const Plan* incumbent_plan,
      GurobiDeadline deadline);
  [[nodiscard]] FrvcpResult solve(const std::vector<int>& customer_ids,
                                  const Plan* incumbent_plan,
                                  double objective_lower_bound_h,
                                  GurobiDeadline deadline);
  void configure(GRBModel& model) const;

  const Instance& instance_;
  RouteEvaluator& evaluator_;
  PlanFactory& factory_;
  GRBEnv& environment_;
  GurobiConfig config_;
  LruCache<std::vector<int>, bool, IntVectorHash> infeasible_cache_;
  std::size_t solve_calls_{0};
  std::size_t cache_hits_{0};
};

[[nodiscard]] std::string gurobi_status_name(int status);

}  // namespace ils_sp
