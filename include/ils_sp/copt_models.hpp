#pragma once

#include <chrono>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "coptcpp_pch.h"

#include "ils_sp/evaluator.hpp"
#include "ils_sp/pool.hpp"

namespace ils_sp {

using CoptDeadline =
    std::optional<std::chrono::steady_clock::time_point>;

struct CoptConfig {
  int threads{1};
  int output_flag{0};
  double mip_gap{1e-4};
  std::optional<double> time_limit_s;
  // Optional SP-integer override only; LP and FRVCP retain threads.
  std::optional<int> sp_integer_threads;
};

struct CoptModelSize {
  std::size_t variables{};
  std::size_t binary_variables{};
  std::size_t linear_constraints{};
  std::size_t sos_constraints{};
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
  std::optional<double> explored_nodes;
  int solution_count{};  // COPT exposes availability (0 or 1), not SolCount.
  bool rounded_target_reached{};
  bool incumbent_stall_limit_reached{};
  CoptModelSize model_size;
};

using SetPartitioningProgressCallback =
    std::function<void(const SetPartitioningProgress&)>;

struct SetPartitioningSolveOptions {
  CoptDeadline deadline;
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
  // May retain an independently validated start on timeout/numerical failure;
  // solution_count remains zero unless a usable backend solution was returned.
  std::vector<RouteColumn> selected_columns;
  std::optional<double> objective_h;
  double runtime_s{};
  double wall_runtime_s{};
  double model_build_runtime_s{};
  std::optional<double> applied_time_limit_s;
  std::optional<double> best_bound_h;
  std::optional<double> relative_gap;
  double explored_nodes{};
  int solution_count{};  // COPT exposes availability (0 or 1), not SolCount.
  std::size_t mip_start_routes{};
  std::optional<double> mip_start_objective_h;
  bool rounded_target_reached{};
  bool incumbent_stall_limit_reached{};
  CoptModelSize model_size;
};

struct LinearRelaxationResult {
  int status{};
  std::string status_name;
  std::optional<double> objective_h;
  std::unordered_map<int, double> duals;
  double runtime_s{};
  double wall_runtime_s{};
  std::optional<double> applied_time_limit_s;
  CoptModelSize model_size;
};

class CoptSetPartitioning {
 public:
  CoptSetPartitioning(const Instance& instance, Envr& environment,
                        CoptConfig config = {});

  [[nodiscard]] SetPartitioningResult solve_integer(
      const std::vector<RouteColumn>& columns,
      SetPartitioningSolveOptions options = {});
  [[nodiscard]] LinearRelaxationResult solve_relaxation(
      const std::vector<RouteColumn>& columns,
      CoptDeadline deadline = std::nullopt);

 private:
  [[nodiscard]] std::vector<std::vector<std::size_t>> customer_incidence(
      const std::vector<RouteColumn>& columns) const;
  void configure(Model& model, bool integer) const;

  const Instance& instance_;
  Envr& environment_;
  CoptConfig config_;
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
  CoptModelSize model_size;
  FrvcpNetworkSize network_size;
};

enum class FrvcpProgressEvent { ModelBuilt, MipStartLoaded, SolveStarted };

// Emitted before entering potentially expensive backend calls, including
// LoadMipStart. Observation only: no changes to the model or stopping policy.
struct FrvcpProgress {
  FrvcpProgressEvent event{FrvcpProgressEvent::ModelBuilt};
  double wall_runtime_s{};
  double path_generation_runtime_s{};
  double model_build_runtime_s{};
  bool mip_start_supplied{};
  CoptModelSize model_size;
  FrvcpNetworkSize network_size;
};
using FrvcpProgressCallback = std::function<void(const FrvcpProgress&)>;

class CoptFrvcp {
 public:
  CoptFrvcp(const Instance& instance, RouteEvaluator& evaluator,
              PlanFactory& factory, Envr& environment,
              CoptConfig config = {},
              std::size_t infeasible_cache_capacity = 20'000);

  [[nodiscard]] FrvcpResult optimize(
      const std::vector<int>& customer_ids,
      CoptDeadline deadline = std::nullopt,
      FrvcpProgressCallback progress = {});
  [[nodiscard]] FrvcpResult optimize(
      const Plan& incumbent_plan,
      CoptDeadline deadline = std::nullopt,
      FrvcpProgressCallback progress = {});
  [[nodiscard]] std::size_t solve_calls() const noexcept { return solve_calls_; }
  [[nodiscard]] std::size_t cache_hits() const noexcept { return cache_hits_; }

 private:
  struct ModelRequest {
    std::vector<int> customers;
    std::optional<std::vector<int>> incumbent_visits;
    std::uint64_t objective_upper_bound_bits{};
    bool operator==(const ModelRequest&) const = default;
  };
  struct ModelRequestHash {
    std::size_t operator()(const ModelRequest& key) const noexcept {
      std::size_t seed = IntVectorHash{}(key.customers);
      hash_combine(seed, key.incumbent_visits.has_value());
      if (key.incumbent_visits) {
        hash_combine(seed, IntVectorHash{}(*key.incumbent_visits));
      }
      hash_combine(seed, key.objective_upper_bound_bits);
      return seed;
    }
  };
  [[nodiscard]] FrvcpResult optimize_impl(
      const std::vector<int>& customer_ids, const Plan* incumbent_plan,
      CoptDeadline deadline, const FrvcpProgressCallback& progress);
  [[nodiscard]] FrvcpResult solve(const std::vector<int>& customer_ids,
                                  const Plan* incumbent_plan,
                                  double objective_lower_bound_h,
                                  CoptDeadline deadline,
                                  const FrvcpProgressCallback& progress);
  void configure(Model& model) const;

  const Instance& instance_;
  RouteEvaluator& evaluator_;
  PlanFactory& factory_;
  Envr& environment_;
  CoptConfig config_;
  // TRC pruning and retained incumbent paths define a restricted model.
  // Its optimal/infeasible certificate is reusable only for this request.
  // Keep the former infeasible-cache capacity as the total model-memo bound.
  LruCache<ModelRequest, FrvcpResult, ModelRequestHash> model_cache_;
  std::size_t solve_calls_{0};
  std::size_t cache_hits_{0};
};

[[nodiscard]] std::string copt_status_name(int status);

}  // namespace ils_sp
