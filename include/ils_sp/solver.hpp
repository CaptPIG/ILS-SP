#pragma once

#include <functional>

#include "ils_sp/gurobi_models.hpp"
#include "ils_sp/search.hpp"

namespace ils_sp {

struct PaperParameters {
  std::size_t penalty_window{50};
  double penalty_factor{1.5};
  double lambda_initial{100.0};
  double lambda_min{0.1};
  double lambda_max{10'000.0};
  std::size_t sp_frequency{1'500};
  std::size_t backup_frequency{1'500};
  std::size_t pool_capacity{5'000};
  std::size_t max_iterations{50'000};
};

struct SolverConfig {
  std::uint64_t seed{36'509};
  PaperParameters paper;
  std::size_t promising_arc_count{40};
  std::size_t granular_threshold{100};
  SearchProfile search_profile{SearchProfile::Paper};
  std::size_t evaluation_cache_capacity{200'000};
  std::size_t exact_cache_capacity{20'000};
  std::size_t infeasible_cache_capacity{20'000};
  GurobiConfig set_partitioning;
  GurobiConfig frvcp;
  std::optional<double> global_time_limit_s;
  std::optional<double> rounded_target_objective_h;
  std::optional<double> sp_incumbent_stall_limit_s{30.0};
  std::optional<double> sp_maximum_time_s{90.0};
  // Explicit non-paper extension. A full route pool may call SP before the
  // regular frequency after both record and coverage quality have stalled.
  std::optional<std::size_t> adaptive_sp_stall_iterations;
  std::size_t adaptive_sp_minimum_interval{500};
  std::size_t progress_interval{100};
};

enum class SolverCheckpointReason {
  Initialization,
  ProgressInterval,
  SetPartitioning,
  RoundedTarget,
  TimeLimit,
  MaxIterations,
};

struct SolverProgress {
  SolverCheckpointReason reason{SolverCheckpointReason::Initialization};
  std::size_t iteration{};
  std::size_t interval_first_iteration{};
  std::size_t interval_iterations{};
  std::optional<double> current_objective_h;
  std::optional<double> current_generalized_cost_h;
  std::optional<double> best_objective_h;
  SearchProfile search_profile{SearchProfile::Paper};
  bool current_feasible{};
  std::size_t current_routes{};
  double penalty_lambda{};
  std::size_t pool_size{};
  std::size_t pool_distinct_coverages{};
  std::size_t backup_pool_size{};
  bool pool_duals_valid{};
  std::uint64_t pool_effective_generation{};
  PoolStatistics pool_statistics;
  std::size_t sp_calls{};
  std::size_t adaptive_sp_calls{};
  std::size_t last_sp_iteration{};
  std::size_t last_best_improvement_iteration{};
  std::size_t last_coverage_improvement_iteration{};
  std::size_t lp_calls{};
  std::size_t frvcp_calls{};
  std::size_t frvcp_solve_calls{};
  std::size_t frvcp_cache_hits{};
  std::size_t exact_plan_cache_size{};
  std::size_t exact_plan_cache_hits{};
  std::size_t exact_plan_cache_misses{};
  std::size_t exact_plan_cache_evictions{};
  std::size_t path_option_cache_size{};
  std::size_t path_option_cache_hits{};
  std::size_t path_option_cache_misses{};
  std::size_t path_option_cache_evictions{};
  std::size_t path_option_first_builds{};
  std::size_t path_option_rebuilds{};
  std::uint64_t path_option_candidates_generated{};
  std::uint64_t path_option_nondominated_generated{};
  std::size_t virtual_completion_cache_size{};
  std::size_t virtual_completion_cache_hits{};
  std::size_t virtual_completion_cache_misses{};
  std::size_t virtual_completion_cache_evictions{};
  std::uint64_t vnd_linked_completion_attempts{};
  std::uint64_t vnd_linked_completion_changed_attempts{};
  std::uint64_t vnd_adaptive_customer_advances{};
  std::uint64_t vnd_adaptive_replace_path_selections{};
  std::uint64_t vnd_adaptive_inserted_replace_path_failures{};
  std::uint64_t vnd_adaptive_final_replace_path_calls{};
  std::uint64_t vnd_adaptive_final_replace_path_gaps_scanned{};
  std::uint64_t vnd_adaptive_final_replace_path_accepted{};
  InsertionCandidateStatistics initial_insertion_statistics;
  InsertionCandidateStatistics repair_insertion_statistics;
  double elapsed_s{};
  double interval_elapsed_s{};
  double initialization_s{};
  double perturbation_s{};
  double vnd_s{};
  double assembler_s{};
  double sp_model_s{};
  double lp_model_s{};
  double frvcp_model_s{};
  std::size_t evaluation_cache_size{};
  std::size_t evaluation_cache_hits{};
  std::size_t evaluation_cache_misses{};
  std::size_t evaluation_cache_evictions{};
  std::size_t interval_evaluation_cache_hits{};
  std::size_t interval_evaluation_cache_misses{};
  std::size_t interval_evaluation_cache_evictions{};
};

enum class ExactStage { SetPartitioning, LinearRelaxation, Frvcp };
enum class ExactStageBoundary { Start, ModelBuilt, Progress, End };

struct ExactStageProgress {
  ExactStage stage{ExactStage::SetPartitioning};
  ExactStageBoundary boundary{ExactStageBoundary::Start};
  std::size_t iteration{};
  std::size_t pool_columns{};
  std::size_t pool_distinct_coverages{};
  std::uint64_t pool_effective_generation{};
  std::size_t selected_routes{};
  std::size_t route_index{};  // One-based for FRVCP; zero otherwise.
  std::size_t route_customer_count{};
  std::size_t route_gap_count{};
  std::optional<double> remaining_time_s;
  std::optional<double> applied_time_limit_s;
  std::optional<double> objective_h;
  std::optional<double> assembled_objective_h;
  std::optional<double> best_bound_h;
  std::optional<double> relative_gap;
  std::optional<double> mip_start_objective_h;
  std::optional<int> status;
  std::string status_name;
  std::string mip_start_source;
  std::string sp_trigger_source;
  double model_runtime_s{};
  double wall_runtime_s{};
  double path_generation_runtime_s{};
  double model_build_runtime_s{};
  double explored_nodes{};
  int solution_count{};
  std::size_t mip_start_routes{};
  bool cache_hit{};
  bool rounded_target_reached{};
  bool incumbent_stall_limit_reached{};
  bool mip_start_supplied{};
  bool analytic_optimum{};
  bool lower_bound_proved_infeasible{};
  std::optional<double> objective_lower_bound_h;
  GurobiModelSize model_size;
  FrvcpNetworkSize network_size;
};

struct SolverResult {
  std::optional<Solution> best_solution;
  std::size_t iterations{};
  double runtime_s{};
  double final_penalty_lambda{};
  std::size_t pool_size{};
  std::size_t pool_distinct_coverages{};
  std::size_t backup_pool_size{};
  bool pool_duals_valid{};
  std::uint64_t pool_effective_generation{};
  PoolStatistics pool_statistics;
  std::size_t sp_calls{};
  std::size_t adaptive_sp_calls{};
  std::size_t lp_calls{};
  std::size_t frvcp_calls{};
  std::size_t frvcp_solve_calls{};
  std::size_t frvcp_cache_hits{};
  std::size_t exact_plan_cache_size{};
  std::size_t exact_plan_cache_hits{};
  std::size_t exact_plan_cache_misses{};
  std::size_t exact_plan_cache_evictions{};
  std::size_t path_option_cache_size{};
  std::size_t path_option_cache_hits{};
  std::size_t path_option_cache_misses{};
  std::size_t path_option_cache_evictions{};
  std::size_t path_option_first_builds{};
  std::size_t path_option_rebuilds{};
  std::uint64_t path_option_candidates_generated{};
  std::uint64_t path_option_nondominated_generated{};
  std::size_t virtual_completion_cache_size{};
  std::size_t virtual_completion_cache_hits{};
  std::size_t virtual_completion_cache_misses{};
  std::size_t virtual_completion_cache_evictions{};
  std::size_t evaluation_cache_hits{};
  std::size_t evaluation_cache_misses{};
  std::size_t evaluation_cache_evictions{};
  double initialization_runtime_s{};
  double perturbation_runtime_s{};
  double vnd_runtime_s{};
  double assembler_runtime_s{};
  double sp_model_runtime_s{};
  double lp_model_runtime_s{};
  double frvcp_model_runtime_s{};
  InsertionCandidateStatistics initial_insertion_statistics;
  InsertionCandidateStatistics repair_insertion_statistics;
  VndStatistics vnd_statistics;
  std::string stop_reason;
};

class IlsSpSolver {
 public:
  IlsSpSolver(const Instance& instance, GRBEnv& environment,
              SolverConfig config = {});

  [[nodiscard]] SolverResult solve(
      std::function<void(const SolverProgress&)> progress = {},
      std::function<void(const ExactStageProgress&)> exact_progress = {});

 private:
  [[nodiscard]] std::optional<Solution> assemble_routes(
      const Solution& local, const std::optional<Solution>& best,
      std::size_t iteration, GurobiDeadline deadline,
      const std::function<void(const ExactStageProgress&)>& exact_progress);
  [[nodiscard]] bool pool_covers_all_customers(
      const std::vector<RouteColumn>& columns) const;
  [[nodiscard]] bool rounded_target_reached(
      const Solution& candidate) const;
  void observe_penalty(bool feasible);
  [[nodiscard]] bool update_best(const Solution& candidate,
                                 std::optional<Solution>& best) const;

  const Instance& instance_;
  SolverConfig config_;
  std::mt19937_64 random_;
  RouteEvaluator evaluator_;
  PlanFactory factory_;
  PathSampler path_sampler_;
  InitialSolutionBuilder initial_builder_;
  Perturbation perturbation_;
  VariableNeighborhoodDescent vnd_;
  RoutePool pool_;
  GurobiSetPartitioning set_partitioning_;
  GurobiFrvcp frvcp_;
  std::vector<std::vector<int>> previous_sp_customer_sets_;
  double penalty_lambda_{};
  std::size_t consecutive_feasible_{0};
  std::size_t consecutive_infeasible_{0};
  std::size_t sp_calls_{0};
  std::size_t adaptive_sp_calls_{0};
  std::size_t last_sp_iteration_{0};
  std::uint64_t pool_effective_generation_at_last_sp_{0};
  std::size_t last_best_improvement_iteration_{0};
  std::size_t last_coverage_improvement_iteration_{0};
  std::size_t lp_calls_{0};
  std::size_t frvcp_calls_{0};
  double sp_model_runtime_s_{0.0};
  double lp_model_runtime_s_{0.0};
  double frvcp_model_runtime_s_{0.0};
};

}  // namespace ils_sp
