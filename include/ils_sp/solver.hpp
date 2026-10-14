#pragma once

#include <functional>

#include "ils_sp/copt_models.hpp"
#include "ils_sp/detail/frvcp_active_budget.hpp"
#include "ils_sp/search.hpp"

namespace ils_sp {

enum class ContinuationPolicy { Operational, Record };

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
  // Record follows TRC Algorithm 1's assembler acceptance. Algorithm 2's S
  // is shared with the caller; its generator assignment keeps the perturbed
  // base if assembly is rejected (authorized Class A state-passing choice).
  // The undefined non-SP return passes through the generator's VND output.
  ContinuationPolicy continuation_policy{ContinuationPolicy::Record};
  bool complete_path_catalog{true};
  bool implicit_replace_path{false};
  // Single proposal only; the four-proposal A experiment was rejected at G1.
  std::size_t replace_path_attempts{1};
  // A-class choice only for the final oversize RandomRoute/Target route.
  PartialRouteSubset partial_route_subset{PartialRouteSubset::Prefix};
  // Class A for Paper construction/repair only; VND retains direct-gap
  // best-improvement ranking. Other profiles define their own assignments.
  std::optional<bool> complete_insertion_paths;
  // Optional A interpretation for initial construction only. Unset preserves
  // the existing insertion/profile policy; it never changes LNS or VND.
  std::optional<bool> complete_initialization_paths;
  // Explicit initialization-only alternative: rebuild all station choices
  // from the direct customer sequence for each insertion candidate.
  bool rebuild_initialization_paths{false};
  // Initialization-only physical insertion; complete only new singleton routes.
  bool preserve_initialization_station_visits{false};
  // Explicit LNS interpretation: customer removal/insertion preserves the
  // physical station visits. No insertion-time path completion; VND keeps
  // its separately specified customer-gap and promising-arc contracts.
  bool preserve_station_visits{false};
  // Explicit A interpretation for Relocate route opening: each child keeps
  // the parent's station visits while retaining only its assigned customers.
  bool project_opening_station_visits{false};
  // Class A: Swap customer IDs in their original physical slots; CS stay put.
  bool preserve_swap_station_visits{false};
  std::size_t evaluation_cache_capacity{200'000};
  std::size_t exact_cache_capacity{20'000};
  // Legacy name: bounds the combined restricted-model optimal/infeasible memo.
  std::size_t infeasible_cache_capacity{20'000};
  CoptConfig set_partitioning;
  CoptConfig frvcp;
  std::optional<double> global_time_limit_s;
  // Optional solve-wide active FRVCP wall allowance; SP/LP are unchanged.
  std::optional<double> frvcp_active_wall_budget_s;
  std::optional<double> rounded_target_objective_h;
  std::optional<double> sp_incumbent_stall_limit_s;
  std::optional<double> sp_maximum_time_s;
  PoolIdentity pool_identity{PoolIdentity::OrderedSequence};
  // Explicit non-paper pool extensions; disabled in reproduction runs.
  bool protect_best_partition{false};
  bool keep_evicted_in_backup{false};
  bool reprice_after_lp{false};
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
  ReplacementSummaryStatistics replacement_summary_statistics;
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
  std::vector<int> route_visits;  // FRVCP input, excluding depot endpoints.
  std::string backend_phase;
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
  bool lp_reused_prepass{};
  double lp_prepass_wall_s{};
  std::optional<double> certificate_gap_upper_h;
  std::optional<double> certificate_acceptance_gap_h;
  std::size_t rc_fixed_columns{};
  std::size_t rc_considered_columns{};
  double model_runtime_s{};
  double wall_runtime_s{};
  double path_generation_runtime_s{};
  double model_build_runtime_s{};
  std::optional<double> explored_nodes;
  int solution_count{};
  std::size_t mip_start_routes{};
  bool cache_hit{};
  bool rounded_target_reached{};
  bool incumbent_stall_limit_reached{};
  bool mip_start_supplied{};
  bool analytic_optimum{};
  bool fixed_order_optimal{};
  bool lower_bound_proved_infeasible{};
  std::optional<double> objective_lower_bound_h;
  CoptModelSize model_size;
  FrvcpNetworkSize network_size;
};

// Emitted only after the existing feasible record-improvement gate succeeds.
// Initialization is not a record: Algorithm 1 starts Sbest empty. Its first
// event is the first accepted candidate, at an iteration >= 1. Elapsed time
// starts at solve() entry, including initialization, excluding construction.
struct RecordImprovementProgress {
  std::size_t iteration{};
  double objective_h{};
  double elapsed_s{};
};

struct SolverResult {
  std::optional<Solution> best_solution;
  // Output-only safety snapshot; never admitted to search or reported as a record.
  std::optional<Solution> initialization_fallback;
  // Terminal search state, which may be infeasible; never a quality record.
  Solution final_solution;
  std::vector<int> pool_uncovered_customers;
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
  std::size_t sp_skipped_incomplete_pool{};
  std::size_t adaptive_sp_calls{};
  std::size_t lp_calls{};
  std::size_t sp_lp_prepass_calls{};
  std::size_t sp_lp_certified_calls{};
  std::size_t sp_lp_fallback_calls{};
  std::size_t sp_lp_exception_calls{};
  double sp_lp_prepass_wall_s{};
  std::size_t sp_rc_fixing_calls{};
  std::size_t sp_rc_fixed_columns{};
  std::size_t sp_rc_considered_columns{};
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
  ReplacementSummaryStatistics replacement_summary_statistics;
  std::size_t continuation_exploration_entries{};
  std::size_t continuation_exploration_iterations{};
  std::size_t continuation_record_exits{};
  std::size_t continuation_epoch_exits{};
  double frvcp_active_wall_spent_s{};
  std::size_t frvcp_budget_skipped_routes{};
  bool frvcp_budget_exhausted{};
  std::string stop_reason;
};

class IlsSpSolver {
 public:
  IlsSpSolver(const Instance& instance, Envr& environment,
              SolverConfig config = {});

  [[nodiscard]] SolverResult solve(
      std::function<void(const SolverProgress&)> progress = {},
      std::function<void(const ExactStageProgress&)> exact_progress = {},
      std::function<void(const RecordImprovementProgress&)> record_progress = {});

 private:
  [[nodiscard]] std::optional<Solution> assemble_routes(
      Solution& local, const std::optional<Solution>& best,
      std::size_t iteration, CoptDeadline deadline,
      detail::FrvcpActiveWallBudget& frvcp_budget,
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
  CoptSetPartitioning set_partitioning_;
  CoptFrvcp frvcp_;
  std::vector<std::vector<int>> previous_sp_customer_sets_;
  double penalty_lambda_{};
  std::size_t consecutive_feasible_{0};
  std::size_t consecutive_infeasible_{0};
  std::size_t sp_calls_{0};
  std::size_t sp_skipped_incomplete_pool_{0};
  std::size_t adaptive_sp_calls_{0};
  std::size_t last_sp_iteration_{0};
  std::uint64_t pool_effective_generation_at_last_sp_{0};
  std::size_t last_best_improvement_iteration_{0};
  std::size_t last_coverage_improvement_iteration_{0};
  std::size_t lp_calls_{0};
  std::size_t sp_lp_prepass_calls_{0};
  std::size_t sp_lp_certified_calls_{0};
  std::size_t sp_lp_fallback_calls_{0};
  std::size_t sp_lp_exception_calls_{0};
  double sp_lp_prepass_wall_s_{0.0};
  std::size_t sp_rc_fixing_calls_{0};
  std::size_t sp_rc_fixed_columns_{0};
  std::size_t sp_rc_considered_columns_{0};
  std::size_t frvcp_calls_{0};
  double sp_model_runtime_s_{0.0};
  double lp_model_runtime_s_{0.0};
  double frvcp_model_runtime_s_{0.0};
};

}  // namespace ils_sp
