#include "ils_sp/solver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_set>

namespace ils_sp {
namespace {

using SolverClock = std::chrono::steady_clock;

[[nodiscard]] double elapsed_seconds(SolverClock::time_point start,
                                     SolverClock::time_point end) {
  return std::chrono::duration<double>(end - start).count();
}

[[nodiscard]] std::optional<double> remaining_seconds(
    GurobiDeadline deadline) {
  if (!deadline.has_value()) return std::nullopt;
  return std::max(
      0.0,
      std::chrono::duration<double>(*deadline - SolverClock::now()).count());
}

[[nodiscard]] bool deadline_expired(GurobiDeadline deadline) {
  return deadline.has_value() && SolverClock::now() >= *deadline;
}

struct EvaluationCacheSnapshot {
  std::size_t hits{};
  std::size_t misses{};
  std::size_t evictions{};
};

struct SpMipStartChoice {
  std::vector<std::size_t> column_indices;
  std::string source{"none"};
  std::optional<double> objective_h;
};

[[nodiscard]] std::optional<SpMipStartChoice> make_sp_mip_start(
    const Instance& instance, const std::vector<RouteColumn>& columns,
    const std::vector<std::vector<int>>& requested_customer_sets,
    std::string source) {
  if (requested_customer_sets.empty()) return std::nullopt;
  SpMipStartChoice choice;
  choice.source = std::move(source);
  choice.column_indices.reserve(requested_customer_sets.size());
  double objective_h = 0.0;
  std::unordered_set<int> covered;
  for (const std::vector<int>& customer_set : requested_customer_sets) {
    const auto column = std::lower_bound(
        columns.begin(), columns.end(), customer_set,
        [](const RouteColumn& candidate, const std::vector<int>& key) {
          return candidate.customer_set < key;
        });
    if (column == columns.end() || column->customer_set != customer_set) {
      return std::nullopt;
    }
    for (const int customer_id : customer_set) {
      if (!covered.insert(customer_id).second) return std::nullopt;
    }
    choice.column_indices.push_back(
        static_cast<std::size_t>(std::distance(columns.begin(), column)));
    objective_h += column->cost_h;
  }
  if (covered.size() != instance.customer_ids().size() ||
      !std::all_of(instance.customer_ids().begin(),
                   instance.customer_ids().end(),
                   [&](int customer_id) {
                     return covered.contains(customer_id);
                   })) {
    return std::nullopt;
  }
  choice.objective_h = objective_h;
  return choice;
}

[[nodiscard]] std::vector<std::vector<int>> solution_customer_sets(
    const Solution& solution) {
  std::vector<std::vector<int>> result;
  result.reserve(solution.plans.size());
  for (const Plan& plan : solution.plans) {
    result.push_back(sorted_customer_set(plan.customer_ids));
  }
  return result;
}

[[nodiscard]] EvaluationCacheSnapshot evaluation_cache_snapshot(
    const RouteEvaluator& evaluator) {
  return EvaluationCacheSnapshot{.hits = evaluator.cache_hits(),
                                 .misses = evaluator.cache_misses(),
                                 .evictions = evaluator.cache_evictions()};
}

}  // namespace

IlsSpSolver::IlsSpSolver(const Instance& instance, GRBEnv& environment,
                         SolverConfig config)
    : instance_(instance),
      config_(std::move(config)),
      random_(config_.seed),
      evaluator_(instance_, config_.evaluation_cache_capacity),
      factory_(instance_, evaluator_, config_.exact_cache_capacity),
      path_sampler_(instance_, factory_, random_),
      initial_builder_(instance_, factory_, path_sampler_, random_,
                       config_.search_profile, config_.promising_arc_count,
                       config_.granular_threshold),
      perturbation_(instance_, factory_, path_sampler_, random_,
                    config_.search_profile, config_.promising_arc_count,
                    config_.granular_threshold),
      vnd_(instance_, factory_, path_sampler_, random_,
           config_.promising_arc_count, config_.search_profile),
      pool_(config_.paper.pool_capacity),
      set_partitioning_(instance_, environment, config_.set_partitioning),
      frvcp_(instance_, evaluator_, factory_, environment, config_.frvcp,
             config_.infeasible_cache_capacity),
      penalty_lambda_(config_.paper.lambda_initial) {
  const PaperParameters& paper = config_.paper;
  if (paper.penalty_window == 0 || paper.penalty_factor <= 1.0 ||
      paper.lambda_min < 0.0 || paper.lambda_min > paper.lambda_initial ||
      paper.lambda_initial > paper.lambda_max || paper.sp_frequency == 0 ||
      paper.backup_frequency == 0 || paper.pool_capacity == 0 ||
      paper.max_iterations == 0 || config_.progress_interval == 0 ||
      (config_.global_time_limit_s.has_value() &&
       (!std::isfinite(*config_.global_time_limit_s) ||
        *config_.global_time_limit_s <= 0.0)) ||
      (config_.rounded_target_objective_h.has_value() &&
       (!std::isfinite(*config_.rounded_target_objective_h) ||
        *config_.rounded_target_objective_h < 0.0)) ||
      (config_.sp_incumbent_stall_limit_s.has_value() &&
       (!std::isfinite(*config_.sp_incumbent_stall_limit_s) ||
        *config_.sp_incumbent_stall_limit_s <= 0.0)) ||
      (config_.sp_maximum_time_s.has_value() &&
       (!std::isfinite(*config_.sp_maximum_time_s) ||
        *config_.sp_maximum_time_s <= 0.0)) ||
      (config_.adaptive_sp_stall_iterations.has_value() &&
       (*config_.adaptive_sp_stall_iterations == 0 ||
        config_.adaptive_sp_minimum_interval == 0))) {
    throw std::invalid_argument("invalid ILS-SP solver configuration");
  }
}

void IlsSpSolver::observe_penalty(bool feasible) {
  if (feasible) {
    ++consecutive_feasible_;
    consecutive_infeasible_ = 0;
  } else {
    ++consecutive_infeasible_;
    consecutive_feasible_ = 0;
  }
  if (consecutive_infeasible_ >= config_.paper.penalty_window) {
    penalty_lambda_ = std::min(config_.paper.lambda_max,
                               penalty_lambda_ * config_.paper.penalty_factor);
    consecutive_infeasible_ = 0;
  } else if (consecutive_feasible_ >= config_.paper.penalty_window) {
    penalty_lambda_ = std::max(config_.paper.lambda_min,
                               penalty_lambda_ / config_.paper.penalty_factor);
    consecutive_feasible_ = 0;
  }
}

bool IlsSpSolver::update_best(const Solution& candidate,
                              std::optional<Solution>& best) const {
  if (!candidate.feasible()) return false;
  if (!best.has_value() ||
      candidate.raw_cost_h() < best->raw_cost_h() - kCostTolerance) {
    best = candidate;
    return true;
  }
  return false;
}

bool IlsSpSolver::pool_covers_all_customers(
    const std::vector<RouteColumn>& columns) const {
  std::unordered_set<int> covered;
  for (const RouteColumn& column : columns) {
    covered.insert(column.customer_set.begin(), column.customer_set.end());
  }
  return std::all_of(instance_.customer_ids().begin(),
                     instance_.customer_ids().end(),
                     [&](int customer_id) { return covered.contains(customer_id); });
}

bool IlsSpSolver::rounded_target_reached(const Solution& candidate) const {
  return config_.rounded_target_objective_h.has_value() &&
         candidate.feasible() &&
         std::llround(candidate.raw_cost_h() * 100.0) <=
             std::llround(*config_.rounded_target_objective_h * 100.0);
}

std::optional<Solution> IlsSpSolver::assemble_routes(
    const Solution& local, const std::optional<Solution>& best,
    std::size_t iteration, GurobiDeadline deadline,
    const std::function<void(const ExactStageProgress&)>& exact_progress) {
  for (const Plan& plan : local.plans) {
    if (plan.evaluation->feasible) {
      const PoolUpdate update = pool_.admit(RouteColumn::from_plan(plan));
      if (update.effective_coverage_improvement) {
        last_coverage_improvement_iteration_ = iteration;
      }
    }
  }
  const bool adaptive_enabled =
      config_.adaptive_sp_stall_iterations.has_value();
  // The paper frequency is an absolute cadence.  An adaptive call is an
  // additional engineering trigger and must never shift the next paper call.
  const bool periodic_sp = iteration % config_.paper.sp_frequency == 0;
  const std::size_t periodic_remainder =
      iteration % config_.paper.sp_frequency;
  const std::size_t iterations_until_periodic =
      periodic_remainder == 0
          ? 0
          : config_.paper.sp_frequency - periodic_remainder;
  const bool pool_advanced_since_last_sp =
      pool_.effective_generation() > pool_effective_generation_at_last_sp_;
  const bool adaptive_sp =
      adaptive_enabled && !periodic_sp && best.has_value() &&
      pool_.main_size() == pool_.capacity() &&
      pool_advanced_since_last_sp &&
      iteration - last_sp_iteration_ >=
          config_.adaptive_sp_minimum_interval &&
      iterations_until_periodic >= config_.adaptive_sp_minimum_interval &&
      iteration - last_best_improvement_iteration_ >=
          *config_.adaptive_sp_stall_iterations &&
      iteration - last_coverage_improvement_iteration_ >=
          *config_.adaptive_sp_stall_iterations;
  const bool run_sp = periodic_sp || adaptive_sp;
  if (iteration % config_.paper.backup_frequency == 0 || run_sp) {
    const PoolRepriceUpdate repriced = pool_.reprice_backup();
    if (repriced.effective_coverage_improvements != 0) {
      last_coverage_improvement_iteration_ = iteration;
    }
  }
  if (!run_sp) return std::nullopt;
  if (deadline_expired(deadline)) return std::nullopt;
  const std::vector<RouteColumn> columns = pool_.main_columns();
  if (deadline_expired(deadline) || !pool_covers_all_customers(columns)) {
    return std::nullopt;
  }

  SpMipStartChoice mip_start;
  const auto consider_mip_start = [&](std::optional<SpMipStartChoice> choice) {
    if (!choice.has_value()) return;
    if (mip_start.column_indices.empty() ||
        (choice->objective_h.has_value() &&
         (!mip_start.objective_h.has_value() ||
          *choice->objective_h <
              *mip_start.objective_h - kCostTolerance))) {
      mip_start = std::move(*choice);
    }
  };
  // Warm starts do not alter the SP model.  Keep the previous-SP partition as
  // the stable tie preference, but do not discard a cheaper representable
  // global record or local feasible partition.
  consider_mip_start(make_sp_mip_start(
      instance_, columns, previous_sp_customer_sets_, "previous_sp"));
  if (best.has_value()) {
    consider_mip_start(make_sp_mip_start(
        instance_, columns, solution_customer_sets(*best), "best_feasible"));
  }
  if (local.feasible()) {
    consider_mip_start(make_sp_mip_start(
        instance_, columns, solution_customer_sets(local), "local_feasible"));
  }
  if (exact_progress) {
    ExactStageProgress checkpoint;
    checkpoint.stage = ExactStage::SetPartitioning;
    checkpoint.boundary = ExactStageBoundary::Start;
    checkpoint.iteration = iteration;
    checkpoint.pool_columns = columns.size();
    checkpoint.pool_distinct_coverages = pool_.distinct_main_coverages();
    checkpoint.pool_effective_generation = pool_.effective_generation();
    checkpoint.mip_start_routes = mip_start.column_indices.size();
    checkpoint.mip_start_objective_h = mip_start.objective_h;
    checkpoint.mip_start_source = mip_start.source;
    checkpoint.sp_trigger_source =
        adaptive_sp && !periodic_sp ? "pool_stall" : "frequency";
    checkpoint.remaining_time_s = remaining_seconds(deadline);
    exact_progress(checkpoint);
  }
  ++sp_calls_;
  last_sp_iteration_ = iteration;
  pool_effective_generation_at_last_sp_ = pool_.effective_generation();
  if (adaptive_sp && !periodic_sp) ++adaptive_sp_calls_;
  SetPartitioningProgressCallback model_progress;
  if (exact_progress) {
    model_progress = [&](const SetPartitioningProgress& progress) {
      ExactStageProgress checkpoint;
      checkpoint.stage = ExactStage::SetPartitioning;
      checkpoint.boundary =
          progress.event == SetPartitioningProgressEvent::ModelBuilt
              ? ExactStageBoundary::ModelBuilt
              : ExactStageBoundary::Progress;
      checkpoint.iteration = iteration;
      checkpoint.pool_columns = columns.size();
      checkpoint.pool_distinct_coverages = pool_.distinct_main_coverages();
      checkpoint.pool_effective_generation = pool_.effective_generation();
      checkpoint.remaining_time_s = remaining_seconds(deadline);
      checkpoint.objective_h = progress.incumbent_objective_h;
      checkpoint.best_bound_h = progress.best_bound_h;
      checkpoint.relative_gap = progress.relative_gap;
      checkpoint.model_runtime_s = progress.runtime_s;
      checkpoint.model_build_runtime_s = progress.model_build_runtime_s;
      checkpoint.explored_nodes = progress.explored_nodes;
      checkpoint.solution_count = progress.solution_count;
      checkpoint.mip_start_routes = mip_start.column_indices.size();
      checkpoint.mip_start_objective_h = mip_start.objective_h;
      checkpoint.mip_start_source = mip_start.source;
      checkpoint.sp_trigger_source =
          adaptive_sp && !periodic_sp ? "pool_stall" : "frequency";
      checkpoint.rounded_target_reached =
          progress.rounded_target_reached;
      checkpoint.incumbent_stall_limit_reached =
          progress.incumbent_stall_limit_reached;
      checkpoint.model_size = progress.model_size;
      exact_progress(checkpoint);
    };
  }
  SetPartitioningResult partition = set_partitioning_.solve_integer(
      columns,
      SetPartitioningSolveOptions{
          .deadline = deadline,
          .mip_start_indices = mip_start.column_indices,
          .rounded_target_objective_h =
              config_.rounded_target_objective_h,
          .incumbent_stall_limit_s =
              config_.sp_incumbent_stall_limit_s,
          .maximum_wall_time_s = config_.sp_maximum_time_s,
          .progress = std::move(model_progress)});
  sp_model_runtime_s_ += partition.runtime_s;
  if (!partition.selected_columns.empty()) {
    previous_sp_customer_sets_.clear();
    previous_sp_customer_sets_.reserve(partition.selected_columns.size());
    for (const RouteColumn& selected : partition.selected_columns) {
      previous_sp_customer_sets_.push_back(selected.customer_set);
    }
  }
  if (exact_progress) {
    ExactStageProgress checkpoint;
    checkpoint.stage = ExactStage::SetPartitioning;
    checkpoint.boundary = ExactStageBoundary::End;
    checkpoint.iteration = iteration;
    checkpoint.pool_columns = columns.size();
    checkpoint.pool_distinct_coverages = pool_.distinct_main_coverages();
    checkpoint.pool_effective_generation = pool_.effective_generation();
    checkpoint.selected_routes = partition.selected_columns.size();
    checkpoint.remaining_time_s = remaining_seconds(deadline);
    checkpoint.applied_time_limit_s = partition.applied_time_limit_s;
    checkpoint.objective_h = partition.objective_h;
    checkpoint.best_bound_h = partition.best_bound_h;
    checkpoint.relative_gap = partition.relative_gap;
    checkpoint.mip_start_objective_h = partition.mip_start_objective_h;
    checkpoint.status = partition.status;
    checkpoint.status_name = partition.status_name;
    checkpoint.mip_start_source = mip_start.source;
    checkpoint.sp_trigger_source =
        adaptive_sp && !periodic_sp ? "pool_stall" : "frequency";
    checkpoint.model_runtime_s = partition.runtime_s;
    checkpoint.wall_runtime_s = partition.wall_runtime_s;
    checkpoint.model_build_runtime_s = partition.model_build_runtime_s;
    checkpoint.explored_nodes = partition.explored_nodes;
    checkpoint.solution_count = partition.solution_count;
    checkpoint.mip_start_routes = partition.mip_start_routes;
    checkpoint.rounded_target_reached =
        partition.rounded_target_reached;
    checkpoint.incumbent_stall_limit_reached =
        partition.incumbent_stall_limit_reached;
    checkpoint.model_size = partition.model_size;
    exact_progress(checkpoint);
  }

  std::optional<Solution> assembled;
  if (!partition.selected_columns.empty()) {
    assembled.emplace();
    assembled->plans.reserve(partition.selected_columns.size());
    for (const RouteColumn& selected : partition.selected_columns) {
      assembled->plans.push_back(selected.plan);
    }
    assembled->validate_partition(instance_);
    // A rounded target is an explicit engineering stop.  Paper-default runs
    // have no target and therefore continue with LP and all FRVCP calls.
    if (rounded_target_reached(*assembled)) return assembled;
  }

  if (deadline_expired(deadline)) return assembled;
  if (exact_progress) {
    ExactStageProgress checkpoint;
    checkpoint.stage = ExactStage::LinearRelaxation;
    checkpoint.boundary = ExactStageBoundary::Start;
    checkpoint.iteration = iteration;
    checkpoint.pool_columns = columns.size();
    checkpoint.selected_routes = partition.selected_columns.size();
    checkpoint.remaining_time_s = remaining_seconds(deadline);
    exact_progress(checkpoint);
  }

  ++lp_calls_;
  LinearRelaxationResult relaxation =
      set_partitioning_.solve_relaxation(columns, deadline);
  lp_model_runtime_s_ += relaxation.runtime_s;
  if (exact_progress) {
    ExactStageProgress checkpoint;
    checkpoint.stage = ExactStage::LinearRelaxation;
    checkpoint.boundary = ExactStageBoundary::End;
    checkpoint.iteration = iteration;
    checkpoint.pool_columns = columns.size();
    checkpoint.selected_routes = partition.selected_columns.size();
    checkpoint.remaining_time_s = remaining_seconds(deadline);
    checkpoint.applied_time_limit_s = relaxation.applied_time_limit_s;
    checkpoint.objective_h = relaxation.objective_h;
    checkpoint.status = relaxation.status;
    checkpoint.status_name = relaxation.status_name;
    checkpoint.model_runtime_s = relaxation.runtime_s;
    checkpoint.wall_runtime_s = relaxation.wall_runtime_s;
    checkpoint.model_size = relaxation.model_size;
    exact_progress(checkpoint);
  }
  if (relaxation.status == GRB_OPTIMAL) {
    pool_.update_duals(std::move(relaxation.duals));
    const PoolRepriceUpdate repriced = pool_.reprice_backup();
    if (repriced.effective_coverage_improvements != 0) {
      last_coverage_improvement_iteration_ = iteration;
    }
  }
  if (!assembled.has_value() || deadline_expired(deadline)) return assembled;

  for (std::size_t index = 0; index < partition.selected_columns.size();
       ++index) {
    if (deadline_expired(deadline) || rounded_target_reached(*assembled)) break;
    const RouteColumn& selected = partition.selected_columns[index];
    if (exact_progress) {
      ExactStageProgress checkpoint;
      checkpoint.stage = ExactStage::Frvcp;
      checkpoint.boundary = ExactStageBoundary::Start;
      checkpoint.iteration = iteration;
      checkpoint.pool_columns = columns.size();
      checkpoint.selected_routes = partition.selected_columns.size();
      checkpoint.route_index = index + 1;
      checkpoint.route_customer_count = selected.plan.customer_ids.size();
      checkpoint.route_gap_count = selected.plan.customer_ids.size() + 1;
      checkpoint.remaining_time_s = remaining_seconds(deadline);
      checkpoint.assembled_objective_h = assembled->raw_cost_h();
      checkpoint.mip_start_objective_h =
          selected.plan.evaluation->raw_cost_h;
      checkpoint.mip_start_source = "sp_selected";
      exact_progress(checkpoint);
    }
    ++frvcp_calls_;
    FrvcpResult optimized = frvcp_.optimize(selected.plan, deadline);
    frvcp_model_runtime_s_ += optimized.runtime_s;
    if (optimized.plan.has_value() &&
        optimized.plan->evaluation->feasible) {
      const double incumbent_cost =
          assembled->plans[index].evaluation->raw_cost_h;
      const double candidate_cost =
          optimized.plan->evaluation->raw_cost_h;
      const bool improves =
          candidate_cost < incumbent_cost - kCostTolerance;
      const bool upgrades_exactness =
          optimized.plan->exact_charging &&
          !assembled->plans[index].exact_charging &&
          candidate_cost <= incumbent_cost + kCostTolerance;
      if (improves || upgrades_exactness) {
        assembled->plans[index] = *optimized.plan;
        const PoolUpdate update =
            pool_.admit(RouteColumn::from_plan(*optimized.plan));
        if (update.effective_coverage_improvement) {
          last_coverage_improvement_iteration_ = iteration;
        }
      }
    }
    const bool target_reached = rounded_target_reached(*assembled);
    if (exact_progress) {
      ExactStageProgress checkpoint;
      checkpoint.stage = ExactStage::Frvcp;
      checkpoint.boundary = ExactStageBoundary::End;
      checkpoint.iteration = iteration;
      checkpoint.pool_columns = columns.size();
      checkpoint.selected_routes = partition.selected_columns.size();
      checkpoint.route_index = index + 1;
      checkpoint.route_customer_count = selected.plan.customer_ids.size();
      checkpoint.route_gap_count = selected.plan.customer_ids.size() + 1;
      checkpoint.remaining_time_s = remaining_seconds(deadline);
      checkpoint.applied_time_limit_s = optimized.applied_time_limit_s;
      checkpoint.objective_h = optimized.model_objective_h;
      checkpoint.assembled_objective_h = assembled->raw_cost_h();
      checkpoint.status = optimized.status;
      checkpoint.status_name = optimized.status_name;
      checkpoint.model_runtime_s = optimized.runtime_s;
      checkpoint.wall_runtime_s = optimized.wall_runtime_s;
      checkpoint.path_generation_runtime_s =
          optimized.path_generation_runtime_s;
      checkpoint.model_build_runtime_s = optimized.model_build_runtime_s;
      checkpoint.cache_hit = optimized.cache_hit;
      checkpoint.rounded_target_reached = target_reached;
      checkpoint.mip_start_objective_h =
          optimized.mip_start_objective_h;
      checkpoint.mip_start_source = "sp_selected";
      checkpoint.mip_start_supplied = optimized.mip_start_supplied;
      checkpoint.analytic_optimum = optimized.analytic_optimum;
      checkpoint.lower_bound_proved_infeasible =
          optimized.lower_bound_proved_infeasible;
      checkpoint.objective_lower_bound_h =
          optimized.objective_lower_bound_h;
      checkpoint.model_size = optimized.model_size;
      checkpoint.network_size = optimized.network_size;
      exact_progress(checkpoint);
    }
    if (target_reached) break;
  }
  assembled->validate_partition(instance_);
  return assembled;
}

SolverResult IlsSpSolver::solve(
    std::function<void(const SolverProgress&)> progress,
    std::function<void(const ExactStageProgress&)> exact_progress) {
  const auto started = SolverClock::now();
  const auto deadline = config_.global_time_limit_s.has_value()
                            ? std::optional<SolverClock::time_point>(
                                  started + std::chrono::duration_cast<
                                                SolverClock::duration>(
                                                std::chrono::duration<double>(
                                                    *config_.global_time_limit_s)))
                            : std::nullopt;
  Solution current = initial_builder_.build(penalty_lambda_);
  const auto initialization_finished = SolverClock::now();
  const double initialization_runtime_s =
      elapsed_seconds(started, initialization_finished);
  std::optional<Solution> best;
  std::size_t iterations = 0;
  std::string stop_reason = "max_iterations";
  double perturbation_runtime_s = 0.0;
  double vnd_runtime_s = 0.0;
  double assembler_runtime_s = 0.0;

  EvaluationCacheSnapshot interval_cache_start{};
  const EvaluationCacheSnapshot initialization_cache =
      evaluation_cache_snapshot(evaluator_);
  if (progress) {
    progress(SolverProgress{
        .reason = SolverCheckpointReason::Initialization,
        .iteration = 0,
        .interval_first_iteration = 0,
        .interval_iterations = 0,
        .current_objective_h = current.raw_cost_h(),
        .current_generalized_cost_h =
            generalized_cost(current, instance_, penalty_lambda_),
        .best_objective_h = std::nullopt,
        .search_profile = config_.search_profile,
        .current_feasible = current.feasible(),
        .current_routes = current.plans.size(),
        .penalty_lambda = penalty_lambda_,
        .pool_size = pool_.main_size(),
        .pool_distinct_coverages = pool_.distinct_main_coverages(),
        .backup_pool_size = pool_.backup_size(),
        .pool_duals_valid = pool_.duals_valid(),
        .pool_effective_generation = pool_.effective_generation(),
        .pool_statistics = pool_.statistics(),
        .sp_calls = sp_calls_,
        .adaptive_sp_calls = adaptive_sp_calls_,
        .last_sp_iteration = last_sp_iteration_,
        .last_best_improvement_iteration =
            last_best_improvement_iteration_,
        .last_coverage_improvement_iteration =
            last_coverage_improvement_iteration_,
        .lp_calls = lp_calls_,
        .frvcp_calls = frvcp_calls_,
        .frvcp_solve_calls = frvcp_.solve_calls(),
        .frvcp_cache_hits = frvcp_.cache_hits(),
        .exact_plan_cache_size = factory_.exact_cache_size(),
        .exact_plan_cache_hits = factory_.exact_cache_hits(),
        .exact_plan_cache_misses = factory_.exact_cache_misses(),
        .exact_plan_cache_evictions = factory_.exact_cache_evictions(),
        .path_option_cache_size = path_sampler_.option_cache_size(),
        .path_option_cache_hits = path_sampler_.option_cache_hits(),
        .path_option_cache_misses = path_sampler_.option_cache_misses(),
        .path_option_cache_evictions = path_sampler_.option_cache_evictions(),
        .path_option_first_builds = path_sampler_.option_first_builds(),
        .path_option_rebuilds = path_sampler_.option_rebuilds(),
        .path_option_candidates_generated =
            path_sampler_.option_candidates_generated(),
        .path_option_nondominated_generated =
            path_sampler_.option_nondominated_generated(),
        .virtual_completion_cache_size =
            path_sampler_.virtual_completion_cache_size(),
        .virtual_completion_cache_hits =
            path_sampler_.virtual_completion_cache_hits(),
        .virtual_completion_cache_misses =
            path_sampler_.virtual_completion_cache_misses(),
        .virtual_completion_cache_evictions =
            path_sampler_.virtual_completion_cache_evictions(),
        .vnd_linked_completion_attempts =
            vnd_.statistics().linked_completion_attempts,
        .vnd_linked_completion_changed_attempts =
            vnd_.statistics().linked_completion_changed_attempts,
        .vnd_adaptive_customer_advances =
            vnd_.statistics().adaptive_customer_advances,
        .vnd_adaptive_replace_path_selections =
            vnd_.statistics().adaptive_replace_path_selections,
        .vnd_adaptive_inserted_replace_path_failures =
            vnd_.statistics().adaptive_inserted_replace_path_failures,
        .vnd_adaptive_final_replace_path_calls =
            vnd_.statistics().adaptive_final_replace_path_calls,
        .vnd_adaptive_final_replace_path_gaps_scanned =
            vnd_.statistics().adaptive_final_replace_path_gaps_scanned,
        .vnd_adaptive_final_replace_path_accepted =
            vnd_.statistics().adaptive_final_replace_path_accepted,
        .initial_insertion_statistics = initial_builder_.statistics(),
        .repair_insertion_statistics = perturbation_.statistics(),
        .elapsed_s = initialization_runtime_s,
        .interval_elapsed_s = initialization_runtime_s,
        .initialization_s = initialization_runtime_s,
        .perturbation_s = 0.0,
        .vnd_s = 0.0,
        .assembler_s = 0.0,
        .sp_model_s = 0.0,
        .lp_model_s = 0.0,
        .frvcp_model_s = 0.0,
        .evaluation_cache_size = evaluator_.cache_size(),
        .evaluation_cache_hits = initialization_cache.hits,
        .evaluation_cache_misses = initialization_cache.misses,
        .evaluation_cache_evictions = initialization_cache.evictions,
        .interval_evaluation_cache_hits = initialization_cache.hits,
        .interval_evaluation_cache_misses = initialization_cache.misses,
        .interval_evaluation_cache_evictions = initialization_cache.evictions});
  }

  auto interval_started = SolverClock::now();
  std::size_t interval_first_iteration = 1;
  double interval_perturbation_s = 0.0;
  double interval_vnd_s = 0.0;
  double interval_assembler_s = 0.0;
  interval_cache_start = initialization_cache;
  double interval_sp_model_start_s = sp_model_runtime_s_;
  double interval_lp_model_start_s = lp_model_runtime_s_;
  double interval_frvcp_model_start_s = frvcp_model_runtime_s_;

  for (std::size_t iteration = 1; iteration <= config_.paper.max_iterations;
       ++iteration) {
    if (deadline.has_value() && SolverClock::now() >= *deadline) {
      stop_reason = "time_limit";
      break;
    }
    const auto perturbation_started = SolverClock::now();
    Solution perturbed = perturbation_.apply(current, penalty_lambda_);
    const double perturbation_s =
        elapsed_seconds(perturbation_started, SolverClock::now());
    perturbation_runtime_s += perturbation_s;
    interval_perturbation_s += perturbation_s;

    const auto vnd_started = SolverClock::now();
    Solution local = vnd_.improve(std::move(perturbed), penalty_lambda_, deadline);
    const double vnd_s = elapsed_seconds(vnd_started, SolverClock::now());
    vnd_runtime_s += vnd_s;
    interval_vnd_s += vnd_s;

    observe_penalty(local.feasible());
    if (update_best(local, best)) {
      last_best_improvement_iteration_ = iteration;
      const std::uint64_t generation_before = pool_.effective_generation();
      pool_.protect_partition(*best);
      if (pool_.effective_generation() != generation_before) {
        last_coverage_improvement_iteration_ = iteration;
      }
    }
    const bool target_before_assembly =
        best.has_value() && rounded_target_reached(*best);
    const bool expired_before_assembly =
        deadline.has_value() && SolverClock::now() >= *deadline;
    if (!expired_before_assembly && !target_before_assembly) {
      const auto assembler_started = SolverClock::now();
      auto assembled =
          assemble_routes(local, best, iteration, deadline, exact_progress);
      const double assembler_s =
          elapsed_seconds(assembler_started, SolverClock::now());
      assembler_runtime_s += assembler_s;
      interval_assembler_s += assembler_s;
      if (assembled.has_value()) {
        current = std::move(*assembled);
        if (update_best(current, best)) {
          last_best_improvement_iteration_ = iteration;
          const std::uint64_t generation_before = pool_.effective_generation();
          pool_.protect_partition(*best);
          if (pool_.effective_generation() != generation_before) {
            last_coverage_improvement_iteration_ = iteration;
          }
        }
      } else {
        // Algorithm 5 only assigns S'' on SP iterations.  On ordinary
        // iterations S' must pass through so that the next perturbation starts
        // from the local solution, as made explicit by Algorithm 1 line 16 of
        // the predecessor HLNS method.
        current = std::move(local);
      }
      iterations = iteration;
    } else {
      current = std::move(local);
      if (target_before_assembly) iterations = iteration;
    }

    const bool deadline_reached =
        deadline.has_value() && SolverClock::now() >= *deadline;
    const bool target_reached =
        best.has_value() && rounded_target_reached(*best);

    std::optional<SolverCheckpointReason> checkpoint_reason;
    if (deadline_reached) {
      checkpoint_reason = SolverCheckpointReason::TimeLimit;
    } else if (target_reached) {
      checkpoint_reason = SolverCheckpointReason::RoundedTarget;
    } else if (iteration == config_.paper.max_iterations) {
      checkpoint_reason = SolverCheckpointReason::MaxIterations;
    } else if (sp_calls_ != 0 && last_sp_iteration_ == iteration) {
      checkpoint_reason = SolverCheckpointReason::SetPartitioning;
    } else if (iteration % config_.progress_interval == 0) {
      checkpoint_reason = SolverCheckpointReason::ProgressInterval;
    }
    if (progress && checkpoint_reason.has_value()) {
      const auto checkpoint_time = SolverClock::now();
      const EvaluationCacheSnapshot cache =
          evaluation_cache_snapshot(evaluator_);
      progress(SolverProgress{
          .reason = *checkpoint_reason,
          .iteration = iterations,
          .interval_first_iteration = interval_first_iteration,
          .interval_iterations =
              iterations >= interval_first_iteration
                  ? iterations - interval_first_iteration + 1
                  : 0,
          .current_objective_h = current.raw_cost_h(),
          .current_generalized_cost_h =
              generalized_cost(current, instance_, penalty_lambda_),
          .best_objective_h =
              best.has_value()
                  ? std::optional<double>(best->raw_cost_h())
                  : std::nullopt,
          .search_profile = config_.search_profile,
          .current_feasible = current.feasible(),
          .current_routes = current.plans.size(),
          .penalty_lambda = penalty_lambda_,
          .pool_size = pool_.main_size(),
          .pool_distinct_coverages = pool_.distinct_main_coverages(),
          .backup_pool_size = pool_.backup_size(),
          .pool_duals_valid = pool_.duals_valid(),
          .pool_effective_generation = pool_.effective_generation(),
          .pool_statistics = pool_.statistics(),
          .sp_calls = sp_calls_,
          .adaptive_sp_calls = adaptive_sp_calls_,
          .last_sp_iteration = last_sp_iteration_,
          .last_best_improvement_iteration =
              last_best_improvement_iteration_,
          .last_coverage_improvement_iteration =
              last_coverage_improvement_iteration_,
          .lp_calls = lp_calls_,
          .frvcp_calls = frvcp_calls_,
          .frvcp_solve_calls = frvcp_.solve_calls(),
          .frvcp_cache_hits = frvcp_.cache_hits(),
          .exact_plan_cache_size = factory_.exact_cache_size(),
          .exact_plan_cache_hits = factory_.exact_cache_hits(),
          .exact_plan_cache_misses = factory_.exact_cache_misses(),
          .exact_plan_cache_evictions = factory_.exact_cache_evictions(),
          .path_option_cache_size = path_sampler_.option_cache_size(),
          .path_option_cache_hits = path_sampler_.option_cache_hits(),
          .path_option_cache_misses = path_sampler_.option_cache_misses(),
          .path_option_cache_evictions =
              path_sampler_.option_cache_evictions(),
          .path_option_first_builds = path_sampler_.option_first_builds(),
          .path_option_rebuilds = path_sampler_.option_rebuilds(),
          .path_option_candidates_generated =
              path_sampler_.option_candidates_generated(),
          .path_option_nondominated_generated =
              path_sampler_.option_nondominated_generated(),
          .virtual_completion_cache_size =
              path_sampler_.virtual_completion_cache_size(),
          .virtual_completion_cache_hits =
              path_sampler_.virtual_completion_cache_hits(),
          .virtual_completion_cache_misses =
              path_sampler_.virtual_completion_cache_misses(),
          .virtual_completion_cache_evictions =
              path_sampler_.virtual_completion_cache_evictions(),
          .vnd_linked_completion_attempts =
              vnd_.statistics().linked_completion_attempts,
          .vnd_linked_completion_changed_attempts =
              vnd_.statistics().linked_completion_changed_attempts,
          .vnd_adaptive_customer_advances =
              vnd_.statistics().adaptive_customer_advances,
          .vnd_adaptive_replace_path_selections =
              vnd_.statistics().adaptive_replace_path_selections,
          .vnd_adaptive_inserted_replace_path_failures =
              vnd_.statistics().adaptive_inserted_replace_path_failures,
          .vnd_adaptive_final_replace_path_calls =
              vnd_.statistics().adaptive_final_replace_path_calls,
          .vnd_adaptive_final_replace_path_gaps_scanned =
              vnd_.statistics().adaptive_final_replace_path_gaps_scanned,
          .vnd_adaptive_final_replace_path_accepted =
              vnd_.statistics().adaptive_final_replace_path_accepted,
          .initial_insertion_statistics = initial_builder_.statistics(),
          .repair_insertion_statistics = perturbation_.statistics(),
          .elapsed_s = elapsed_seconds(started, checkpoint_time),
          .interval_elapsed_s =
              elapsed_seconds(interval_started, checkpoint_time),
          .initialization_s = initialization_runtime_s,
          .perturbation_s = interval_perturbation_s,
          .vnd_s = interval_vnd_s,
          .assembler_s = interval_assembler_s,
          .sp_model_s = sp_model_runtime_s_ - interval_sp_model_start_s,
          .lp_model_s = lp_model_runtime_s_ - interval_lp_model_start_s,
          .frvcp_model_s =
              frvcp_model_runtime_s_ - interval_frvcp_model_start_s,
          .evaluation_cache_size = evaluator_.cache_size(),
          .evaluation_cache_hits = cache.hits,
          .evaluation_cache_misses = cache.misses,
          .evaluation_cache_evictions = cache.evictions,
          .interval_evaluation_cache_hits =
              cache.hits - interval_cache_start.hits,
          .interval_evaluation_cache_misses =
              cache.misses - interval_cache_start.misses,
          .interval_evaluation_cache_evictions =
              cache.evictions - interval_cache_start.evictions});
      interval_started = SolverClock::now();
      interval_first_iteration = iterations + 1;
      interval_perturbation_s = 0.0;
      interval_vnd_s = 0.0;
      interval_assembler_s = 0.0;
      interval_cache_start = cache;
      interval_sp_model_start_s = sp_model_runtime_s_;
      interval_lp_model_start_s = lp_model_runtime_s_;
      interval_frvcp_model_start_s = frvcp_model_runtime_s_;
    }
    if (target_reached) {
      stop_reason = "rounded_target";
      break;
    }
    if (deadline_reached) {
      stop_reason = "time_limit";
      break;
    }
  }
  const double runtime = elapsed_seconds(started, SolverClock::now());
  return SolverResult{
      .best_solution = std::move(best),
      .iterations = iterations,
      .runtime_s = runtime,
      .final_penalty_lambda = penalty_lambda_,
      .pool_size = pool_.main_size(),
      .pool_distinct_coverages = pool_.distinct_main_coverages(),
      .backup_pool_size = pool_.backup_size(),
      .pool_duals_valid = pool_.duals_valid(),
      .pool_effective_generation = pool_.effective_generation(),
      .pool_statistics = pool_.statistics(),
      .sp_calls = sp_calls_,
      .adaptive_sp_calls = adaptive_sp_calls_,
      .lp_calls = lp_calls_,
      .frvcp_calls = frvcp_calls_,
      .frvcp_solve_calls = frvcp_.solve_calls(),
      .frvcp_cache_hits = frvcp_.cache_hits(),
      .exact_plan_cache_size = factory_.exact_cache_size(),
      .exact_plan_cache_hits = factory_.exact_cache_hits(),
      .exact_plan_cache_misses = factory_.exact_cache_misses(),
      .exact_plan_cache_evictions = factory_.exact_cache_evictions(),
      .path_option_cache_size = path_sampler_.option_cache_size(),
      .path_option_cache_hits = path_sampler_.option_cache_hits(),
      .path_option_cache_misses = path_sampler_.option_cache_misses(),
      .path_option_cache_evictions = path_sampler_.option_cache_evictions(),
      .path_option_first_builds = path_sampler_.option_first_builds(),
      .path_option_rebuilds = path_sampler_.option_rebuilds(),
      .path_option_candidates_generated =
          path_sampler_.option_candidates_generated(),
      .path_option_nondominated_generated =
          path_sampler_.option_nondominated_generated(),
      .virtual_completion_cache_size =
          path_sampler_.virtual_completion_cache_size(),
      .virtual_completion_cache_hits =
          path_sampler_.virtual_completion_cache_hits(),
      .virtual_completion_cache_misses =
          path_sampler_.virtual_completion_cache_misses(),
      .virtual_completion_cache_evictions =
          path_sampler_.virtual_completion_cache_evictions(),
      .evaluation_cache_hits = evaluator_.cache_hits(),
      .evaluation_cache_misses = evaluator_.cache_misses(),
      .evaluation_cache_evictions = evaluator_.cache_evictions(),
      .initialization_runtime_s = initialization_runtime_s,
      .perturbation_runtime_s = perturbation_runtime_s,
      .vnd_runtime_s = vnd_runtime_s,
      .assembler_runtime_s = assembler_runtime_s,
      .sp_model_runtime_s = sp_model_runtime_s_,
      .lp_model_runtime_s = lp_model_runtime_s_,
      .frvcp_model_runtime_s = frvcp_model_runtime_s_,
      .initial_insertion_statistics = initial_builder_.statistics(),
      .repair_insertion_statistics = perturbation_.statistics(),
      .vnd_statistics = vnd_.statistics(),
      .stop_reason = std::move(stop_reason)};
}

}  // namespace ils_sp
