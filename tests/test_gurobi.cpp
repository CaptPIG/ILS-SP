#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "ils_sp/gurobi_models.hpp"
#include "ils_sp/solver.hpp"

namespace {

void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

void require_close(double actual, double expected, double tolerance,
                   const std::string& message) {
  if (std::abs(actual - expected) > tolerance) {
    throw std::runtime_error(message + ": actual=" + std::to_string(actual) +
                             " expected=" + std::to_string(expected));
  }
}

ils_sp::Instance make_partition_instance() {
  return ils_sp::Instance{
      "tiny_partition",
      "synthetic",
      {{.id = 0, .kind = ils_sp::NodeKind::Depot, .x_km = 0.0, .y_km = 0.0},
       {.id = 1,
        .kind = ils_sp::NodeKind::Customer,
        .x_km = 1.0,
        .y_km = 0.0},
       {.id = 2,
        .kind = ils_sp::NodeKind::Customer,
        .x_km = 2.0,
        .y_km = 0.0}},
      {.speed_kmph = 10.0,
       .consumption_wh_per_km = 1.0,
       .battery_capacity_wh = 100.0,
       .max_route_duration_h = 10.0},
      {}};
}

ils_sp::Instance make_charging_instance() {
  std::unordered_map<std::string, ils_sp::ChargingCurve> curves;
  curves.emplace(
      "linear",
      ils_sp::ChargingCurve{"linear", {{.time_h = 0.0, .energy_wh = 0.0},
                                        {.time_h = 1.0,
                                         .energy_wh = 10.0}}});
  return ils_sp::Instance{
      "tiny_charging",
      "synthetic",
      {{.id = 0, .kind = ils_sp::NodeKind::Depot, .x_km = 0.0, .y_km = 0.0},
       {.id = 1,
        .kind = ils_sp::NodeKind::Customer,
        .x_km = 10.0,
        .y_km = 0.0},
       {.id = 2,
        .kind = ils_sp::NodeKind::Station,
        .x_km = 5.0,
        .y_km = 0.0,
        .station_type = "linear"}},
      {.speed_kmph = 10.0,
       .consumption_wh_per_km = 1.0,
       .battery_capacity_wh = 10.0,
       .max_route_duration_h = 10.0},
      std::move(curves)};
}

ils_sp::Instance make_pruning_instance() {
  std::unordered_map<std::string, ils_sp::ChargingCurve> curves;
  curves.emplace(
      "fast",
      ils_sp::ChargingCurve{"fast", {{.time_h = 0.0, .energy_wh = 0.0},
                                      {.time_h = 1.0,
                                       .energy_wh = 10.0}}});
  curves.emplace(
      "slow",
      ils_sp::ChargingCurve{"slow", {{.time_h = 0.0, .energy_wh = 0.0},
                                      {.time_h = 2.0,
                                       .energy_wh = 10.0}}});
  return ils_sp::Instance{
      "tiny_pruning",
      "synthetic",
      {{.id = 0, .kind = ils_sp::NodeKind::Depot, .x_km = 0.0, .y_km = 0.0},
       {.id = 1,
        .kind = ils_sp::NodeKind::Customer,
        .x_km = 10.0,
        .y_km = 0.0},
       {.id = 2,
        .kind = ils_sp::NodeKind::Station,
        .x_km = 5.0,
        .y_km = 0.0,
        .station_type = "slow"},
       {.id = 3,
        .kind = ils_sp::NodeKind::Station,
        .x_km = 5.0,
        .y_km = 0.0,
        .station_type = "fast"},
       {.id = 4,
        .kind = ils_sp::NodeKind::Station,
        .x_km = 100.0,
        .y_km = 0.0,
        .station_type = "fast"}},
      {.speed_kmph = 10.0,
       .consumption_wh_per_km = 1.0,
       .battery_capacity_wh = 10.0,
       .max_route_duration_h = 10.0},
      std::move(curves)};
}

}  // namespace

int main() {
  try {
    GRBEnv environment(true);
    environment.set(GRB_IntParam_OutputFlag, 0);
    environment.start();

    const ils_sp::Instance partition_instance = make_partition_instance();
    ils_sp::RouteEvaluator partition_evaluator(partition_instance, 32);
    ils_sp::PlanFactory partition_factory(partition_instance,
                                          partition_evaluator, 32);
    const ils_sp::Plan one =
        partition_factory.make_plan({1}, ils_sp::direct_gaps(1));
    const ils_sp::Plan two =
        partition_factory.make_plan({2}, ils_sp::direct_gaps(1));
    const ils_sp::Plan combined =
        partition_factory.make_plan({1, 2}, ils_sp::direct_gaps(2));
    std::vector<ils_sp::RouteColumn> columns{
        ils_sp::RouteColumn::from_plan(one),
        ils_sp::RouteColumn::from_plan(two),
        ils_sp::RouteColumn::from_plan(combined)};

    require_close(ils_sp::GurobiConfig{}.mip_gap, 1e-4, 1e-12,
                  "paper-default MIP gap is not 1e-4");
    require(ils_sp::SolverConfig{}.sp_incumbent_stall_limit_s ==
                    std::optional<double>(30.0) &&
                ils_sp::SolverConfig{}.sp_maximum_time_s ==
                    std::optional<double>(90.0) &&
                !ils_sp::SolverConfig{}
                     .adaptive_sp_stall_iterations.has_value() &&
                !ils_sp::SolverConfig{}
                     .rounded_target_objective_h.has_value(),
            "SP engineering stop defaults or opt-in target are wrong");

    ils_sp::GurobiSetPartitioning partitioning(partition_instance,
                                                environment);
    std::vector<ils_sp::SetPartitioningProgress> partition_progress;
    const ils_sp::SetPartitioningResult integer =
        partitioning.solve_integer(
            columns,
            ils_sp::SetPartitioningSolveOptions{
                .mip_start_indices = {0, 1},
                .progress = [&](const ils_sp::SetPartitioningProgress& value) {
                  partition_progress.push_back(value);
                }});
    require(integer.status == GRB_OPTIMAL,
            "integer set partitioning did not solve optimally");
    require(integer.selected_columns.size() == 1 &&
                integer.selected_columns.front().customer_set ==
                    std::vector<int>({1, 2}),
            "integer set partitioning did not select the combined route");
    require_close(*integer.objective_h, combined.evaluation->raw_cost_h, 1e-9,
                  "integer set partitioning objective is wrong");
    require(!integer.applied_time_limit_s.has_value() &&
                integer.model_size.variables == columns.size() &&
                integer.model_size.binary_variables == columns.size() &&
                integer.model_size.linear_constraints == 2 &&
                integer.model_size.linear_nonzeros == 4 &&
                integer.model_build_runtime_s >= 0.0 &&
                integer.mip_start_routes == 2 &&
                integer.mip_start_objective_h.has_value(),
            "paper-default SP unexpectedly has a limit or wrong model size");
    require(!partition_progress.empty() &&
                partition_progress.front().event ==
                    ils_sp::SetPartitioningProgressEvent::ModelBuilt &&
                partition_progress.front().model_size.linear_nonzeros == 4,
            "SP model-build checkpoint is missing or incomplete");

    const double warm_start_target_h = one.evaluation->raw_cost_h +
                                       two.evaluation->raw_cost_h;
    const ils_sp::SetPartitioningResult target_stopped =
        partitioning.solve_integer(
            columns,
            ils_sp::SetPartitioningSolveOptions{
                .mip_start_indices = {0, 1},
                .rounded_target_objective_h = warm_start_target_h});
    require((target_stopped.status == GRB_USER_OBJ_LIMIT ||
             target_stopped.status == GRB_OPTIMAL) &&
                target_stopped.objective_h.has_value() &&
                target_stopped.rounded_target_reached &&
                std::llround(*target_stopped.objective_h * 100.0) <=
                    std::llround(warm_start_target_h * 100.0),
            "explicit rounded target did not stop with a qualifying SP "
            "incumbent: status=" + target_stopped.status_name +
                " objective=" +
                (target_stopped.objective_h.has_value()
                     ? std::to_string(*target_stopped.objective_h)
                     : std::string("none")) +
                " reached=" +
                std::to_string(target_stopped.rounded_target_reached));

    const ils_sp::LinearRelaxationResult relaxation =
        partitioning.solve_relaxation(columns);
    require(relaxation.status == GRB_OPTIMAL && relaxation.duals.size() == 2,
            "LP relaxation did not return both customer duals");
    require_close(*relaxation.objective_h, combined.evaluation->raw_cost_h,
                  1e-9, "LP relaxation objective is wrong");
    require(!relaxation.applied_time_limit_s.has_value() &&
                relaxation.model_size.variables == columns.size() &&
                relaxation.model_size.binary_variables == 0 &&
                relaxation.model_size.linear_constraints == 2,
            "paper-default LP unexpectedly has a limit or wrong model size");

    ils_sp::RouteEvaluator direct_evaluator(partition_instance, 32);
    ils_sp::PlanFactory direct_factory(partition_instance, direct_evaluator,
                                       32);
    ils_sp::GurobiFrvcp direct_frvcp(partition_instance, direct_evaluator,
                                    direct_factory, environment);
    const ils_sp::FrvcpResult direct_optimum =
        direct_frvcp.optimize(std::vector<int>{1, 2});
    require(direct_optimum.status == GRB_OPTIMAL &&
                direct_optimum.plan.has_value() &&
                direct_optimum.plan->exact_charging &&
                direct_optimum.analytic_optimum &&
                direct_frvcp.solve_calls() == 0 &&
                direct_optimum.model_size.variables == 0 &&
                direct_optimum.objective_lower_bound_h.has_value(),
            "feasible all-direct route was not certified without Gurobi");
    require_close(*direct_optimum.objective_lower_bound_h,
                  *direct_optimum.model_objective_h, 1e-9,
                  "direct-route certificate and objective disagree");

    ils_sp::RouteEvaluator label_evaluator(partition_instance, 8);
    ils_sp::PlanFactory label_factory(partition_instance, label_evaluator, 1);
    ils_sp::GurobiFrvcp label_frvcp(partition_instance, label_evaluator,
                                   label_factory, environment);
    const ils_sp::FrvcpResult first_label =
        label_frvcp.optimize(std::vector<int>{1, 2});
    const ils_sp::FrvcpResult second_label =
        label_frvcp.optimize(std::vector<int>{2, 1});
    require(first_label.plan.has_value() && second_label.plan.has_value() &&
                !label_factory.exact_plan({1, 2}).has_value(),
            "exact-label eviction fixture did not evict the first route");
    const std::size_t solves_before_label_reuse = label_frvcp.solve_calls();
    const ils_sp::FrvcpResult restored_label =
        label_frvcp.optimize(*first_label.plan);
    require(restored_label.status == GRB_OPTIMAL && restored_label.cache_hit &&
                !restored_label.analytic_optimum &&
                label_frvcp.solve_calls() == solves_before_label_reuse &&
                label_factory.exact_plan({1, 2}).has_value(),
            "optimized FRVCP label did not survive exact-cache eviction");

    const ils_sp::Instance charging_instance = make_charging_instance();
    ils_sp::RouteEvaluator charging_evaluator(charging_instance, 32);
    ils_sp::PlanFactory charging_factory(charging_instance, charging_evaluator,
                                         32);
    ils_sp::GurobiFrvcp frvcp(charging_instance, charging_evaluator,
                              charging_factory, environment);
    const ils_sp::Plan charging_incumbent =
        charging_factory.make_plan({1}, {{2}, {2}});
    const ils_sp::FrvcpResult optimized =
        frvcp.optimize(charging_incumbent);
    require(optimized.status == GRB_OPTIMAL && optimized.plan.has_value(),
            "FRVCP did not find the feasible charging route");
    require(optimized.plan->route.visits == std::vector<int>({2, 1, 2}),
            "FRVCP did not use the required station on both customer gaps");
    require_close(*optimized.model_objective_h, 3.0, 1e-6,
                  "FRVCP travel-plus-charging objective is wrong");
    require(optimized.plan->evaluation->feasible &&
                optimized.plan->exact_charging,
            "FRVCP result was not replayed as an exact feasible plan");
    require(!optimized.applied_time_limit_s.has_value() &&
                optimized.model_size.variables > 0 &&
                optimized.model_size.binary_variables > 0 &&
                optimized.model_size.general_constraints > 0 &&
                optimized.mip_start_supplied &&
                optimized.mip_start_objective_h.has_value() &&
                optimized.objective_lower_bound_h.has_value() &&
                optimized.network_size.customer_count == 1 &&
                optimized.network_size.gap_count == 2,
            "paper-default FRVCP limit or model/network statistics are wrong");
    require_close(*optimized.mip_start_objective_h, 3.0, 1e-6,
                  "FRVCP did not preserve the selected-route MIP start");
    require_close(*optimized.objective_lower_bound_h, 3.0, 1e-6,
                  "FRVCP route lower bound is wrong");

    const ils_sp::FrvcpResult cached = frvcp.optimize({1});
    require(cached.cache_hit && frvcp.solve_calls() == 1 &&
                frvcp.cache_hits() == 1,
            "exact FRVCP cache did not suppress the repeated Gurobi solve");

    ils_sp::RouteEvaluator limited_evaluator(charging_instance, 32);
    ils_sp::PlanFactory limited_factory(charging_instance, limited_evaluator,
                                        32);
    ils_sp::GurobiConfig limited_config;
    limited_config.time_limit_s = 1e-9;
    ils_sp::GurobiFrvcp limited_frvcp(
        charging_instance, limited_evaluator, limited_factory, environment,
        limited_config);
    const ils_sp::Plan limited_incumbent =
        limited_factory.make_plan({1}, {{2}, {2}});
    const ils_sp::FrvcpResult limited =
        limited_frvcp.optimize(limited_incumbent);
    require(limited.plan.has_value() && limited.plan->evaluation->feasible &&
                limited.applied_time_limit_s.has_value() &&
                limited.mip_start_supplied,
            "time-limited FRVCP discarded its feasible incumbent");
    if (limited.status != GRB_OPTIMAL) {
      require(!limited.plan->exact_charging,
              "non-optimal FRVCP incumbent entered the exact cache contract");
    }

    const ils_sp::Instance pruning_instance = make_pruning_instance();
    ils_sp::RouteEvaluator pruning_evaluator(pruning_instance, 32);
    ils_sp::PlanFactory pruning_factory(pruning_instance, pruning_evaluator,
                                        32);
    ils_sp::GurobiFrvcp pruning_frvcp(pruning_instance, pruning_evaluator,
                                     pruning_factory, environment);
    const ils_sp::FrvcpResult pruned = pruning_frvcp.optimize({1});
    require(pruned.status == GRB_OPTIMAL && pruned.plan.has_value() &&
                pruned.plan->route.visits == std::vector<int>({3, 1, 3}),
            "safe FRVCP pruning removed the faster feasible station path");
    require_close(*pruned.model_objective_h, 3.0, 1e-6,
                  "pruned FRVCP objective is wrong");
    require(pruned.network_size.enumerated_paths >
                    pruned.network_size.active_paths &&
                pruned.network_size.active_paths == 4 &&
                pruned.network_size.active_station_visits == 2 &&
                pruned.network_size.candidate_station_copies ==
                    pruned.network_size.enumerated_station_visits &&
                pruned.network_size.active_station_copies ==
                    pruned.network_size.active_station_visits &&
                pruned.network_size.candidate_arcs == 0 &&
                pruned.network_size.active_arcs == 0,
            "path FRVCP did not remove dominated or unreachable CSPs");

    ils_sp::SolverConfig checkpoint_config;
    checkpoint_config.paper.max_iterations = 2;
    checkpoint_config.paper.sp_frequency = 100;
    checkpoint_config.paper.backup_frequency = 100;
    checkpoint_config.progress_interval = 2;
    std::vector<ils_sp::SolverProgress> checkpoints;
    ils_sp::IlsSpSolver checkpoint_solver(partition_instance, environment,
                                          checkpoint_config);
    const ils_sp::SolverResult checkpoint_result = checkpoint_solver.solve(
        [&](const ils_sp::SolverProgress& checkpoint) {
          checkpoints.push_back(checkpoint);
        });
    require(checkpoint_result.iterations == 2 && checkpoints.size() == 2,
            "solver did not emit initialization and final checkpoints");
    const ils_sp::SolverProgress& initialization = checkpoints.front();
    require(initialization.reason ==
                ils_sp::SolverCheckpointReason::Initialization &&
                initialization.iteration == 0 &&
                initialization.interval_iterations == 0,
            "initialization checkpoint boundary is wrong");
    require(initialization.current_objective_h.has_value() &&
                initialization.current_generalized_cost_h.has_value() &&
                initialization.current_routes > 0 &&
                initialization.initialization_s >= 0.0,
            "initialization checkpoint omitted solution or timing state");
    require(initialization.interval_evaluation_cache_hits ==
                initialization.evaluation_cache_hits &&
                initialization.interval_evaluation_cache_misses ==
                    initialization.evaluation_cache_misses &&
                initialization.interval_evaluation_cache_evictions ==
                    initialization.evaluation_cache_evictions,
            "initialization checkpoint cache delta is inconsistent");

    const ils_sp::SolverProgress& final_checkpoint = checkpoints.back();
    require(final_checkpoint.reason ==
                ils_sp::SolverCheckpointReason::MaxIterations &&
                final_checkpoint.iteration == 2 &&
                final_checkpoint.interval_first_iteration == 1 &&
                final_checkpoint.interval_iterations == 2,
            "final iteration checkpoint boundary is wrong");
    require(final_checkpoint.perturbation_s >= 0.0 &&
                final_checkpoint.vnd_s >= 0.0 &&
                final_checkpoint.assembler_s >= 0.0 &&
                final_checkpoint.interval_elapsed_s + 1e-9 >=
                    final_checkpoint.perturbation_s +
                        final_checkpoint.vnd_s +
                        final_checkpoint.assembler_s,
            "checkpoint phase timings are inconsistent");
    require(initialization.interval_evaluation_cache_hits +
                    final_checkpoint.interval_evaluation_cache_hits ==
                checkpoint_result.evaluation_cache_hits &&
                initialization.interval_evaluation_cache_misses +
                        final_checkpoint.interval_evaluation_cache_misses ==
                    checkpoint_result.evaluation_cache_misses &&
                initialization.interval_evaluation_cache_evictions +
                        final_checkpoint.interval_evaluation_cache_evictions ==
                    checkpoint_result.evaluation_cache_evictions,
            "checkpoint cache deltas do not reconstruct final totals");
    require_close(checkpoint_result.initialization_runtime_s,
                  initialization.initialization_s, 1e-9,
                  "final result lost initialization timing");

    ils_sp::SolverConfig adaptive_config;
    adaptive_config.paper.max_iterations = 4;
    adaptive_config.paper.sp_frequency = 100;
    adaptive_config.paper.backup_frequency = 100;
    adaptive_config.paper.pool_capacity = 2;
    adaptive_config.adaptive_sp_stall_iterations = 1;
    adaptive_config.adaptive_sp_minimum_interval = 1;
    adaptive_config.progress_interval = 1;
    std::vector<ils_sp::ExactStageProgress> adaptive_exact_checkpoints;
    ils_sp::IlsSpSolver adaptive_solver(partition_instance, environment,
                                        adaptive_config);
    const ils_sp::SolverResult adaptive_result = adaptive_solver.solve(
        [](const ils_sp::SolverProgress&) {},
        [&](const ils_sp::ExactStageProgress& checkpoint) {
          adaptive_exact_checkpoints.push_back(checkpoint);
        });
    const auto adaptive_start = std::find_if(
        adaptive_exact_checkpoints.begin(), adaptive_exact_checkpoints.end(),
        [](const ils_sp::ExactStageProgress& checkpoint) {
          return checkpoint.stage == ils_sp::ExactStage::SetPartitioning &&
                 checkpoint.boundary ==
                     ils_sp::ExactStageBoundary::Start &&
                 checkpoint.sp_trigger_source == "pool_stall";
        });
    require(adaptive_result.adaptive_sp_calls > 0 &&
                adaptive_start != adaptive_exact_checkpoints.end() &&
                adaptive_start->iteration <
                    adaptive_config.paper.sp_frequency,
            "pool-stall trigger did not call SP before the regular frequency");

    ils_sp::SolverConfig additive_config = adaptive_config;
    additive_config.paper.max_iterations = 4;
    additive_config.paper.sp_frequency = 2;
    additive_config.paper.backup_frequency = 2;
    std::vector<ils_sp::ExactStageProgress> additive_exact_checkpoints;
    ils_sp::IlsSpSolver additive_solver(partition_instance, environment,
                                        additive_config);
    const ils_sp::SolverResult additive_result = additive_solver.solve(
        [](const ils_sp::SolverProgress&) {},
        [&](const ils_sp::ExactStageProgress& checkpoint) {
          additive_exact_checkpoints.push_back(checkpoint);
        });
    std::vector<std::size_t> periodic_starts;
    for (const ils_sp::ExactStageProgress& checkpoint :
         additive_exact_checkpoints) {
      if (checkpoint.stage == ils_sp::ExactStage::SetPartitioning &&
          checkpoint.boundary == ils_sp::ExactStageBoundary::Start &&
          checkpoint.sp_trigger_source == "frequency") {
        periodic_starts.push_back(checkpoint.iteration);
      }
    }
    require(additive_result.iterations == 4 &&
                periodic_starts == std::vector<std::size_t>({2, 4}),
            "an adaptive SP shifted the absolute paper SP cadence");

    ils_sp::SolverConfig exact_config;
    exact_config.paper.max_iterations = 2;
    exact_config.paper.sp_frequency = 1;
    exact_config.paper.backup_frequency = 1;
    exact_config.progress_interval = 1;
    exact_config.global_time_limit_s = 5.0;
    std::vector<ils_sp::ExactStageProgress> exact_checkpoints;
    ils_sp::IlsSpSolver exact_solver(partition_instance, environment,
                                     exact_config);
    const ils_sp::SolverResult exact_result = exact_solver.solve(
        [](const ils_sp::SolverProgress&) {},
        [&](const ils_sp::ExactStageProgress& checkpoint) {
          exact_checkpoints.push_back(checkpoint);
        });
    require(exact_result.iterations == 2 && exact_result.sp_calls == 2 &&
                exact_result.lp_calls == 2 && exact_result.frvcp_calls > 0,
            "two-iteration solver did not traverse SP, LP and FRVCP");
    const auto has_boundary = [&](ils_sp::ExactStage stage,
                                  ils_sp::ExactStageBoundary boundary) {
      return std::any_of(exact_checkpoints.begin(), exact_checkpoints.end(),
                         [&](const ils_sp::ExactStageProgress& checkpoint) {
                           return checkpoint.stage == stage &&
                                  checkpoint.boundary == boundary;
                         });
    };
    for (const ils_sp::ExactStage stage :
         {ils_sp::ExactStage::SetPartitioning,
          ils_sp::ExactStage::LinearRelaxation,
          ils_sp::ExactStage::Frvcp}) {
      require(has_boundary(stage, ils_sp::ExactStageBoundary::Start) &&
                  has_boundary(stage, ils_sp::ExactStageBoundary::End),
              "exact-stage start/end checkpoints are incomplete");
    }
    require(has_boundary(ils_sp::ExactStage::SetPartitioning,
                         ils_sp::ExactStageBoundary::ModelBuilt),
            "SP model-built checkpoint is missing");
    std::vector<ils_sp::ExactStageProgress> sp_starts;
    std::copy_if(exact_checkpoints.begin(), exact_checkpoints.end(),
                 std::back_inserter(sp_starts),
                 [](const ils_sp::ExactStageProgress& checkpoint) {
                   return checkpoint.stage ==
                              ils_sp::ExactStage::SetPartitioning &&
                          checkpoint.boundary ==
                              ils_sp::ExactStageBoundary::Start;
                 });
    require(sp_starts.size() == 2 &&
                sp_starts[0].mip_start_source == "best_feasible" &&
                sp_starts[1].mip_start_source == "previous_sp" &&
                sp_starts[0].mip_start_routes > 0 &&
                sp_starts[1].mip_start_routes > 0 &&
                sp_starts[0].pool_distinct_coverages > 0 &&
                sp_starts[0].pool_distinct_coverages <=
                    sp_starts[0].pool_columns &&
                sp_starts[1].pool_distinct_coverages > 0 &&
                sp_starts[1].pool_distinct_coverages <=
                    sp_starts[1].pool_columns,
            "SP did not use best-feasible then previous-SP MIP starts");
    for (const ils_sp::ExactStageProgress& checkpoint : exact_checkpoints) {
      require(checkpoint.remaining_time_s.has_value() &&
                  *checkpoint.remaining_time_s >= 0.0 &&
                  *checkpoint.remaining_time_s <= 5.0,
              "global remaining time is missing from exact checkpoint");
      if (checkpoint.boundary != ils_sp::ExactStageBoundary::End ||
          (checkpoint.stage == ils_sp::ExactStage::Frvcp &&
           (checkpoint.cache_hit || checkpoint.analytic_optimum ||
            checkpoint.lower_bound_proved_infeasible))) {
        continue;
      }
      require(checkpoint.applied_time_limit_s.has_value() &&
                  *checkpoint.applied_time_limit_s > 0.0 &&
                  *checkpoint.applied_time_limit_s <= 5.0 &&
                  checkpoint.status.has_value() &&
                  checkpoint.model_size.variables > 0,
              "global deadline was not applied to an exact Gurobi model");
    }
    require(std::any_of(
                exact_checkpoints.begin(), exact_checkpoints.end(),
                [](const ils_sp::ExactStageProgress& checkpoint) {
                  return checkpoint.stage == ils_sp::ExactStage::Frvcp &&
                         checkpoint.boundary ==
                             ils_sp::ExactStageBoundary::End &&
                         checkpoint.analytic_optimum &&
                         checkpoint.objective_lower_bound_h.has_value() &&
                         checkpoint.model_size.variables == 0;
                }),
            "FRVCP analytic certificate is missing from checkpoints");

    ils_sp::SolverConfig expired_config;
    expired_config.paper.max_iterations = 2;
    expired_config.paper.sp_frequency = 1;
    expired_config.paper.backup_frequency = 1;
    expired_config.global_time_limit_s = 1e-12;
    ils_sp::IlsSpSolver expired_solver(partition_instance, environment,
                                       expired_config);
    const ils_sp::SolverResult expired_result = expired_solver.solve();
    require(expired_result.stop_reason == "time_limit" &&
                expired_result.iterations == 0 &&
                expired_result.pool_size == 0 &&
                expired_result.sp_calls == 0 && expired_result.lp_calls == 0 &&
                expired_result.frvcp_calls == 0,
            "pre-assembler timeout was reported as a completed ILS iteration");

    std::cout << "all Gurobi integration tests passed\n";
    return 0;
  } catch (const GRBException& error) {
    std::cerr << "Gurobi test failure " << error.getErrorCode() << ": "
              << error.getMessage() << '\n';
    return 1;
  } catch (const std::exception& error) {
    std::cerr << "test failure: " << error.what() << '\n';
    return 1;
  }
}
