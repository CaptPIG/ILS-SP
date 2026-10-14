#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

#include "ils_sp/solver.hpp"
#include "ils_sp/detail/input_contract.hpp"

#ifndef ILS_SP_TYPED_DELETE_ENTRY
#define ILS_SP_TYPED_DELETE_ENTRY 0
#endif

#ifndef ILS_SP_FORWARD_ONLY_TWO_OPT
#define ILS_SP_FORWARD_ONLY_TWO_OPT 0
#endif

#ifndef ILS_SP_PHYSICAL_REVERSE_EMPTY_CHILD
#define ILS_SP_PHYSICAL_REVERSE_EMPTY_CHILD 0
#endif
#include "ils_sp/xml.hpp"

namespace {

struct Arguments {
  std::filesystem::path instance;
  ils_sp::SolverConfig config;
};

[[nodiscard]] std::string_view checkpoint_reason_name(
    ils_sp::SolverCheckpointReason reason) {
  switch (reason) {
    case ils_sp::SolverCheckpointReason::Initialization:
      return "initialization";
    case ils_sp::SolverCheckpointReason::ProgressInterval:
      return "progress_interval";
    case ils_sp::SolverCheckpointReason::SetPartitioning:
      return "set_partitioning";
    case ils_sp::SolverCheckpointReason::RoundedTarget:
      return "rounded_target";
    case ils_sp::SolverCheckpointReason::TimeLimit:
      return "time_limit";
    case ils_sp::SolverCheckpointReason::MaxIterations:
      return "max_iterations";
  }
  throw std::logic_error("unknown solver checkpoint reason");
}

[[nodiscard]] std::string_view exact_stage_name(ils_sp::ExactStage stage) {
  switch (stage) {
    case ils_sp::ExactStage::SetPartitioning:
      return "sp";
    case ils_sp::ExactStage::LinearRelaxation:
      return "lp";
    case ils_sp::ExactStage::Frvcp:
      return "frvcp";
  }
  throw std::logic_error("unknown exact stage");
}

[[nodiscard]] std::string_view exact_boundary_name(
    ils_sp::ExactStageBoundary boundary) {
  switch (boundary) {
    case ils_sp::ExactStageBoundary::Start:
      return "start";
    case ils_sp::ExactStageBoundary::ModelBuilt:
      return "model_built";
    case ils_sp::ExactStageBoundary::Progress:
      return "progress";
    case ils_sp::ExactStageBoundary::End:
      return "end";
  }
  throw std::logic_error("unknown exact-stage boundary");
}

[[nodiscard]] std::string_view search_profile_name(
    ils_sp::SearchProfile profile) {
  switch (profile) {
    case ils_sp::SearchProfile::Paper:
      return "paper";
    case ils_sp::SearchProfile::LinkedWinner:
      return "linked-winner";
    case ils_sp::SearchProfile::LinkedInsertionCandidates:
      return "linked-insertion-candidates";
    case ils_sp::SearchProfile::GranularAdaptive:
      return "granular-adaptive";
  }
  throw std::logic_error("unknown search profile");
}

void log_optional(std::string_view name, std::optional<double> value,
                  std::string_view missing = "none") {
  std::cerr << ' ' << name << '=';
  if (value.has_value()) {
    std::cerr << *value;
  } else {
    std::cerr << missing;
  }
}

void log_exact_stage(const ils_sp::ExactStageProgress& progress) {
  std::cerr << std::fixed << std::setprecision(6)
            << "exact_stage=" << exact_stage_name(progress.stage)
            << " event=" << exact_boundary_name(progress.boundary)
            << " iteration=" << progress.iteration
            << " pool_columns=" << progress.pool_columns;
  if (progress.stage == ils_sp::ExactStage::SetPartitioning) {
    std::cerr << " pool_distinct_coverages="
              << progress.pool_distinct_coverages
              << " pool_effective_generation="
              << progress.pool_effective_generation;
  }
  if (progress.stage == ils_sp::ExactStage::Frvcp) {
    std::cerr << " route=" << progress.route_index << '/'
              << progress.selected_routes
              << " route_customers=" << progress.route_customer_count
              << " route_gaps=" << progress.route_gap_count;
    log_optional("assembled", progress.assembled_objective_h);
  } else if (progress.selected_routes != 0) {
    std::cerr << " selected_routes=" << progress.selected_routes;
  }
  if (progress.stage == ils_sp::ExactStage::SetPartitioning) {
    std::cerr << " sp_trigger="
              << (progress.sp_trigger_source.empty()
                      ? "frequency"
                      : progress.sp_trigger_source)
              << " mip_start_source="
              << (progress.mip_start_source.empty()
                      ? "none"
                      : progress.mip_start_source)
              << " mip_start_routes=" << progress.mip_start_routes;
    log_optional("mip_start_objective", progress.mip_start_objective_h);
  } else if (progress.stage == ils_sp::ExactStage::Frvcp &&
             progress.boundary == ils_sp::ExactStageBoundary::Start) {
    std::cerr << " mip_start_source="
              << (progress.mip_start_source.empty()
                      ? "none"
                      : progress.mip_start_source);
    log_optional("mip_start_objective", progress.mip_start_objective_h);
    std::cerr << " physical_visits=";
    for (std::size_t i = 0; i < progress.route_visits.size(); ++i) {
      if (i != 0) std::cerr << ',';
      std::cerr << progress.route_visits[i];
    }
  }
  log_optional("remaining_s", progress.remaining_time_s, "unlimited");
  if (!progress.backend_phase.empty()) {
    std::cerr << " backend_phase=" << progress.backend_phase
              << " wall_s=" << progress.wall_runtime_s
              << " mip_start_supplied=" << progress.mip_start_supplied;
  }
  if (progress.boundary == ils_sp::ExactStageBoundary::ModelBuilt) {
    std::cerr << " build_s=" << progress.model_build_runtime_s
              << " variables=" << progress.model_size.variables
              << " binaries=" << progress.model_size.binary_variables
              << " sos_constraints=" << progress.model_size.sos_constraints
              << " linear_constraints="
              << progress.model_size.linear_constraints
              << " linear_nonzeros="
              << progress.model_size.linear_nonzeros;
    if (progress.stage == ils_sp::ExactStage::Frvcp) {
      std::cerr << " path_generation_s=" << progress.path_generation_runtime_s
                << " enumerated_paths=" << progress.network_size.enumerated_paths
                << " active_paths=" << progress.network_size.active_paths
                << " active_station_visits="
                << progress.network_size.active_station_visits;
    }
  } else if (progress.boundary == ils_sp::ExactStageBoundary::Progress &&
             progress.backend_phase.empty()) {
    log_optional("incumbent", progress.objective_h);
    log_optional("bound", progress.best_bound_h);
    log_optional("gap", progress.relative_gap);
    log_optional("nodes", progress.explored_nodes);
    std::cerr << " model_s=" << progress.model_runtime_s
              << " solutions=" << progress.solution_count
              << " rounded_target_reached="
              << progress.rounded_target_reached
              << " incumbent_stall_reached="
              << progress.incumbent_stall_limit_reached;
  } else if (progress.boundary == ils_sp::ExactStageBoundary::End) {
    std::cerr << " status="
              << (progress.status_name.empty() ? "none"
                                               : progress.status_name);
    log_optional("objective", progress.objective_h);
    std::cerr << " rc_fixed_columns=" << progress.rc_fixed_columns
              << " rc_considered_columns=" << progress.rc_considered_columns
              << " lp_prepass_s=" << progress.lp_prepass_wall_s
              << " lp_reused_prepass=" << progress.lp_reused_prepass;
    log_optional("certificate_gap_upper_h", progress.certificate_gap_upper_h);
    log_optional("certificate_acceptance_gap_h", progress.certificate_acceptance_gap_h);
    if (progress.stage == ils_sp::ExactStage::SetPartitioning) {
      log_optional("bound", progress.best_bound_h);
      log_optional("gap", progress.relative_gap);
      log_optional("nodes", progress.explored_nodes);
      std::cerr << " solutions=" << progress.solution_count
                << " rounded_target_reached="
                << progress.rounded_target_reached
                << " incumbent_stall_reached="
                << progress.incumbent_stall_limit_reached;
    }
    log_optional("model_limit_s", progress.applied_time_limit_s,
                 "unlimited");
    std::cerr << " build_s=" << progress.model_build_runtime_s
              << " model_s=" << progress.model_runtime_s
              << " wall_s=" << progress.wall_runtime_s
              << " variables=" << progress.model_size.variables
              << " binaries=" << progress.model_size.binary_variables
              << " linear_constraints="
              << progress.model_size.linear_constraints
              << " sos_constraints="
              << progress.model_size.sos_constraints
              << " linear_nonzeros="
              << progress.model_size.linear_nonzeros;
    if (progress.stage == ils_sp::ExactStage::Frvcp) {
      std::cerr << " cache_hit=" << progress.cache_hit
                << " mip_start_supplied=" << progress.mip_start_supplied;
      log_optional("mip_start_objective",
                   progress.mip_start_objective_h);
      log_optional("objective_lower_bound",
                   progress.objective_lower_bound_h);
      std::cerr << " analytic_optimum=" << progress.analytic_optimum
                << " fixed_order_optimal=" << progress.fixed_order_optimal
                << " lower_bound_infeasible="
                << progress.lower_bound_proved_infeasible
                << " path_generation_s="
                << progress.path_generation_runtime_s
                << " enumerated_paths="
                << progress.network_size.enumerated_paths
                << " active_paths=" << progress.network_size.active_paths
                << " enumerated_station_visits="
                << progress.network_size.enumerated_station_visits
                << " active_station_visits="
                << progress.network_size.active_station_visits
                << " candidate_station_copies="
                << progress.network_size.candidate_station_copies
                << " active_station_copies="
                << progress.network_size.active_station_copies
                << " candidate_arcs="
                << progress.network_size.candidate_arcs
                << " active_arcs=" << progress.network_size.active_arcs
                << " rounded_target_reached="
                << progress.rounded_target_reached;
    }
  }
  std::cerr << '\n';
}

[[noreturn]] void usage(int status) {
  std::ostream& output = status == 0 ? std::cout : std::cerr;
  output
      << "Usage: ils-sp --instance FILE [options]\n"
      << "  --iterations N          ILS iterations (paper: 50000)\n"
      << "  --seed N                random seed\n"
      << "  --threads N             COPT threads per model\n"
      << "  --sp-threads N          SP-integer threads only (1..128; default: --threads)\n"
      << "  --sp-frequency N        SP/LP frequency (paper: 1500)\n"
      << "  --backup-frequency N    backup repricing frequency (paper: 1500)\n"
      << "  --pool-capacity N       main route-pool size (paper: 5000)\n"
      << "  --promising-arcs N      Vidal |Gamma| (default: 40)\n"
      << "  --granular-threshold N  disable granular topology filtering when customers <= N "
         "(default: 100)\n"
      << "  --search-profile P      paper, linked-winner, "
         "linked-insertion-candidates, or granular-adaptive\n"
      << "  --cache-capacity N      physical route-evaluation LRU capacity\n"
      << "  --time-limit SEC        cooperative soft deadline (stages may overrun)\n"
      << "  --rounded-target X      stop when objective rounded to 2 decimals <= X\n"
      << "  --sp-stall-limit SEC    engineering SP stall cap (default: 0/off)\n"
      << "  --sp-time-limit SEC     engineering SP wall cap (default: 0/off)\n"
      << "  --pool-identity P       ordered (default) or coverage experiment\n"
      << "  --path-catalog P        complete (default) or compressed experiment\n"
      << "  --implicit-replace-path  opt-in full-support random path proposals\n"
      << "  --replace-path-attempts N  random VND stage: 1 only (four-trial experiment withdrawn)\n"
      << "  --partial-route-subset P  final oversize route: prefix (default) or random (A)\n"
      << "  --insertion-paths P     Construction/repair: direct, completed, preserve-cs\n"
      << "  --initialization-paths P  Paper initialization: direct, completed, rebuilt or physical (A)\n"
      << "  --route-opening-paths P  Relocate: inherited (default), projected (A)\n"
      << "  --preserve-swap-stations  opt-in physical-slot customer Swap (Class A)\n"
      << "  --continuation P        record (default/TRC) or operational experiment\n"
      << "  --protect-best-partition  opt-in incumbent pool protection\n"
      << "  --keep-evicted-in-backup  opt-in retention of evicted main columns\n"
      << "  --reprice-after-lp       opt-in extra backup repricing after LP\n"
      << "  --adaptive-sp-stall N   opt-in early SP after N stalled ILS rounds (0 off)\n"
      << "  --adaptive-sp-min-interval N  minimum rounds between adaptive SP calls (default: 500)\n"
      << "  --exact-time-limit SEC  optional per-model cap (global is tighter)\n"
      << "  --frvcp-wall-budget SEC  optional solve-wide active FRVCP wall budget\n"
      << "  --copt-output         show COPT logs\n"
      << "  --progress N            checkpoint interval in completed iterations\n";
  ils_sp::detail::checked_flush(output);
  std::exit(status);
}

template <typename Integer>
Integer parse_integer(std::string_view value, std::string_view option) {
  return ils_sp::detail::parse_nonnegative_integer<Integer>(value, option);
}

double parse_double(std::string_view value, std::string_view option) {
  std::size_t consumed = 0;
  const double parsed = std::stod(std::string(value), &consumed);
  if (consumed != value.size() || !std::isfinite(parsed)) {
    throw std::invalid_argument("invalid value for " + std::string(option));
  }
  return parsed;
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments result;
  auto value_after = [&](int& index, std::string_view option) -> std::string_view {
    if (++index >= argc) {
      throw std::invalid_argument("missing value after " + std::string(option));
    }
    return argv[index];
  };
  for (int index = 1; index < argc; ++index) {
    const std::string_view option = argv[index];
    if (option == "--help" || option == "-h") usage(0);
    if (option == "--instance") {
      result.instance = value_after(index, option);
    } else if (option == "--iterations") {
      result.config.paper.max_iterations =
          parse_integer<std::size_t>(value_after(index, option), option);
    } else if (option == "--seed") {
      result.config.seed =
          parse_integer<std::uint64_t>(value_after(index, option), option);
    } else if (option == "--threads") {
      const int threads = parse_integer<int>(value_after(index, option), option);
      result.config.set_partitioning.threads = threads;
      result.config.frvcp.threads = threads;
    } else if (option == "--sp-threads") {
      const int threads = parse_integer<int>(value_after(index, option), option);
      if (threads < 1 || threads > 128) {
        throw std::invalid_argument("--sp-threads must be between 1 and 128");
      }
      result.config.set_partitioning.sp_integer_threads = static_cast<int>(threads);
    } else if (option == "--sp-frequency") {
      result.config.paper.sp_frequency =
          parse_integer<std::size_t>(value_after(index, option), option);
    } else if (option == "--backup-frequency") {
      result.config.paper.backup_frequency =
          parse_integer<std::size_t>(value_after(index, option), option);
    } else if (option == "--pool-capacity") {
      result.config.paper.pool_capacity =
          parse_integer<std::size_t>(value_after(index, option), option);
    } else if (option == "--promising-arcs") {
      result.config.promising_arc_count =
          parse_integer<std::size_t>(value_after(index, option), option);
    } else if (option == "--granular-threshold") {
      result.config.granular_threshold =
          parse_integer<std::size_t>(value_after(index, option), option);
    } else if (option == "--search-profile") {
      const std::string_view profile = value_after(index, option);
      if (profile == "paper") {
        result.config.search_profile = ils_sp::SearchProfile::Paper;
      } else if (profile == "linked-winner") {
        result.config.search_profile =
            ils_sp::SearchProfile::LinkedWinner;
      } else if (profile == "linked-insertion-candidates") {
        result.config.search_profile =
            ils_sp::SearchProfile::LinkedInsertionCandidates;
      } else if (profile == "granular-adaptive") {
        result.config.search_profile =
            ils_sp::SearchProfile::GranularAdaptive;
      } else {
        throw std::invalid_argument(
            "--search-profile must be paper, linked-winner, "
            "linked-insertion-candidates, or granular-adaptive");
      }
    } else if (option == "--cache-capacity") {
      result.config.evaluation_cache_capacity =
          parse_integer<std::size_t>(value_after(index, option), option);
    } else if (option == "--time-limit") {
      result.config.global_time_limit_s =
          parse_double(value_after(index, option), option);
    } else if (option == "--rounded-target") {
      result.config.rounded_target_objective_h =
          parse_double(value_after(index, option), option);
    } else if (option == "--insertion-paths") {
      const auto policy = value_after(index, option);
      result.config.preserve_station_visits = policy == "preserve-cs";
      if (policy == "direct" || policy == "preserve-cs")
        result.config.complete_insertion_paths = false;
      else if (policy == "completed") result.config.complete_insertion_paths = true;
      else throw std::invalid_argument("--insertion-paths must be direct, completed or preserve-cs");
    } else if (option == "--initialization-paths") {
      const auto policy = value_after(index, option);
      if (policy != "direct" && policy != "completed" && policy != "rebuilt" && policy != "physical") {
        throw std::invalid_argument("--initialization-paths must be direct, completed, rebuilt or physical");
      }
      result.config.complete_initialization_paths = policy != "direct";
      result.config.rebuild_initialization_paths = policy == "rebuilt";
      result.config.preserve_initialization_station_visits = policy == "physical";
    } else if (option == "--route-opening-paths") {
      const auto policy = value_after(index, option);
      if (policy != "inherited" && policy != "projected") {
        throw std::invalid_argument("--route-opening-paths must be inherited or projected");
      }
      result.config.project_opening_station_visits = policy == "projected";
    } else if (option == "--preserve-swap-stations") {
      result.config.preserve_swap_station_visits = true;
    } else if (option == "--continuation") {
      const auto policy = value_after(index, option);
      if (policy == "operational") {
        result.config.continuation_policy = ils_sp::ContinuationPolicy::Operational;
      } else if (policy == "record") {
        result.config.continuation_policy = ils_sp::ContinuationPolicy::Record;
      } else {
        throw std::invalid_argument("--continuation must be operational or record");
      }
    } else if (option == "--implicit-replace-path") {
      result.config.implicit_replace_path = true;
    } else if (option == "--replace-path-attempts") {
      result.config.replace_path_attempts =
          parse_integer<std::size_t>(value_after(index, option), option);
      if (result.config.replace_path_attempts != 1) {
        throw std::invalid_argument("--replace-path-attempts must be 1; the four-trial G1 experiment was rejected");
      }
    } else if (option == "--partial-route-subset") {
      const auto subset = value_after(index, option);
      if (subset == "prefix") result.config.partial_route_subset = ils_sp::PartialRouteSubset::Prefix;
      else if (subset == "random") result.config.partial_route_subset = ils_sp::PartialRouteSubset::Random;
      else throw std::invalid_argument("--partial-route-subset must be prefix or random");
    } else if (option == "--path-catalog") {
      const auto catalog = value_after(index, option);
      if (catalog == "complete") {
        result.config.complete_path_catalog = true;
      } else if (catalog == "compressed") {
        result.config.complete_path_catalog = false;
      } else {
        throw std::invalid_argument("--path-catalog must be complete or compressed");
      }
    } else if (option == "--pool-identity") {
      const auto identity = value_after(index, option);
      if (identity == "ordered") {
        result.config.pool_identity = ils_sp::PoolIdentity::OrderedSequence;
      } else if (identity == "coverage") {
        result.config.pool_identity = ils_sp::PoolIdentity::CustomerCoverage;
      } else {
        throw std::invalid_argument("--pool-identity must be ordered or coverage");
      }
    } else if (option == "--protect-best-partition") {
      result.config.protect_best_partition = true;
    } else if (option == "--keep-evicted-in-backup") {
      result.config.keep_evicted_in_backup = true;
    } else if (option == "--reprice-after-lp") {
      result.config.reprice_after_lp = true;
    } else if (option == "--sp-stall-limit") {
      const double limit = parse_double(value_after(index, option), option);
      result.config.sp_incumbent_stall_limit_s =
          limit == 0.0 ? std::nullopt : std::optional<double>(limit);
    } else if (option == "--sp-time-limit") {
      const double limit = parse_double(value_after(index, option), option);
      result.config.sp_maximum_time_s =
          limit == 0.0 ? std::nullopt : std::optional<double>(limit);
    } else if (option == "--adaptive-sp-stall") {
      const std::size_t stall =
          parse_integer<std::size_t>(value_after(index, option), option);
      result.config.adaptive_sp_stall_iterations =
          stall == 0 ? std::nullopt : std::optional<std::size_t>(stall);
    } else if (option == "--adaptive-sp-min-interval") {
      result.config.adaptive_sp_minimum_interval =
          parse_integer<std::size_t>(value_after(index, option), option);
    } else if (option == "--frvcp-wall-budget") {
      const double limit = parse_double(value_after(index, option), option);
      if (!std::isfinite(limit) || limit <= 0.0) {
        throw std::invalid_argument("--frvcp-wall-budget must be finite and positive");
      }
      result.config.frvcp_active_wall_budget_s = limit;
    } else if (option == "--exact-time-limit") {
      const double limit = parse_double(value_after(index, option), option);
      result.config.set_partitioning.time_limit_s = limit;
      result.config.frvcp.time_limit_s = limit;
    } else if (option == "--copt-output" || option == "--gurobi-output") {
      result.config.set_partitioning.output_flag = 1;
      result.config.frvcp.output_flag = 1;
    } else if (option == "--progress") {
      result.config.progress_interval =
          parse_integer<std::size_t>(value_after(index, option), option);
    } else {
      throw std::invalid_argument("unknown option " + std::string(option));
    }
  }
  if (result.instance.empty()) {
    throw std::invalid_argument("--instance is required");
  }
  if (result.config.global_time_limit_s) {
    (void)ils_sp::detail::checked_deadline_after<std::chrono::steady_clock>(
        std::chrono::steady_clock::now(), *result.config.global_time_limit_s);
  }
  if (result.config.sp_maximum_time_s) {
    (void)ils_sp::detail::checked_deadline_after<std::chrono::steady_clock>(
        std::chrono::steady_clock::now(), *result.config.sp_maximum_time_s);
  }
  if (result.config.rounded_target_objective_h) {
    ils_sp::detail::validate_rounded_target(*result.config.rounded_target_objective_h);
  }
  if (result.config.search_profile == ils_sp::SearchProfile::Paper &&
      result.config.preserve_station_visits &&
      !result.config.implicit_replace_path &&
      !result.config.project_opening_station_visits) {
    throw std::invalid_argument(
        "cooperative preserve-cs requires --route-opening-paths projected");
  }
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    const ils_sp::Instance instance =
        ils_sp::load_montoya_instance(arguments.instance);
    EnvrConfig environment_config;
    environment_config.Set("nobanner", "1");
    Envr environment(environment_config);
    std::cerr << "configuration external_solver=COPT external_solver_version="
              << COPT_VERSION_MAJOR << '.' << COPT_VERSION_MINOR << '.'
              << COPT_VERSION_TECHNICAL << " search_profile="
              << search_profile_name(arguments.config.search_profile)
              << " continuation="
              << (arguments.config.continuation_policy == ils_sp::ContinuationPolicy::Record
                      ? "record" : "operational")
              << " preserve_swap_stations=" << arguments.config.preserve_swap_station_visits
              << " path_catalog="
              << (arguments.config.complete_path_catalog ? "complete" : "compressed")
              << " implicit_replace_path=" << arguments.config.implicit_replace_path
              << " replace_path_attempts=" << arguments.config.replace_path_attempts
              << " partial_route_subset="
              << (arguments.config.partial_route_subset == ils_sp::PartialRouteSubset::Prefix ? "prefix" : "random")
              << " paper_insertion_paths="
              << (arguments.config.search_profile != ils_sp::SearchProfile::Paper
                      ? "profile" : arguments.config.preserve_station_visits
                                        ? "preserve-cs"
                                        : arguments.config.complete_insertion_paths.value_or(false)
                                        ? "completed" : "direct")
              << " paper_initialization_paths="
              << (arguments.config.search_profile != ils_sp::SearchProfile::Paper
                      ? "profile"
                      : arguments.config.rebuild_initialization_paths ? "rebuilt"
                      : arguments.config.preserve_initialization_station_visits ? "physical"
                      : arguments.config.complete_initialization_paths.value_or(
                            arguments.config.complete_insertion_paths.value_or(false))
                            ? "completed" : "direct")
              << " route_opening_paths="
              << (arguments.config.project_opening_station_visits ? "projected" : "inherited")
              << " pool_identity="
              << (arguments.config.pool_identity == ils_sp::PoolIdentity::CustomerCoverage
                      ? "coverage" : "ordered")
              << " protect_best_partition=" << arguments.config.protect_best_partition
              << " keep_evicted_in_backup=" << arguments.config.keep_evicted_in_backup
              << " reprice_after_lp=" << arguments.config.reprice_after_lp
              << " granular_threshold="
              << arguments.config.granular_threshold
              << " adaptive_sp_stall=";
    if (arguments.config.adaptive_sp_stall_iterations.has_value()) {
      std::cerr << *arguments.config.adaptive_sp_stall_iterations;
    } else {
      std::cerr << "off";
    }
    std::cerr << " adaptive_sp_min_interval="
              << arguments.config.adaptive_sp_minimum_interval;
    if (arguments.config.set_partitioning.sp_integer_threads) {
      std::cerr << " sp_threads="
                << *arguments.config.set_partitioning.sp_integer_threads
                << " lp_threads=" << arguments.config.set_partitioning.threads
                << " frvcp_threads=" << arguments.config.frvcp.threads;
    }
    std::cerr << " frvcp_wall_budget_s=";
    if (arguments.config.frvcp_active_wall_budget_s) {
      std::cerr << *arguments.config.frvcp_active_wall_budget_s;
    } else {
      std::cerr << "unlimited";
    }
    std::cerr << '\n';
    ils_sp::IlsSpSolver solver(instance, environment, arguments.config);
    std::optional<ils_sp::RecordImprovementProgress> last_search_record;
    const ils_sp::SolverResult result = solver.solve(
        [](const ils_sp::SolverProgress& progress) {
          std::cerr << std::fixed << std::setprecision(6)
                    << "iteration=" << progress.iteration
                    << " checkpoint="
                    << checkpoint_reason_name(progress.reason)
                    << " interval_first="
                    << progress.interval_first_iteration
                    << " interval_iterations=" << progress.interval_iterations
                    << " current=";
          if (progress.current_objective_h.has_value()) {
            std::cerr << *progress.current_objective_h;
          } else {
            std::cerr << "none";
          }
          std::cerr << " generalized=";
          if (progress.current_generalized_cost_h.has_value()) {
            std::cerr << *progress.current_generalized_cost_h;
          } else {
            std::cerr << "none";
          }
          std::cerr << " best=";
          if (progress.best_objective_h.has_value()) {
            std::cerr << *progress.best_objective_h;
          } else {
            std::cerr << "none";
          }
          std::cerr << " feasible=" << progress.current_feasible
                    << " search_profile="
                    << search_profile_name(progress.search_profile)
                    << " routes=" << progress.current_routes
                    << " lambda=" << progress.penalty_lambda
                    << " pool=" << progress.pool_size
                    << " pool_distinct_coverages="
                    << progress.pool_distinct_coverages
                    << " backup=" << progress.backup_pool_size
                    << " pool_duals_valid=" << progress.pool_duals_valid
                    << " pool_effective_generation="
                    << progress.pool_effective_generation
                    << " pool_demotions="
                    << progress.pool_statistics.main_demotions
                    << " pool_repricing_promotions="
                    << progress.pool_statistics.repricing_promotions
                    << " sp=" << progress.sp_calls
                    << " adaptive_sp=" << progress.adaptive_sp_calls
                    << " last_sp_iteration="
                    << progress.last_sp_iteration
                    << " last_best_improvement_iteration="
                    << progress.last_best_improvement_iteration
                    << " last_coverage_improvement_iteration="
                    << progress.last_coverage_improvement_iteration
                    << " lp=" << progress.lp_calls
                    << " frvcp=" << progress.frvcp_calls
                    << " frvcp_solves=" << progress.frvcp_solve_calls
                    << " frvcp_cache_hits=" << progress.frvcp_cache_hits
                    << " exact_plan_cache_size="
                    << progress.exact_plan_cache_size
                    << " exact_plan_cache_hits="
                    << progress.exact_plan_cache_hits
                    << " exact_plan_cache_misses="
                    << progress.exact_plan_cache_misses
                    << " exact_plan_cache_evictions="
                    << progress.exact_plan_cache_evictions
                    << " path_option_cache_size="
                    << progress.path_option_cache_size
                    << " path_option_cache_hits="
                    << progress.path_option_cache_hits
                    << " path_option_cache_misses="
                    << progress.path_option_cache_misses
                    << " path_option_cache_evictions="
                    << progress.path_option_cache_evictions
                    << " path_option_first_builds="
                    << progress.path_option_first_builds
                    << " path_option_rebuilds="
                    << progress.path_option_rebuilds
                    << " path_option_candidates_generated="
                    << progress.path_option_candidates_generated
                    << " path_option_nondominated_generated="
                    << progress.path_option_nondominated_generated
                    << " replacement_summary_hits="
                    << progress.replacement_summary_statistics.hits
                    << " replacement_summary_misses="
                    << progress.replacement_summary_statistics.misses
                    << " replacement_labels_built="
                    << progress.replacement_summary_statistics.labels_built
                    << " replacement_bound_certified="
                    << progress.replacement_summary_statistics.bound_certified
                    << " replacement_bound_checks="
                    << progress.replacement_summary_statistics.bound_checks
                    << " replacement_bound_pruned="
                    << progress.replacement_summary_statistics.bound_pruned
                    << " replacement_travel_hits="
                    << progress.replacement_summary_statistics.travel_hits
                    << " replacement_summary_fp_fallback="
                    << progress.replacement_summary_statistics.fp_fallback
                    << " replacement_tail_checks="
                    << progress.replacement_summary_statistics.tail_checks
                    << " replacement_raw_tail_descriptors_skipped="
                    << progress.replacement_summary_statistics.raw_tail_descriptors_skipped
                    << " replacement_partial_bound_hits="
                    << progress.replacement_summary_statistics.partial_bound_hits
                    << " replacement_propagation_pruned="
                    << progress.replacement_summary_statistics.propagation_pruned
                    << " virtual_completion_cache_size="
                    << progress.virtual_completion_cache_size
                    << " virtual_completion_cache_hits="
                    << progress.virtual_completion_cache_hits
                    << " virtual_completion_cache_misses="
                    << progress.virtual_completion_cache_misses
                    << " virtual_completion_cache_evictions="
                    << progress.virtual_completion_cache_evictions
                    << " vnd_linked_completion_attempts="
                    << progress.vnd_linked_completion_attempts
                    << " vnd_linked_completion_changed_attempts="
                    << progress.vnd_linked_completion_changed_attempts
                    << " vnd_adaptive_customer_advances="
                    << progress.vnd_adaptive_customer_advances
                    << " vnd_adaptive_replace_path_selections="
                    << progress.vnd_adaptive_replace_path_selections
                    << " vnd_adaptive_inserted_rp_failures="
                    << progress.vnd_adaptive_inserted_replace_path_failures
                    << " vnd_adaptive_final_rp_calls="
                    << progress.vnd_adaptive_final_replace_path_calls
                    << " vnd_adaptive_final_rp_gaps="
                    << progress.vnd_adaptive_final_replace_path_gaps_scanned
                    << " vnd_adaptive_final_rp_accepted="
                    << progress.vnd_adaptive_final_replace_path_accepted
                    << " initial_linked_candidates_evaluated="
                    << progress.initial_insertion_statistics
                           .candidates_evaluated
                    << " initial_rough_candidates_ranked="
                    << progress.initial_insertion_statistics
                           .candidates_ranked_by_rough_bound
                    << " initial_linked_customers_pruned_by_regret="
                    << progress.initial_insertion_statistics
                           .customers_pruned_by_regret
                    << " initial_linked_candidates_pruned_by_regret="
                    << progress.initial_insertion_statistics
                           .candidates_pruned_by_regret
                    << " initial_linked_completion_memo_hits="
                    << progress.initial_insertion_statistics
                           .completion_memo_hits
                    << " initial_linked_strict_winner_flips="
                    << progress.initial_insertion_statistics
                           .strict_winner_flips
                    << " initial_virtual_candidates="
                    << progress.initial_insertion_statistics
                           .virtual_candidates_evaluated
                    << " initial_virtual_options="
                    << progress.initial_insertion_statistics
                           .virtual_options_scored
                    << " initial_virtual_options_lb_pruned="
                    << progress.initial_insertion_statistics
                           .virtual_options_pruned_by_lower_bound
                    << " initial_virtual_strict_flips="
                    << progress.initial_insertion_statistics
                           .virtual_strict_winner_flips
                    << " repair_linked_candidates_evaluated="
                    << progress.repair_insertion_statistics
                           .candidates_evaluated
                    << " repair_granular_rejections="
                    << progress.repair_insertion_statistics
                           .candidates_rejected_by_granular
                    << " repair_rough_candidates_ranked="
                    << progress.repair_insertion_statistics
                           .candidates_ranked_by_rough_bound
                    << " repair_linked_customers_pruned_by_regret="
                    << progress.repair_insertion_statistics
                           .customers_pruned_by_regret
                    << " repair_linked_candidates_pruned_by_regret="
                    << progress.repair_insertion_statistics
                           .candidates_pruned_by_regret
                    << " repair_linked_completion_memo_hits="
                    << progress.repair_insertion_statistics
                           .completion_memo_hits
                    << " repair_linked_strict_winner_flips="
                    << progress.repair_insertion_statistics
                           .strict_winner_flips
                    << " repair_virtual_candidates="
                    << progress.repair_insertion_statistics
                           .virtual_candidates_evaluated
                    << " repair_virtual_options="
                    << progress.repair_insertion_statistics
                           .virtual_options_scored
                    << " repair_virtual_options_lb_pruned="
                    << progress.repair_insertion_statistics
                           .virtual_options_pruned_by_lower_bound
                    << " repair_virtual_candidates_lb_pruned="
                    << progress.repair_insertion_statistics
                           .virtual_candidates_pruned_by_bound
                    << " repair_virtual_customers_regret_pruned="
                    << progress.repair_insertion_statistics
                           .virtual_customers_pruned_by_regret
                    << " repair_granular_escape_selected="
                    << progress.repair_insertion_statistics
                           .granular_escape_candidates_selected
                    << " repair_virtual_strict_flips="
                    << progress.repair_insertion_statistics
                           .virtual_strict_winner_flips
                    << " repair_post_hint_attempts="
                    << progress.repair_insertion_statistics
                           .post_repair_hint_attempts
                    << " repair_post_hint_accepted="
                    << progress.repair_insertion_statistics
                           .post_repair_hint_accepted
                    << " repair_destroy_s="
                    << progress.repair_insertion_statistics.repair_timing
                           .destroy_s
                    << " repair_granular_direct_s="
                    << progress.repair_insertion_statistics.repair_timing
                           .granular_direct_scoring_s
                    << " repair_frontier_selection_s="
                    << progress.repair_insertion_statistics.repair_timing
                           .frontier_selection_s
                    << " repair_virtual_completion_s="
                    << progress.repair_insertion_statistics.repair_timing
                           .virtual_completion_s
                    << " repair_path_option_miss_build_s="
                    << progress.repair_insertion_statistics.repair_timing
                           .path_option_cache_miss_build_s
                    << " repair_nonlinear_option_scan_s="
                    << progress.repair_insertion_statistics.repair_timing
                           .nonlinear_option_scanning_s
                    << " repair_winner_materialization_s="
                    << progress.repair_insertion_statistics.repair_timing
                           .winner_materialization_s
                    << " repair_post_hint_s="
                    << progress.repair_insertion_statistics.repair_timing
                           .post_hint_s
                    << " repair_affected_completion_calls="
                    << progress.repair_insertion_statistics
                           .affected_path_completion.calls
                    << " repair_affected_completion_total_s="
                    << progress.repair_insertion_statistics
                           .affected_path_completion.total_s
                    << " repair_affected_completion_catalog_s="
                    << progress.repair_insertion_statistics
                           .affected_path_completion.option_catalog_s
                    << " repair_affected_completion_fast_s="
                    << progress.repair_insertion_statistics
                           .affected_path_completion.fast_scoring_s
                    << " repair_affected_completion_materialization_s="
                    << progress.repair_insertion_statistics
                           .affected_path_completion.full_materialization_s
                    << " interval_s=" << progress.interval_elapsed_s
                    << " initialization_s=" << progress.initialization_s
                    << " perturbation_s=" << progress.perturbation_s
                    << " vnd_s=" << progress.vnd_s
                    << " assembler_s=" << progress.assembler_s
                    << " sp_model_s=" << progress.sp_model_s
                    << " lp_model_s=" << progress.lp_model_s
                    << " frvcp_model_s=" << progress.frvcp_model_s
                    << " eval_cache_size="
                    << progress.evaluation_cache_size
                    << " eval_hits=" << progress.evaluation_cache_hits
                    << " eval_misses=" << progress.evaluation_cache_misses
                    << " eval_evictions="
                    << progress.evaluation_cache_evictions
                    << " interval_eval_hits="
                    << progress.interval_evaluation_cache_hits
                    << " interval_eval_misses="
                    << progress.interval_evaluation_cache_misses
                    << " interval_eval_evictions="
                    << progress.interval_evaluation_cache_evictions
                    << " elapsed_s=" << progress.elapsed_s << '\n';
        },
        [](const ils_sp::ExactStageProgress& progress) {
          log_exact_stage(progress);
        },
        [&last_search_record](const ils_sp::RecordImprovementProgress& record) {
          last_search_record = record;
          // Round-trip precision retains the actual double objective, not the
          // two-decimal paper display. Do not change subsequent log formatting.
          const auto flags = std::cerr.flags();
          const auto precision = std::cerr.precision();
          std::cerr << std::defaultfloat
                    << std::setprecision(std::numeric_limits<double>::max_digits10)
                    << "record_improvement iteration=" << record.iteration
                    << " objective_h=" << record.objective_h
                    << " elapsed_s=" << record.elapsed_s << '\n';
          std::cerr.flags(flags);
          std::cerr.precision(precision);
        });

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "instance=" << instance.name() << '\n'
              << "search_profile="
              << search_profile_name(arguments.config.search_profile)
              << '\n'
              << "iterations=" << result.iterations << '\n'
              << "stop_reason=" << result.stop_reason << '\n'
              << "runtime_s=" << result.runtime_s << '\n'
              << "final_lambda=" << result.final_penalty_lambda << '\n'
              << "pool_size=" << result.pool_size << '\n'
              << "pool_distinct_coverages="
              << result.pool_distinct_coverages << '\n'
              << "backup_pool_size=" << result.backup_pool_size << '\n'
              << "pool_duals_valid=" << result.pool_duals_valid << '\n'
              << "pool_effective_generation="
              << result.pool_effective_generation << '\n'
              << "pool_admission_attempts="
              << result.pool_statistics.admission_attempts << '\n'
              << "pool_dominated_rejections="
              << result.pool_statistics.dominated_rejections << '\n'
              << "pool_same_sequence_replacements="
              << result.pool_statistics.same_sequence_replacements << '\n'
              << "pool_same_coverage_replacements="
              << result.pool_statistics.same_coverage_replacements << '\n'
              << "pool_negative_rc_admissions="
              << result.pool_statistics.negative_reduced_cost_admissions << '\n'
              << "pool_main_insertions="
              << result.pool_statistics.main_insertions << '\n'
              << "pool_backup_insertions="
              << result.pool_statistics.backup_insertions << '\n'
              << "pool_main_demotions="
              << result.pool_statistics.main_demotions << '\n'
              << "pool_repricing_candidates="
              << result.pool_statistics.repricing_candidates << '\n'
              << "pool_repricing_promotions="
              << result.pool_statistics.repricing_promotions << '\n'
              << "pool_effective_coverage_improvements="
              << result.pool_statistics.effective_coverage_improvements << '\n'
              << "sp_calls=" << result.sp_calls << '\n'
              << "adaptive_sp_calls=" << result.adaptive_sp_calls << '\n'
              << "lp_calls=" << result.lp_calls << '\n'
              << "sp_lp_variant=n40_dual_start_and_rc_ub0_v2" << '\n'
              << "sp_lp_prepass_calls=" << result.sp_lp_prepass_calls << '\n'
              << "sp_lp_certified_calls=" << result.sp_lp_certified_calls << '\n'
              << "sp_lp_fallback_calls=" << result.sp_lp_fallback_calls << '\n'
              << "sp_lp_exception_calls=" << result.sp_lp_exception_calls << '\n'
              << "sp_rc_fixing_calls=" << result.sp_rc_fixing_calls << '\n'
              << "sp_rc_fixed_columns=" << result.sp_rc_fixed_columns << '\n'
              << "sp_rc_considered_columns=" << result.sp_rc_considered_columns << '\n'
              << "sp_lp_prepass_wall_s=" << result.sp_lp_prepass_wall_s << '\n'  
              << "frvcp_calls=" << result.frvcp_calls << '\n'
              << "frvcp_solve_calls=" << result.frvcp_solve_calls << '\n'
              << "frvcp_cache_hits=" << result.frvcp_cache_hits << '\n'
              << "exact_plan_cache_size=" << result.exact_plan_cache_size
              << '\n'
              << "exact_plan_cache_hits=" << result.exact_plan_cache_hits
              << '\n'
              << "exact_plan_cache_misses="
              << result.exact_plan_cache_misses << '\n'
              << "exact_plan_cache_evictions="
              << result.exact_plan_cache_evictions << '\n'
              << "path_option_cache_size=" << result.path_option_cache_size
              << '\n'
              << "path_option_cache_hits=" << result.path_option_cache_hits
              << '\n'
              << "path_option_cache_misses=" << result.path_option_cache_misses
              << '\n'
              << "path_option_cache_evictions="
              << result.path_option_cache_evictions << '\n'
              << "path_option_first_builds="
              << result.path_option_first_builds << '\n'
              << "path_option_rebuilds=" << result.path_option_rebuilds << '\n'
              << "path_option_candidates_generated="
              << result.path_option_candidates_generated << '\n'
              << "path_option_nondominated_generated="
              << result.path_option_nondominated_generated << '\n'
              << "virtual_completion_cache_size="
              << result.virtual_completion_cache_size << '\n'
              << "virtual_completion_cache_hits="
              << result.virtual_completion_cache_hits << '\n'
              << "virtual_completion_cache_misses="
              << result.virtual_completion_cache_misses << '\n'
              << "virtual_completion_cache_evictions="
              << result.virtual_completion_cache_evictions << '\n'
              << "replacement_summary_hits="
              << result.replacement_summary_statistics.hits << '\n'
              << "replacement_summary_misses="
              << result.replacement_summary_statistics.misses << '\n'
              << "replacement_labels_built="
              << result.replacement_summary_statistics.labels_built << '\n'
              << "replacement_bound_certified="
              << result.replacement_summary_statistics.bound_certified << '\n'
              << "replacement_bound_checks="
              << result.replacement_summary_statistics.bound_checks << '\n'
              << "replacement_bound_pruned="
              << result.replacement_summary_statistics.bound_pruned << '\n'
              << "replacement_travel_hits="
              << result.replacement_summary_statistics.travel_hits << '\n'
              << "replacement_summary_fp_fallback="
              << result.replacement_summary_statistics.fp_fallback << '\n'
              << "replacement_tail_checks="
              << result.replacement_summary_statistics.tail_checks << '\n'
              << "replacement_raw_tail_descriptors_skipped="
              << result.replacement_summary_statistics.raw_tail_descriptors_skipped << '\n'
              << "replacement_partial_bound_hits="
              << result.replacement_summary_statistics.partial_bound_hits << '\n'
              << "replacement_propagation_pruned="
              << result.replacement_summary_statistics.propagation_pruned << '\n'
              << "replacement_summary_groups="
              << result.replacement_summary_statistics.groups << '\n'
              << "replacement_summary_bytes="
              << result.replacement_summary_statistics.bytes << '\n'
              << "replacement_summary_peak_bytes="
              << result.replacement_summary_statistics.peak_bytes << '\n'
              << "replacement_summary_evictions="
              << result.replacement_summary_statistics.evictions << '\n'
              << "evaluation_cache_hits=" << result.evaluation_cache_hits
              << '\n'
              << "evaluation_cache_misses=" << result.evaluation_cache_misses
              << '\n'
              << "evaluation_cache_evictions="
              << result.evaluation_cache_evictions << '\n'
              << "initialization_runtime_s="
              << result.initialization_runtime_s << '\n'
              << "perturbation_runtime_s="
              << result.perturbation_runtime_s << '\n'
              << "vnd_runtime_s=" << result.vnd_runtime_s << '\n'
              << "assembler_runtime_s=" << result.assembler_runtime_s << '\n'
              << "sp_model_runtime_s=" << result.sp_model_runtime_s << '\n'
              << "lp_model_runtime_s=" << result.lp_model_runtime_s << '\n'
              << "frvcp_model_runtime_s=" << result.frvcp_model_runtime_s
              << '\n'
              << "vnd_customer_neighborhood_calls="
              << result.vnd_statistics.customer_neighborhood_calls << '\n'
              << "vnd_descriptors_considered="
              << result.vnd_statistics.descriptors_considered << '\n'
              << "vnd_route_support_groups="
              << result.vnd_statistics.route_support_groups << '\n'
              << "vnd_route_support_cache_hits="
              << result.vnd_statistics.route_support_cache_hits << '\n'
              << "vnd_route_support_buckets_skipped="
              << result.vnd_statistics.route_support_buckets_skipped << '\n'
              << "vnd_descriptors_skipped_by_memory="
              << result.vnd_statistics.descriptors_skipped_by_memory << '\n'
              << "vnd_granular_rejections="
              << result.vnd_statistics.granular_rejections << '\n'
              << "vnd_moves_evaluated="
              << result.vnd_statistics.moves_evaluated << '\n'
              << "vnd_slice_routes_evaluated="
              << result.vnd_statistics.slice_routes_evaluated << '\n'
              << "vnd_relocate_source_summary_hits="
              << result.vnd_statistics.relocate_source_summary_hits << '\n'
              << "vnd_relocate_source_summary_misses="
              << result.vnd_statistics.relocate_source_summary_misses << '\n'
              << "vnd_physical_summary_hits="
              << result.vnd_statistics.physical_summary_hits << '\n'
              << "vnd_physical_summary_misses="
              << result.vnd_statistics.physical_summary_misses << '\n'
              << "vnd_physical_summary_cache_groups="
              << result.vnd_statistics.physical_summary_cache_groups << '\n'
              << "vnd_physical_summary_cache_bytes="
              << result.vnd_statistics.physical_summary_cache_bytes << '\n'
              << "vnd_physical_summary_cache_peak_bytes="
              << result.vnd_statistics.physical_summary_cache_peak_bytes << '\n'
              << "vnd_physical_summary_cache_evictions="
              << result.vnd_statistics.physical_summary_cache_evictions << '\n'
              << "vnd_accepted_moves="
              << result.vnd_statistics.accepted_moves << '\n'
              << "vnd_linked_completion_attempts="
              << result.vnd_statistics.linked_completion_attempts << '\n'
              << "vnd_linked_completion_changed_attempts="
              << result.vnd_statistics.linked_completion_changed_attempts
              << '\n'
              << "vnd_adaptive_customer_advances="
              << result.vnd_statistics.adaptive_customer_advances << '\n'
              << "vnd_adaptive_replace_path_selections="
              << result.vnd_statistics.adaptive_replace_path_selections
              << '\n'
              << "vnd_adaptive_inserted_replace_path_failures="
              << result.vnd_statistics.adaptive_inserted_replace_path_failures
              << '\n'
              << "vnd_adaptive_final_replace_path_calls="
              << result.vnd_statistics.adaptive_final_replace_path_calls
              << '\n'
              << "vnd_adaptive_final_replace_path_gaps_scanned="
              << result.vnd_statistics.adaptive_final_replace_path_gaps_scanned
              << '\n'
              << "vnd_adaptive_final_replace_path_accepted="
              << result.vnd_statistics.adaptive_final_replace_path_accepted
              << '\n'
              << "vnd_replace_path_calls="
              << result.vnd_statistics.replace_path_calls << '\n'
              << "vnd_replace_path_accepted="
              << result.vnd_statistics.replace_path_accepted << '\n'
              << "vnd_replace_path_proposals="
              << result.vnd_statistics.replace_path_proposals << '\n'
              << "vnd_replace_path_null_proposals="
              << result.vnd_statistics.replace_path_null_proposals << '\n'
              << "vnd_replace_path_nonimproving_proposals="
              << result.vnd_statistics.replace_path_nonimproving_proposals << '\n'
              << "vnd_cooperative_vnd_calls="
              << result.vnd_statistics.cooperative_vnd_calls << '\n'
              << "vnd_cooperative_proposal_budget="
              << result.vnd_statistics.cooperative_proposal_budget << '\n'
              << "vnd_cooperative_entry_proposals="
              << result.vnd_statistics.cooperative_entry_proposals << '\n'
              << "vnd_rvnd_neighborhood_selections="
              << result.vnd_statistics.rvnd_neighborhood_selections << '\n'
              << "vnd_rvnd_improvement_resets="
              << result.vnd_statistics.rvnd_improvement_resets << '\n' 
              << "initial_linked_insertion_decisions="
              << result.initial_insertion_statistics.decisions << '\n'
              << "initial_linked_insertion_greedy_decisions="
              << result.initial_insertion_statistics.greedy_decisions << '\n'
              << "initial_linked_insertion_candidates_generated="
              << result.initial_insertion_statistics.candidates_generated
              << '\n'
              << "initial_linked_insertion_candidates_evaluated="
              << result.initial_insertion_statistics.candidates_evaluated
              << '\n'
              << "initial_insertion_candidates_ranked_by_rough_bound="
              << result.initial_insertion_statistics
                     .candidates_ranked_by_rough_bound
              << '\n'
              << "initial_linked_insertion_completion_memo_hits="
              << result.initial_insertion_statistics.completion_memo_hits
              << '\n'
              << "initial_linked_insertion_completion_memo_misses="
              << result.initial_insertion_statistics.completion_memo_misses
              << '\n'
              << "initial_linked_insertion_candidates_pruned="
              << result.initial_insertion_statistics.candidates_pruned_by_bound
              << '\n'
              << "initial_linked_insertion_candidates_pruned_by_budget="
              << result.initial_insertion_statistics
                     .candidates_pruned_by_budget
              << '\n'
              << "initial_linked_insertion_customers_pruned_by_regret="
              << result.initial_insertion_statistics.customers_pruned_by_regret
              << '\n'
              << "initial_linked_insertion_candidates_pruned_by_regret="
              << result.initial_insertion_statistics.candidates_pruned_by_regret
              << '\n'
              << "initial_linked_insertion_candidates_changed="
              << result.initial_insertion_statistics.candidates_changed
              << '\n'
              << "initial_linked_insertion_infeasible_to_feasible="
              << result.initial_insertion_statistics.infeasible_to_feasible
              << '\n'
              << "initial_linked_insertion_winner_flips="
              << result.initial_insertion_statistics.winner_flips << '\n'
              << "initial_linked_insertion_strict_winner_flips="
              << result.initial_insertion_statistics.strict_winner_flips
              << '\n'
              << "initial_linked_insertion_selected_infeasible_to_feasible="
              << result.initial_insertion_statistics
                     .selected_infeasible_to_feasible
              << '\n'
              << "initial_virtual_insertion_candidates_evaluated="
              << result.initial_insertion_statistics
                     .virtual_candidates_evaluated
              << '\n'
              << "initial_virtual_insertion_options_scored="
              << result.initial_insertion_statistics.virtual_options_scored
              << '\n'
              << "initial_virtual_insertion_options_pruned_by_lower_bound="
              << result.initial_insertion_statistics
                     .virtual_options_pruned_by_lower_bound
              << '\n'
              << "initial_virtual_insertion_memo_hits="
              << result.initial_insertion_statistics
                     .virtual_completion_memo_hits
              << '\n'
              << "initial_virtual_insertion_memo_misses="
              << result.initial_insertion_statistics
                     .virtual_completion_memo_misses
              << '\n'
              << "initial_virtual_insertion_candidates_changed="
              << result.initial_insertion_statistics
                     .virtual_candidates_changed
              << '\n'
              << "initial_virtual_insertion_infeasible_to_feasible="
              << result.initial_insertion_statistics
                     .virtual_infeasible_to_feasible
              << '\n'
              << "initial_virtual_insertion_winner_flips="
              << result.initial_insertion_statistics.virtual_winner_flips
              << '\n'
              << "initial_virtual_insertion_strict_winner_flips="
              << result.initial_insertion_statistics
                     .virtual_strict_winner_flips
              << '\n'
              << "initial_virtual_insertion_candidates_pruned_by_bound="
              << result.initial_insertion_statistics
                     .virtual_candidates_pruned_by_bound
              << '\n'
              << "initial_virtual_insertion_customers_pruned_by_regret="
              << result.initial_insertion_statistics
                     .virtual_customers_pruned_by_regret
              << '\n'
              << "initial_virtual_insertion_candidates_pruned_by_regret="
              << result.initial_insertion_statistics
                     .virtual_candidates_pruned_by_regret
              << '\n'
              << "initial_granular_escape_candidates_seen="
              << result.initial_insertion_statistics
                     .granular_escape_candidates_seen
              << '\n'
              << "initial_granular_escape_candidates_selected="
              << result.initial_insertion_statistics
                     .granular_escape_candidates_selected
              << '\n'
              << "destroy_partial_route_events="
              << result.repair_insertion_statistics.partial_route_events << '\n'
              << "destroy_randomized_partial_routes="
              << result.repair_insertion_statistics.randomized_partial_routes << '\n'
              << "repair_station_shells_retained=" << result.repair_insertion_statistics.station_shells_retained << '\n'
              << "repair_station_shells_activated=" << result.repair_insertion_statistics.station_shells_activated << '\n'
              << "repair_station_shells_discarded=" << result.repair_insertion_statistics.station_shells_discarded << '\n'
              << "repair_station_conservation_checks=" << result.repair_insertion_statistics.station_conservation_checks << '\n'
              << "repair_station_conservation_rollbacks=" << result.repair_insertion_statistics.station_conservation_rollbacks << '\n'
              << "repair_station_conservation_noops=" << result.repair_insertion_statistics.station_conservation_noops << '\n'
              << "repair_station_conservation_effective=" << result.repair_insertion_statistics.station_conservation_effective << '\n'
              << "repair_variant=" << (arguments.config.preserve_station_visits
                    ? "joint_affected_station_occurrences_v1" : "legacy") << '\n'
              << "repair_decision_replay_variant=" << (arguments.config.preserve_station_visits
                    ? "disabled_joint_station_tokens" : "epoch16_weak_identity_tape_128k") << '\n'
              << "repair_joint_repair_calls=" << result.repair_insertion_statistics.joint_repair_calls << '\n'
              << "repair_joint_removed_station_tokens=" << result.repair_insertion_statistics.joint_removed_station_tokens << '\n'
              << "repair_joint_customer_insertion_attempts=" << result.repair_insertion_statistics.joint_customer_insertion_attempts << '\n'
              << "repair_joint_station_insertion_attempts=" << result.repair_insertion_statistics.joint_station_insertion_attempts << '\n'
              << "repair_joint_blocked_token_scans=" << result.repair_insertion_statistics.joint_blocked_token_scans << '\n'
              << "repair_joint_deadlock_rollbacks=" << result.repair_insertion_statistics.joint_deadlock_rollbacks << '\n'
              << "repair_joint_invariant_rollbacks=" << result.repair_insertion_statistics.joint_invariant_rollbacks << '\n'
              << "repair_joint_max_pending_tokens=" << result.repair_insertion_statistics.joint_max_pending_tokens << '\n'
              << "repair_joint_max_removed_station_tokens=" << result.repair_insertion_statistics.joint_max_removed_station_tokens << '\n'
              << "repair_memo_hits=" << result.repair_insertion_statistics.repair_memo_hits << '\n'
              << "repair_memo_misses=" << result.repair_insertion_statistics.repair_memo_misses << '\n'
              << "repair_memo_epoch_changes=" << result.repair_insertion_statistics.repair_memo_epoch_changes << '\n'
              << "repair_memo_invalid_entries=" << result.repair_insertion_statistics.repair_memo_invalid_entries << '\n'
              << "repair_memo_evictions=" << result.repair_insertion_statistics.repair_memo_evictions << '\n'
              << "repair_memo_byte_bypasses=" << result.repair_insertion_statistics.repair_memo_byte_bypasses << '\n'
              << "repair_memo_candidates_saved=" << result.repair_insertion_statistics.repair_memo_candidates_saved << '\n'
              << "repair_memo_profiles_saved=" << result.repair_insertion_statistics.repair_memo_profiles_saved << '\n'
              << "repair_physical_insertion_decisions="
              << result.repair_insertion_statistics.physical_insertion_decisions << '\n'
              << "repair_physical_insertion_candidates="
              << result.repair_insertion_statistics.physical_insertion_candidates << '\n'
              << "repair_physical_route_profiles_built="
              << result.repair_insertion_statistics.physical_route_profiles_built << '\n'
              << "repair_linked_insertion_decisions="
              << result.repair_insertion_statistics.decisions << '\n'
              << "repair_linked_insertion_greedy_decisions="
              << result.repair_insertion_statistics.greedy_decisions << '\n'
              << "repair_linked_insertion_two_regret_decisions="
              << result.repair_insertion_statistics.two_regret_decisions
              << '\n'
              << "repair_linked_insertion_candidates_generated="
              << result.repair_insertion_statistics.candidates_generated
              << '\n'
              << "repair_linked_insertion_candidates_evaluated="
              << result.repair_insertion_statistics.candidates_evaluated
              << '\n'
              << "repair_linked_insertion_completion_memo_hits="
              << result.repair_insertion_statistics.completion_memo_hits
              << '\n'
              << "repair_linked_insertion_completion_memo_misses="
              << result.repair_insertion_statistics.completion_memo_misses
              << '\n'
              << "repair_linked_insertion_candidates_pruned="
              << result.repair_insertion_statistics.candidates_pruned_by_bound
              << '\n'
              << "repair_linked_insertion_candidates_pruned_by_budget="
              << result.repair_insertion_statistics
                     .candidates_pruned_by_budget
              << '\n'
              << "repair_linked_insertion_customers_pruned_by_regret="
              << result.repair_insertion_statistics.customers_pruned_by_regret
              << '\n'
              << "repair_linked_insertion_candidates_pruned_by_regret="
              << result.repair_insertion_statistics.candidates_pruned_by_regret
              << '\n'
              << "repair_insertion_candidates_rejected_by_granular="
              << result.repair_insertion_statistics
                     .candidates_rejected_by_granular
              << '\n'
              << "repair_insertion_candidates_ranked_by_rough_bound="
              << result.repair_insertion_statistics
                     .candidates_ranked_by_rough_bound
              << '\n'
              << "repair_linked_insertion_candidates_changed="
              << result.repair_insertion_statistics.candidates_changed
              << '\n'
              << "repair_linked_insertion_infeasible_to_feasible="
              << result.repair_insertion_statistics.infeasible_to_feasible
              << '\n'
              << "repair_linked_insertion_winner_flips="
              << result.repair_insertion_statistics.winner_flips << '\n'
              << "repair_linked_insertion_strict_winner_flips="
              << result.repair_insertion_statistics.strict_winner_flips
              << '\n'
              << "repair_linked_insertion_selected_infeasible_to_feasible="
              << result.repair_insertion_statistics
                     .selected_infeasible_to_feasible
              << '\n'
              << "repair_virtual_insertion_candidates_evaluated="
              << result.repair_insertion_statistics
                     .virtual_candidates_evaluated
              << '\n'
              << "repair_virtual_insertion_options_scored="
              << result.repair_insertion_statistics.virtual_options_scored
              << '\n'
              << "repair_virtual_insertion_options_pruned_by_lower_bound="
              << result.repair_insertion_statistics
                     .virtual_options_pruned_by_lower_bound
              << '\n'
              << "repair_virtual_insertion_memo_hits="
              << result.repair_insertion_statistics
                     .virtual_completion_memo_hits
              << '\n'
              << "repair_virtual_insertion_memo_misses="
              << result.repair_insertion_statistics
                     .virtual_completion_memo_misses
              << '\n'
              << "repair_virtual_insertion_candidates_changed="
              << result.repair_insertion_statistics
                     .virtual_candidates_changed
              << '\n'
              << "repair_virtual_insertion_infeasible_to_feasible="
              << result.repair_insertion_statistics
                     .virtual_infeasible_to_feasible
              << '\n'
              << "repair_virtual_insertion_winner_flips="
              << result.repair_insertion_statistics.virtual_winner_flips
              << '\n'
              << "repair_virtual_insertion_strict_winner_flips="
              << result.repair_insertion_statistics
                     .virtual_strict_winner_flips
              << '\n'
              << "repair_virtual_insertion_candidates_pruned_by_bound="
              << result.repair_insertion_statistics
                     .virtual_candidates_pruned_by_bound
              << '\n'
              << "repair_virtual_insertion_customers_pruned_by_regret="
              << result.repair_insertion_statistics
                     .virtual_customers_pruned_by_regret
              << '\n'
              << "repair_virtual_insertion_candidates_pruned_by_regret="
              << result.repair_insertion_statistics
                     .virtual_candidates_pruned_by_regret
              << '\n'
              << "repair_granular_escape_candidates_seen="
              << result.repair_insertion_statistics
                     .granular_escape_candidates_seen
              << '\n'
              << "repair_granular_escape_candidates_selected="
              << result.repair_insertion_statistics
                     .granular_escape_candidates_selected
              << '\n'
              << "repair_post_insertion_hint_attempts="
              << result.repair_insertion_statistics.post_repair_hint_attempts
              << '\n'
              << "repair_post_insertion_hint_accepted="
              << result.repair_insertion_statistics.post_repair_hint_accepted
              << '\n'
              << "repair_destroy_runtime_s="
              << result.repair_insertion_statistics.repair_timing.destroy_s
              << '\n'
              << "repair_granular_direct_scoring_runtime_s="
              << result.repair_insertion_statistics.repair_timing
                     .granular_direct_scoring_s
              << '\n'
              << "repair_frontier_selection_runtime_s="
              << result.repair_insertion_statistics.repair_timing
                     .frontier_selection_s
              << '\n'
              << "repair_virtual_completion_runtime_s="
              << result.repair_insertion_statistics.repair_timing
                     .virtual_completion_s
              << '\n'
              << "repair_virtual_completion_cross_cache_hits="
              << result.repair_insertion_statistics.repair_timing
                     .virtual_completion_cross_cache_hits
              << '\n'
              << "repair_virtual_completion_cross_cache_misses="
              << result.repair_insertion_statistics.repair_timing
                     .virtual_completion_cross_cache_misses
              << '\n'
              << "repair_virtual_completion_cross_cache_evictions="
              << result.repair_insertion_statistics.repair_timing
                     .virtual_completion_cross_cache_evictions
              << '\n'
              << "repair_path_option_cache_miss_build_runtime_s="
              << result.repair_insertion_statistics.repair_timing
                     .path_option_cache_miss_build_s
              << '\n'
              << "repair_path_option_cache_hits="
              << result.repair_insertion_statistics.repair_timing
                     .path_option_cache_hits
              << '\n'
              << "repair_path_option_cache_misses="
              << result.repair_insertion_statistics.repair_timing
                     .path_option_cache_misses
              << '\n'
              << "repair_path_option_cache_evictions="
              << result.repair_insertion_statistics.repair_timing
                     .path_option_cache_evictions
              << '\n'
              << "repair_path_option_candidates_generated="
              << result.repair_insertion_statistics.repair_timing
                     .path_option_candidates_generated
              << '\n'
              << "repair_path_option_nondominated="
              << result.repair_insertion_statistics.repair_timing
                     .path_option_nondominated
              << '\n'
              << "repair_nonlinear_option_scanning_runtime_s="
              << result.repair_insertion_statistics.repair_timing
                     .nonlinear_option_scanning_s
              << '\n'
              << "repair_winner_materialization_runtime_s="
              << result.repair_insertion_statistics.repair_timing
                     .winner_materialization_s
              << '\n'
              << "repair_post_hint_runtime_s="
              << result.repair_insertion_statistics.repair_timing.post_hint_s
              << '\n'
              << "repair_affected_completion_calls="
              << result.repair_insertion_statistics
                     .affected_path_completion.calls
              << '\n'
              << "repair_affected_completion_eligible_gaps="
              << result.repair_insertion_statistics
                     .affected_path_completion.eligible_gaps
              << '\n'
              << "repair_affected_completion_search_passes="
              << result.repair_insertion_statistics
                     .affected_path_completion.search_passes
              << '\n'
              << "repair_affected_completion_options_scored="
              << result.repair_insertion_statistics
                     .affected_path_completion.options_scored
              << '\n'
              << "repair_affected_completion_option_cache_hits="
              << result.repair_insertion_statistics
                     .affected_path_completion.option_cache_hits
              << '\n'
              << "repair_affected_completion_option_cache_misses="
              << result.repair_insertion_statistics
                     .affected_path_completion.option_cache_misses
              << '\n'
              << "repair_affected_completion_full_materializations="
              << result.repair_insertion_statistics
                     .affected_path_completion.full_materializations
              << '\n'
              << "repair_affected_completion_accepted_steps="
              << result.repair_insertion_statistics
                     .affected_path_completion.accepted_steps
              << '\n'
              << "repair_affected_completion_total_s="
              << result.repair_insertion_statistics
                     .affected_path_completion.total_s
              << '\n'
              << "repair_affected_completion_eligibility_s="
              << result.repair_insertion_statistics
                     .affected_path_completion.eligibility_s
              << '\n'
              << "repair_affected_completion_setup_s="
              << result.repair_insertion_statistics
                     .affected_path_completion.setup_s
              << '\n'
              << "repair_affected_completion_option_catalog_s="
              << result.repair_insertion_statistics
                     .affected_path_completion.option_catalog_s
              << '\n'
              << "repair_affected_completion_fast_scoring_s="
              << result.repair_insertion_statistics
                     .affected_path_completion.fast_scoring_s
              << '\n'
              << "repair_affected_completion_full_materialization_s="
              << result.repair_insertion_statistics
                     .affected_path_completion.full_materialization_s
              << '\n'
              << "repair_affected_completion_acceptance_s="
              << result.repair_insertion_statistics
                     .affected_path_completion.acceptance_s
              << '\n';
    std::cout << "equivalent_scoring_variant=joint_best2_forward_profile_pass_ready_v1\n"
              << "vnd_two_opt_reverse_label_slots_saved=" << result.vnd_statistics.two_opt_reverse_label_slots_saved << '\n'
              << "vnd_physical_profile_repeat_checks_saved=" << result.vnd_statistics.physical_profile_repeat_checks_saved << '\n';
    std::cout << "joint_repair_evaluation_variant=local_route_best2_v1\n"
              << "repair_joint_rank_cache_hits=" << result.repair_insertion_statistics.joint_rank_cache_hits << '\n'
              << "repair_joint_rank_cache_misses=" << result.repair_insertion_statistics.joint_rank_cache_misses << '\n'
              << "repair_joint_rank_fresh_descriptors=" << result.repair_insertion_statistics.joint_rank_fresh_descriptors << '\n'
              << "repair_joint_rank_reused_descriptors=" << result.repair_insertion_statistics.joint_rank_reused_descriptors << '\n'
              << "repair_joint_rank_score_calls=" << result.repair_insertion_statistics.joint_rank_score_calls << '\n';
    std::cout << "frvcp_active_wall_spent_s=" << result.frvcp_active_wall_spent_s << '\n'
              << "frvcp_budget_skipped_routes=" << result.frvcp_budget_skipped_routes << '\n'
              << "frvcp_budget_exhausted=" << result.frvcp_budget_exhausted << '\n';
    std::cout << "continuation_variant="
              << (arguments.config.search_profile == ils_sp::SearchProfile::Paper &&
                  arguments.config.preserve_station_visits &&
                  !arguments.config.implicit_replace_path &&
                  arguments.config.continuation_policy == ils_sp::ContinuationPolicy::Record
                      ? "trc_record_generator_fixed_vnd_anchor_no_adjacent_station_runs_no_epoch" : "legacy") << '\n'
              << "continuation_exploration_entries=" << result.continuation_exploration_entries << '\n'
              << "continuation_exploration_iterations=" << result.continuation_exploration_iterations << '\n'
              << "continuation_record_exits=" << result.continuation_record_exits << '\n'
              << "continuation_epoch_exits=" << result.continuation_epoch_exits << '\n';
    std::cout << "vnd_order_policy="
              << (arguments.config.search_profile == ils_sp::SearchProfile::Paper &&
                  arguments.config.preserve_station_visits &&
                  !arguments.config.implicit_replace_path
#if ILS_SP_TYPED_DELETE_ENTRY
                      ? "remaining_mask_rvnd_forward_only_typed_delete_entry" : "fixed_forward_only_two_opt") << '\n'
#elif ILS_SP_FORWARD_ONLY_TWO_OPT
                      ? "remaining_mask_rvnd_entry_opening_energy_proxy_uniform_ties" : "fixed_forward_only_two_opt") << '\n'
#elif ILS_SP_PHYSICAL_REVERSE_EMPTY_CHILD
                      ? "remaining_mask_rvnd_physical_reverse_empty_child" : "fixed") << '\n'
#else
                      ? "remaining_mask_rvnd" : "fixed") << '\n'
#endif
              << "search_feasible_record=" << result.best_solution.has_value() << '\n';
    // The retained feasible initialization is output-only. It cannot affect
    // search record acceptance, pool admission, warm starts or termination.
    const bool returned_initialization = result.initialization_fallback.has_value();
    const ils_sp::Solution* returned_solution = returned_initialization
        ? &*result.initialization_fallback
        : result.best_solution.has_value() ? &*result.best_solution : nullptr;
    const char* returned_origin = returned_initialization ? "initialization"
        : returned_solution ? "search_record" : "none";
    std::cout << "return_selection_variant=feasible_initialization_or_search_best_v1\n"
              << "solution_origin=" << returned_origin << '\n';
    if (returned_solution) {
      const auto flags = std::cout.flags();
      const auto precision = std::cout.precision();
      std::cout << std::defaultfloat
                << std::setprecision(std::numeric_limits<double>::max_digits10)
                << "returned_best_iteration="
                << (returned_initialization ? 0 : last_search_record->iteration) << '\n'
                << "returned_best_elapsed_s="
                << (returned_initialization ? result.initialization_runtime_s
                                            : last_search_record->elapsed_s) << '\n';
      std::cout.flags(flags);
      std::cout.precision(precision);
    } else {
      std::cout << "returned_best_iteration=none\nreturned_best_elapsed_s=none\n";
      std::cout << "objective_h=none\n"
                << "sp_skipped_incomplete_pool="
                << result.sp_skipped_incomplete_pool << '\n'
                << "pool_uncovered_customer_count="
                << result.pool_uncovered_customers.size() << '\n'
                << "pool_uncovered_customer_ids=";
      for (std::size_t index = 0;
           index < result.pool_uncovered_customers.size(); ++index) {
        if (index != 0) std::cout << ',';
        std::cout << result.pool_uncovered_customers[index];
      }
      std::cout << '\n'
                << "current_feasible=" << result.final_solution.feasible()
                << '\n'
                << "current_raw_cost_h=" << result.final_solution.raw_cost_h()
                << '\n';
      for (std::size_t index = 0;
           index < result.final_solution.plans.size(); ++index) {
        const auto& plan = result.final_solution.plans[index];
        const auto& evaluation = *plan.evaluation;
        std::cout << "current_route_" << index << '=' << instance.depot().id;
        for (const int node_id : plan.route.visits) {
          std::cout << "->" << node_id;
        }
        std::cout << "->" << instance.depot().id << '\n'
                  << "current_route_feasible_" << index << '='
                  << evaluation.feasible << '\n'
                  << "current_route_energy_shortfall_wh_" << index << '='
                  << evaluation.energy_shortfall_wh << '\n'
                  << "current_route_duration_excess_h_" << index << '='
                  << evaluation.duration_excess_h << '\n';
      }
      return 2;
    }
    std::cout << "objective_h=" << returned_solution->raw_cost_h() << '\n';
    for (std::size_t route_index = 0;
         route_index < returned_solution->plans.size(); ++route_index) {
      std::cout << "route_" << route_index << '=' << instance.depot().id;
      for (const int node_id :
           returned_solution->plans[route_index].route.visits) {
        std::cout << "->" << node_id;
      }
      std::cout << "->" << instance.depot().id << '\n';
    }
    ils_sp::detail::checked_flush(std::cout);
    return 0;
  } catch (const CoptException& error) {
    std::cerr << "COPT error " << error.GetCode() << ": "
              << error.what() << '\n';
    return 3;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
