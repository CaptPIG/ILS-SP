#pragma once

#include <array>
#include <chrono>
#include <cstdint>
#include <random>
#include <unordered_set>

#include "ils_sp/evaluator.hpp"

namespace ils_sp {

enum class DestroyOperator {
  RandomRemoval,
  RandomRoute,
  TargetRemoval,
  ClosestRemoval
};
enum class RepairOperator { GreedyInsertion, TwoRegretInsertion };
enum class Neighborhood {
  Swap,
  Relocate,
  TwoOptStar,
  ReplacePath
};
enum class SearchProfile {
  Paper,
  LinkedWinner,
  LinkedInsertionCandidates,
  GranularAdaptive
};

struct QuasiExhaustivePathResult {
  Solution solution;
  std::uint64_t gaps_scanned{};
  std::uint64_t accepted{};
};

struct LimitedRandomPathResult {
  Solution solution;
  std::uint64_t eligible_gaps{};
  std::uint64_t failure_limit{};
  std::uint64_t attempts{};
  std::uint64_t failed_attempts{};
  std::uint64_t accepted{};
};

// A lightweight Algorithm-3 endpoint used only to rank one customer
// insertion.  It records the two paths incident to the inserted customer but
// never constructs a Plan or publishes a physical-route evaluation.
struct VirtualInsertionCompletion {
  RouteCostSummary summary;
  double generalized_cost_h{};
  std::array<AnchorPair, 2> gaps;
  std::array<std::vector<int>, 2> paths;
  std::uint64_t options_scored{};
  std::uint64_t options_pruned_by_lower_bound{};
  bool changed{};
  bool infeasible_to_feasible{};
};

struct InsertionPathHint {
  AnchorPair gap;
  std::vector<int> stations;
};

struct HintedPathResult {
  Solution solution;
  std::uint64_t attempts{};
  std::uint64_t accepted{};
};

struct RepairTimingStatistics {
  std::uint64_t virtual_completion_cross_cache_hits{};
  std::uint64_t virtual_completion_cross_cache_misses{};
  std::uint64_t virtual_completion_cross_cache_evictions{};
  std::uint64_t path_option_cache_hits{};
  std::uint64_t path_option_cache_misses{};
  std::uint64_t path_option_cache_evictions{};
  std::uint64_t path_option_candidates_generated{};
  std::uint64_t path_option_nondominated{};
  double destroy_s{};
  double granular_direct_scoring_s{};
  double frontier_selection_s{};
  double virtual_completion_s{};
  double path_option_cache_miss_build_s{};
  double nonlinear_option_scanning_s{};
  double winner_materialization_s{};
  double post_hint_s{};
};

// Aggregate timings for deterministic affected-gap completion.  Timers are
// sampled once per coarse phase, never once per station-path option.
struct PathCompletionStatistics {
  std::uint64_t calls{};
  std::uint64_t eligible_gaps{};
  std::uint64_t search_passes{};
  std::uint64_t options_scored{};
  std::uint64_t option_cache_hits{};
  std::uint64_t option_cache_misses{};
  std::uint64_t full_materializations{};
  std::uint64_t accepted_steps{};
  double total_s{};
  double eligibility_s{};
  double setup_s{};
  double option_catalog_s{};
  double fast_scoring_s{};
  double full_materialization_s{};
  double acceptance_s{};
};

class GranularNeighborhood {
 public:
  using EdgeSet = std::unordered_set<std::uint64_t>;

  explicit GranularNeighborhood(const Instance& instance,
                                std::size_t maximum_neighbors = 40);

  [[nodiscard]] std::size_t gamma_size() const noexcept { return gamma_size_; }
  [[nodiscard]] const std::vector<int>& neighbors_of(int customer_id) const;
  [[nodiscard]] const std::vector<int>& promising_neighbors_of(
      int customer_id) const;
  [[nodiscard]] const std::vector<std::pair<int, int>>& promising_pairs()
      const noexcept {
    return promising_pairs_;
  }
  [[nodiscard]] std::size_t customer_rank(int customer_id) const;
  [[nodiscard]] bool promising_by_rank(std::size_t left_rank,
                                       std::size_t right_rank) const;
  [[nodiscard]] bool promising(int left_id, int right_id) const;
  // Repair insertion is granular when at least one of its newly-created
  // customer-customer arcs belongs to directed Gamma.  Route opening is
  // handled separately because it creates no customer-customer arc.
  [[nodiscard]] bool allows_insertion(std::span<const int> route,
                                      std::size_t position,
                                      int customer_id) const;
  [[nodiscard]] bool allows(
      const std::vector<std::vector<int>>& current,
      const std::vector<std::vector<int>>& candidate) const;
  [[nodiscard]] EdgeSet edge_set(
      const std::vector<std::vector<int>>& sequences) const;
  [[nodiscard]] bool allows(
      const EdgeSet& current_edges,
      const std::vector<std::vector<int>>& candidate) const;

 private:
  [[nodiscard]] static std::uint64_t edge_key(int left_id, int right_id);

  std::size_t gamma_size_{};
  std::unordered_map<int, std::size_t> customer_rank_by_id_;
  std::vector<unsigned char> promising_matrix_;
  std::unordered_map<int, std::vector<int>> directed_neighbors_;
  std::unordered_map<int, std::vector<int>> promising_neighbors_;
  std::vector<std::pair<int, int>> promising_pairs_;
  EdgeSet promising_edges_;
};

class PathSampler {
 public:
  PathSampler(const Instance& instance, PlanFactory& factory,
              std::mt19937_64& random);

  [[nodiscard]] std::optional<Solution> random_replacement(
      const Solution& solution);
  [[nodiscard]] LimitedRandomPathResult limited_random_replacement(
      Solution solution, double penalty_lambda,
      std::span<const AnchorPair> eligible_anchor_pairs,
      std::optional<std::chrono::steady_clock::time_point> deadline =
          std::nullopt);
  [[nodiscard]] QuasiExhaustivePathResult quasi_exhaustive_replacement(
      Solution solution, double penalty_lambda,
      std::span<const AnchorPair> eligible_anchor_pairs,
      std::optional<std::chrono::steady_clock::time_point> deadline =
          std::nullopt);
  [[nodiscard]] std::optional<Plan> best_replacement(
      const Plan& plan, double penalty_lambda);
  // Starting from direct paths for every new directed anchor pair, repeatedly
  // apply the best strictly improving single-gap replacement.  Paths inherited
  // from the incumbent are never eligible.
  [[nodiscard]] Plan repair_new_gaps(
      Plan plan, const PlanFactory::GapMap& inherited_gaps,
      double penalty_lambda,
      std::optional<std::chrono::steady_clock::time_point> deadline =
          std::nullopt);
  // Repair the introduced/reversed gaps and allow every inherited nonempty
  // charging path on the changed route to be shortened, removed, or replaced.
  [[nodiscard]] Plan repair_affected_paths(
      Plan plan, const PlanFactory::GapMap& inherited_gaps,
      double penalty_lambda,
      std::optional<std::chrono::steady_clock::time_point> deadline =
          std::nullopt,
      PathCompletionStatistics* statistics = nullptr);
  [[nodiscard]] std::vector<std::vector<int>> paths_between(
      int origin_id, int target_id);
  [[nodiscard]] VirtualInsertionCompletion virtual_complete_insertion(
      const Plan* source_plan, std::size_t position, int customer_id,
      double penalty_lambda,
      RepairTimingStatistics* timing = nullptr);
  [[nodiscard]] HintedPathResult apply_insertion_hints(
      Solution solution, std::span<const InsertionPathHint> hints,
      double penalty_lambda,
      std::optional<std::chrono::steady_clock::time_point> deadline =
          std::nullopt);
  [[nodiscard]] std::size_t option_cache_size() const noexcept {
    return option_cache_.size();
  }
  [[nodiscard]] std::size_t option_cache_hits() const noexcept {
    return option_cache_.hits();
  }
  [[nodiscard]] std::size_t option_cache_misses() const noexcept {
    return option_cache_.misses();
  }
  [[nodiscard]] std::size_t option_cache_evictions() const noexcept {
    return option_cache_.evictions();
  }
  [[nodiscard]] std::size_t option_first_builds() const noexcept {
    return option_first_builds_;
  }
  [[nodiscard]] std::size_t option_rebuilds() const noexcept {
    return option_rebuilds_;
  }
  [[nodiscard]] std::uint64_t option_candidates_generated() const noexcept {
    return option_candidates_generated_;
  }
  [[nodiscard]] std::uint64_t option_nondominated_generated() const noexcept {
    return option_nondominated_generated_;
  }
  [[nodiscard]] std::size_t virtual_completion_cache_size() const noexcept {
    return virtual_completion_cache_.size();
  }
  [[nodiscard]] std::size_t virtual_completion_cache_hits() const noexcept {
    return virtual_completion_cache_.hits();
  }
  [[nodiscard]] std::size_t virtual_completion_cache_misses() const noexcept {
    return virtual_completion_cache_.misses();
  }
  [[nodiscard]] std::size_t virtual_completion_cache_evictions()
      const noexcept {
    return virtual_completion_cache_.evictions();
  }

 private:
  struct PathOption {
    std::vector<int> stations;
    double total_distance_km{};
    ChargingSequenceLabel label;
  };
  struct VirtualOptionScanWorkspace {
    std::vector<double> lower_bounds;
    std::vector<double> suffix_minimum_lower_bound;
    std::vector<std::size_t> suffix_candidate_count;
  };
  struct VirtualOptionCandidate {
    double cost{};
    std::size_t side{};
    std::size_t option_index{};
  };
  struct ExcludedPathOption {
    std::size_t gap{};
    std::vector<int> stations;
  };
  using PathOptions = std::vector<PathOption>;

  [[nodiscard]] std::vector<std::vector<int>> alternatives(
      int origin_id, int target_id, const std::vector<int>& current_path);
  [[nodiscard]] std::shared_ptr<const PathOptions> path_options(
      int origin_id, int target_id,
      RepairTimingStatistics* timing = nullptr);
  [[nodiscard]] std::optional<Plan> improving_atomic_replacement(
      const Plan& plan, std::size_t gap, const std::vector<int>& stations,
      double penalty_lambda,
      std::optional<std::chrono::steady_clock::time_point> deadline);
  [[nodiscard]] std::optional<Plan> best_replacement_in_gaps(
      const Plan& plan, double penalty_lambda,
      std::span<const std::size_t> eligible_gaps,
      std::span<const ExcludedPathOption> excluded,
      std::optional<std::chrono::steady_clock::time_point> deadline,
      PathCompletionStatistics* statistics);
  [[nodiscard]] Plan improve_path_subset(
      Plan plan, double penalty_lambda,
      std::span<const std::size_t> eligible_gaps,
      std::optional<std::chrono::steady_clock::time_point> deadline,
      PathCompletionStatistics* statistics = nullptr);
  void build_curve_order();
  void build_reachable_station_index();
  void build_shortest_station_paths();

  struct StationPath {
    double distance{};
    std::vector<int> stations;
    std::size_t fastest_curve_rank{};
    std::size_t slowest_curve_rank{};
  };
  struct ReachableStation {
    std::size_t station_index{};
    double distance_km{};
  };
  struct VirtualCompletionCacheKey {
    std::vector<int> source_physical_visits;
    int customer_id{};
    std::size_t position{};
    std::uint64_t penalty_lambda_bits{};
    std::uint32_t semantics_version{1};
    bool opens_route{};
    bool operator==(const VirtualCompletionCacheKey&) const = default;
  };
  struct VirtualCompletionCacheKeyHash {
    std::size_t operator()(const VirtualCompletionCacheKey& key) const noexcept {
      std::size_t seed = IntVectorHash{}(key.source_physical_visits);
      hash_combine(seed, key.customer_id);
      hash_combine(seed, key.position);
      hash_combine(seed, key.penalty_lambda_bits);
      hash_combine(seed, key.semantics_version);
      hash_combine(seed, key.opens_route);
      return seed;
    }
  };

  const Instance& instance_;
  PlanFactory& factory_;
  std::mt19937_64& random_;
  std::unordered_map<std::string, std::size_t> curve_rank_by_type_;
  std::vector<std::size_t> station_curve_ranks_;
  std::unordered_map<int, std::vector<ReachableStation>>
      stations_reachable_from_anchor_;
  std::unordered_map<int, std::vector<ReachableStation>>
      stations_reaching_anchor_;
  std::vector<std::vector<std::vector<StationPath>>>
      shortest_station_paths_;
  std::uint64_t option_candidates_generated_{};
  std::uint64_t option_nondominated_generated_{};
  std::unordered_set<AnchorPair, AnchorPairHash> option_keys_built_;
  std::size_t option_first_builds_{};
  std::size_t option_rebuilds_{};
  LruCache<AnchorPair, std::shared_ptr<const PathOptions>, AnchorPairHash>
      option_cache_{8'192};
  LruCache<VirtualCompletionCacheKey,
           std::shared_ptr<const VirtualInsertionCompletion>,
           VirtualCompletionCacheKeyHash>
      virtual_completion_cache_{8'192};
  std::array<VirtualOptionScanWorkspace, 2> virtual_option_workspaces_;
  std::vector<VirtualOptionCandidate> virtual_tolerance_band_;
};

struct InsertionCandidateStatistics {
  std::uint64_t decisions{};
  std::uint64_t greedy_decisions{};
  std::uint64_t two_regret_decisions{};
  std::uint64_t candidates_generated{};
  std::uint64_t candidates_evaluated{};
  std::uint64_t candidates_pruned_by_bound{};
  std::uint64_t candidates_pruned_by_budget{};
  std::uint64_t customers_pruned_by_regret{};
  std::uint64_t candidates_pruned_by_regret{};
  std::uint64_t candidates_rejected_by_granular{};
  std::uint64_t candidates_ranked_by_rough_bound{};
  std::uint64_t completion_memo_hits{};
  std::uint64_t completion_memo_misses{};
  std::uint64_t candidates_changed{};
  std::uint64_t infeasible_to_feasible{};
  std::uint64_t winner_flips{};
  std::uint64_t strict_winner_flips{};
  std::uint64_t selected_infeasible_to_feasible{};
  std::uint64_t virtual_candidates_evaluated{};
  std::uint64_t virtual_options_scored{};
  std::uint64_t virtual_options_pruned_by_lower_bound{};
  std::uint64_t virtual_completion_memo_hits{};
  std::uint64_t virtual_completion_memo_misses{};
  std::uint64_t virtual_candidates_changed{};
  std::uint64_t virtual_infeasible_to_feasible{};
  std::uint64_t virtual_winner_flips{};
  std::uint64_t virtual_strict_winner_flips{};
  std::uint64_t virtual_candidates_pruned_by_bound{};
  std::uint64_t virtual_customers_pruned_by_regret{};
  std::uint64_t virtual_candidates_pruned_by_regret{};
  std::uint64_t granular_escape_candidates_seen{};
  std::uint64_t granular_escape_candidates_selected{};
  std::uint64_t post_repair_hint_attempts{};
  std::uint64_t post_repair_hint_accepted{};
  RepairTimingStatistics repair_timing;
  PathCompletionStatistics affected_path_completion;
};

class InitialSolutionBuilder {
 public:
  InitialSolutionBuilder(const Instance& instance, PlanFactory& factory,
                         PathSampler& path_sampler, std::mt19937_64& random,
                         SearchProfile profile = SearchProfile::Paper,
                         std::size_t promising_arc_count = 40,
                         std::size_t granular_threshold = 100);
  [[nodiscard]] Solution build(double penalty_lambda);
  [[nodiscard]] const InsertionCandidateStatistics& statistics() const
      noexcept {
    return statistics_;
  }

 private:
  const Instance& instance_;
  PlanFactory& factory_;
  PathSampler& path_sampler_;
  std::mt19937_64& random_;
  SearchProfile profile_;
  std::unique_ptr<GranularNeighborhood> granular_;
  double maximum_charging_rate_wh_per_h_{};
  InsertionCandidateStatistics statistics_;
};

struct DestroyedSolution {
  Solution remaining;
  std::vector<int> removed_customers;
};

class Perturbation {
 public:
  Perturbation(const Instance& instance, PlanFactory& factory,
               PathSampler& path_sampler, std::mt19937_64& random,
               SearchProfile profile = SearchProfile::Paper,
               std::size_t promising_arc_count = 40,
               std::size_t granular_threshold = 100);
  [[nodiscard]] Solution apply(const Solution& solution,
                               double penalty_lambda);
  [[nodiscard]] const InsertionCandidateStatistics& statistics() const
      noexcept {
    return statistics_;
  }

 private:
  friend struct PerturbationTestAccess;

  [[nodiscard]] DestroyedSolution destroy(
      const Solution& solution, DestroyOperator operation,
      std::size_t removal_count, double penalty_lambda);
  [[nodiscard]] Solution repair(DestroyedSolution destroyed,
                                RepairOperator operation,
                                double penalty_lambda,
                                std::vector<InsertionPathHint>* hints =
                                    nullptr);

  const Instance& instance_;
  PlanFactory& factory_;
  PathSampler& path_sampler_;
  std::mt19937_64& random_;
  SearchProfile profile_;
  std::unique_ptr<GranularNeighborhood> granular_;
  double maximum_charging_rate_wh_per_h_{};
  InsertionCandidateStatistics statistics_;
  std::unique_ptr<RouteSequenceScorer> sequence_scorer_;
};

struct VndStatistics {
  std::uint64_t customer_neighborhood_calls{};
  std::uint64_t descriptors_considered{};
  std::uint64_t route_support_groups{};
  std::uint64_t route_support_cache_hits{};
  std::uint64_t route_support_buckets_skipped{};
  std::uint64_t descriptors_skipped_by_memory{};
  std::uint64_t granular_rejections{};
  std::uint64_t moves_evaluated{};
  std::uint64_t slice_routes_evaluated{};
  std::uint64_t relocate_source_summary_hits{};
  std::uint64_t relocate_source_summary_misses{};
  std::uint64_t accepted_moves{};
  // VND winner attempts only; initialization and perturbation are not counted.
  std::uint64_t linked_completion_attempts{};
  std::uint64_t linked_completion_changed_attempts{};
  std::uint64_t replace_path_calls{};
  std::uint64_t replace_path_accepted{};
  std::uint64_t adaptive_customer_advances{};
  std::uint64_t adaptive_replace_path_selections{};
  std::uint64_t adaptive_inserted_replace_path_failures{};
  std::uint64_t adaptive_final_replace_path_calls{};
  std::uint64_t adaptive_final_replace_path_gaps_scanned{};
  std::uint64_t adaptive_final_replace_path_accepted{};
};

class VariableNeighborhoodDescent {
 public:
  VariableNeighborhoodDescent(const Instance& instance, PlanFactory& factory,
                              PathSampler& path_sampler,
                              std::mt19937_64& random,
                              std::size_t promising_arc_count = 40,
                              SearchProfile profile = SearchProfile::Paper);

  [[nodiscard]] Solution improve(
      Solution solution, double penalty_lambda,
      std::optional<std::chrono::steady_clock::time_point> deadline =
          std::nullopt);
  [[nodiscard]] const GranularNeighborhood& granular() const noexcept {
    return granular_;
  }
  [[nodiscard]] const VndStatistics& statistics() const noexcept {
    return statistics_;
  }

 private:
  static constexpr std::size_t kMaximumMoveSlices = 5;

  struct DescribedRouteChange {
    std::size_t route{};
    std::array<IncumbentRouteSlice, kMaximumMoveSlices> slices{};
    std::size_t slice_count{};
  };

  struct MoveDescriptor {
    std::array<DescribedRouteChange, 2> changes{};
    std::size_t change_count{};
  };

  struct MoveMemoryKey {
    Neighborhood neighborhood{};
    // Move descriptors store absolute route slots.  A content version can
    // reappear in another slot after inter-route moves, so both parts are
    // required to identify a support whose cached descriptor is applicable.
    std::array<std::size_t, 2> route_indices{};
    std::array<std::uint64_t, 2> route_versions{};
    std::size_t route_count{};
    bool opens_route{};
    bool operator==(const MoveMemoryKey&) const = default;
  };

  struct MoveMemoryKeyHash {
    std::size_t operator()(const MoveMemoryKey& key) const noexcept {
      std::size_t seed = static_cast<std::size_t>(key.neighborhood);
      hash_combine(seed, key.opens_route);
      hash_combine(seed, key.route_count);
      for (std::size_t index = 0; index < key.route_count; ++index) {
        hash_combine(seed, key.route_indices[index]);
        hash_combine(seed, key.route_versions[index]);
      }
      return seed;
    }
  };

  struct MoveMemoryEntry {
    // Cache the fast cost of the candidate support, not its delta.  The same
    // physical visits can later carry a different charging evaluation, so the
    // incumbent support cost is subtracted again on every lookup.
    double candidate_cost{};
    MoveDescriptor move;
  };

  enum class SupportState : unsigned char { Unknown, Open, Closed };

  struct SupportWorkingMemory {
    std::uint64_t epoch{};
    SupportState state{SupportState::Unknown};
    bool active{};
    MoveMemoryKey key{};
    MoveMemoryEntry entry{};
  };

  [[nodiscard]] std::optional<Solution> best_neighbor(
      const Solution& solution, Neighborhood neighborhood,
      double penalty_lambda,
      std::optional<std::chrono::steady_clock::time_point> deadline);
  [[nodiscard]] std::optional<Solution> best_customer_neighbor(
      const Solution& solution, Neighborhood neighborhood,
      double penalty_lambda,
      std::optional<std::chrono::steady_clock::time_point> deadline);
  [[nodiscard]] std::optional<Solution> replace_path_neighbor(
      const Solution& solution, double penalty_lambda);
  [[nodiscard]] Solution inserted_replace_path(
      Solution solution, double penalty_lambda,
      std::span<const AnchorPair> dirty_anchor_pairs,
      std::optional<std::chrono::steady_clock::time_point> deadline);
  [[nodiscard]] Solution final_replace_path(
      Solution solution, double penalty_lambda,
      std::optional<std::chrono::steady_clock::time_point> deadline);
  [[nodiscard]] bool expired(
      std::optional<std::chrono::steady_clock::time_point> deadline) const;

  const Instance& instance_;
  PlanFactory& factory_;
  GranularNeighborhood granular_;
  PathSampler& path_sampler_;
  std::mt19937_64& random_;
  SearchProfile profile_;
  LruCache<MoveMemoryKey, MoveMemoryEntry, MoveMemoryKeyHash> move_memory_{
      100'000};
  std::unordered_map<std::vector<int>, std::uint64_t, IntVectorHash>
      route_versions_;
  std::uint64_t next_route_version_{1};
  std::unique_ptr<RouteSequenceScorer> sequence_scorer_;
  std::vector<SupportWorkingMemory> support_working_memory_;
  std::uint64_t support_epoch_{};
  std::vector<std::pair<int, int>> swap_customer_pairs_;
  std::vector<std::tuple<int, std::size_t, std::size_t>> relocate_moves_;
  std::vector<
      std::tuple<std::size_t, std::size_t, std::size_t, std::size_t, bool>>
      two_opt_moves_;
  VndStatistics statistics_;
  const std::array<Neighborhood, 4> order_{Neighborhood::Swap,
                                           Neighborhood::Relocate,
                                           Neighborhood::TwoOptStar,
                                           Neighborhood::ReplacePath};
};

[[nodiscard]] Solution build_solution_from_sequences(
    const std::vector<std::vector<int>>& sequences, PlanFactory& factory,
    const Solution* reference);

}  // namespace ils_sp
