#pragma once

#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <new>
#include <type_traits>
#include <random>
#include <unordered_set>
#include <variant>

#include "ils_sp/evaluator.hpp"

namespace ils_sp {

enum class DestroyOperator {
  RandomRemoval,
  RandomRoute,
  TargetRemoval,
  ClosestRemoval
};
enum class RepairOperator { GreedyInsertion, TwoRegretInsertion };
enum class PartialRouteSubset { Prefix, Random };
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

// Optional coarse timings for the standalone catalog diagnostic. No clocks
// are sampled for these phases when the caller does not request them.
struct PathCatalogStatistics {
  bool cache_hit{};
  double setup_s{};
  double enumeration_s{};
  double candidate_views_s{};
  double candidate_sort_s{};
  double dominance_filter_s{};
  double workspace_release_s{};
  double descriptor_materialization_s{};
  double label_materialization_s{};
  double option_sort_s{};
  double cache_admission_s{};
  double cleanup_s{};
  double export_s{};
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
  bool consecutive_customer_ids_{};
  std::vector<unsigned char> promising_matrix_;
  std::unordered_map<int, std::vector<int>> directed_neighbors_;
  std::unordered_map<int, std::vector<int>> promising_neighbors_;
  std::vector<std::pair<int, int>> promising_pairs_;
  EdgeSet promising_edges_;
};

// Keep the public diagnostic schema. With cross-call summary storage removed,
// hits, groups, bytes, peak_bytes, evictions, travel_hits and partial_bound_hits
// remain zero; the remaining fields still report scoring and bound work.
struct ReplacementSummaryStatistics {
  std::uint64_t hits{};
  std::uint64_t misses{};
  std::uint64_t labels_built{};
  std::size_t groups{};
  std::size_t bytes{};
  std::size_t peak_bytes{};
  std::size_t evictions{};
  bool bound_certified{};
  std::uint64_t bound_checks{};
  std::uint64_t bound_pruned{};
  std::uint64_t travel_hits{};
  bool fp_fallback{};
  std::uint64_t tail_checks{};
  // Raw descriptors, including current/excluded paths; not eligible scores.
  std::uint64_t raw_tail_descriptors_skipped{};
  std::uint64_t partial_bound_hits{};
  std::uint64_t propagation_pruned{};
};

class PathSampler {
 public:
  // Lower bound for a new customer inserted into a source route. Unsupported
  // numerical/route domains return -infinity, disabling linked pruning.
  [[nodiscard]] double linked_insertion_lower_bound(
      const RouteCostSummary& summary, double penalty_lambda,
      std::size_t source_visits, std::size_t source_customers) const;
  PathSampler(const Instance& instance, PlanFactory& factory,
              std::mt19937_64& random, bool complete_paths = true,
              bool implicit_replacement = false);

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
  // Experimental single_step selects one existing summary-best proposal,
  // full-verifies it, and stops even if that proposal fails full verification.
  // The default retains the original fixed-point completion everywhere.
  [[nodiscard]] Plan repair_new_gaps(
      Plan plan, const PlanFactory::GapMap& inherited_gaps,
      double penalty_lambda,
      std::optional<std::chrono::steady_clock::time_point> deadline =
          std::nullopt,
      bool single_step = false);
  // Repair the introduced/reversed gaps and allow every inherited nonempty
  // charging path on the changed route to be shortened, removed, or replaced.
  [[nodiscard]] Plan repair_affected_paths(
      Plan plan, const PlanFactory::GapMap& inherited_gaps,
      double penalty_lambda,
      std::optional<std::chrono::steady_clock::time_point> deadline =
          std::nullopt,
      PathCompletionStatistics* statistics = nullptr);
  [[nodiscard]] std::vector<std::vector<int>> paths_between(
      int origin_id, int target_id,
      PathCatalogStatistics* statistics = nullptr);
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
    return catalog_cache_.size();
  }
  [[nodiscard]] std::size_t option_cache_hits() const noexcept {
    return catalog_cache_.hits();
  }
  [[nodiscard]] std::size_t option_cache_misses() const noexcept {
    return catalog_cache_.misses();
  }
  [[nodiscard]] std::size_t option_cache_evictions() const noexcept {
    return catalog_cache_.evictions();
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
  [[nodiscard]] const ReplacementSummaryStatistics&
  replacement_summary_statistics() const noexcept {
    return replacement_summary_statistics_;
  }

 private:
  friend class VariableNeighborhoodDescent;
  struct BoundedProxyChoice {
    std::size_t position{};
    std::vector<int> stations;
    double proxy_cost{};
  };
  struct BoundedProxyScan {
    std::optional<BoundedProxyChoice> choice;
    std::size_t raw_descriptors{}, work{};
    // completed=false means discard every provisional choice and fall back.
    bool completed{};
  };
  [[nodiscard]] BoundedProxyScan bounded_proxy_choice(
      const Solution& solution,
      std::span<const std::pair<std::size_t, std::size_t>> gaps,
      std::size_t first, double penalty_lambda,
      std::size_t descriptor_cap = 256, std::size_t work_cap = 4096);
  struct RandomReplacementProposal {
    std::size_t route_index{};
    Plan plan;
  };
  [[nodiscard]] std::optional<RandomReplacementProposal>
  random_replacement_proposal(const Solution& solution);

  [[nodiscard]] static double catalog_travel_lower_bound(double distance, double speed);
  struct PathDescriptor {
    // Packed catalogs store MSB-first nonzero four-bit station codes; flat
    // catalogs store the original offset in the low 32 bits and size above it.
    std::uint64_t station_storage{};
    double total_distance_km{};
  };
  static_assert(sizeof(PathDescriptor) == 16);
  struct CatalogDiskBacking;
  // Keep the frozen vector layout separately so the experimental reader does
  // not change the original logical LRU charge or catalogue identity schedule.
  struct PathCatalogPayload {
    using DecodeBuffer = std::array<int, 15>;
    std::vector<PathDescriptor> paths;
    std::vector<int> station_ids;
    std::vector<std::uint32_t> lexical_order;
    DecodeBuffer packed_station_ids{};
    std::uint8_t packed_station_count{};
    bool packed{};
  };
  struct PathCatalog : PathCatalogPayload {
    std::shared_ptr<const CatalogDiskBacking> disk_backing;
    [[nodiscard]] std::uint64_t encode_path(
        std::span<const int> stations, std::uint32_t flat_offset) const;
    [[nodiscard]] std::size_t path_size(const PathDescriptor& path) const noexcept;
    [[nodiscard]] std::span<const int> stations(
        const PathDescriptor& path, DecodeBuffer& decoded) const;
    [[nodiscard]] bool lexically_less(
        const PathDescriptor& left, const PathDescriptor& right) const;
    [[nodiscard]] std::size_t retained_bytes() const noexcept;
    [[nodiscard]] std::size_t size() const noexcept;
  };
  struct CatalogBlockKey {
    std::uint64_t record{}, offset{};
    bool operator==(const CatalogBlockKey&) const = default;
  };
  struct CatalogBlockKeyHash {
    std::size_t operator()(const CatalogBlockKey& key) const noexcept {
      std::size_t seed = std::hash<std::uint64_t>{}(key.record);
      hash_combine(seed, key.offset);
      return seed;
    }
  };
  struct CatalogPageState {
    using Page = std::vector<unsigned char>;
    LruCache<CatalogBlockKey, std::shared_ptr<const Page>, CatalogBlockKeyHash>
        pages{std::numeric_limits<std::size_t>::max()};
    std::size_t budget{64 * 1024 * 1024}, bytes{}, peak_bytes{};
    std::size_t active_readers{}, peak_readers{};
    std::uint64_t hits{}, misses{}, bytes_read{}, failures{}, audits{};
    double read_s{};
  };
  struct CatalogDiskBacking {
    std::array<std::uint64_t, 32> header{};
    std::string file;
    using BlockChecksum = std::array<std::uint64_t, 2>;
    std::vector<BlockChecksum> block_hashes;
    std::shared_ptr<CatalogPageState> state;
    std::uint64_t record_id{};
  };
  // Short-lived, noncopyable reader. Scalar values are copied out of verified
  // pages, so callers cannot retain a pointer into an evictable cache block.
  struct CatalogReader {
    explicit CatalogReader(const PathCatalog& catalog);
    ~CatalogReader();
    CatalogReader(const CatalogReader&) = delete;
    CatalogReader& operator=(const CatalogReader&) = delete;
    [[nodiscard]] PathDescriptor descriptor(std::size_t index) const;
    [[nodiscard]] std::uint32_t lexical_index(std::size_t ordinal) const;
    [[nodiscard]] std::pair<std::size_t, std::size_t> equal_range(
        const std::vector<int>& path) const;
    void verify_all() const;
    const PathCatalog& catalog;
   private:
    [[noreturn]] void fail(const char* reason) const;
    void copy_value(void* destination, std::size_t width, std::uint64_t section,
                    std::uint64_t section_bytes, std::uint64_t relative,
                    std::size_t hash_base) const;
    std::shared_ptr<const CatalogDiskBacking> backing;
    int fd{-1};
  };
  // Optional, run-local backing store. A fresh explicit workspace directory
  // owns these files; ordinary runs retain the original RAM-only behavior.
  struct CatalogDiskCache {
    using Header = std::array<std::uint64_t, 32>;
    struct Record {
      Header header{};
      bool complete{};
      std::shared_ptr<const CatalogDiskBacking> backing;
    };
    std::filesystem::path root;
    std::unordered_map<AnchorPair, Record, AnchorPairHash> records;
    std::uint64_t byte_limit{8ULL * 1024 * 1024 * 1024};
    std::uint64_t metadata_limit{128ULL * 1024 * 1024}, metadata_bytes{};
    std::shared_ptr<CatalogPageState> page_state;
    std::uint64_t charged_bytes{}, bytes_read{}, bytes_written{};
    std::uint64_t hits{}, misses{}, stores{}, failures{}, quota_skips{}, fp_skips{};
    double read_s{}, write_s{};
    bool requested{}, enabled{};
    ~CatalogDiskCache();
  };
  void initialize_catalog_disk_cache();
  [[nodiscard]] std::shared_ptr<const PathCatalog> disk_cache_load(AnchorPair key);
  void disk_cache_store(AnchorPair key, const PathCatalog& catalog);
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

  [[nodiscard]] std::optional<std::vector<int>> sample_explicit_path(
      int origin_id, int target_id, const std::vector<int>& current_path);
  [[nodiscard]] std::shared_ptr<const PathCatalog> path_catalog(
      int origin_id, int target_id,
      RepairTimingStatistics* timing = nullptr,
      PathCatalogStatistics* statistics = nullptr);
  [[nodiscard]] std::shared_ptr<const PathOptions> path_options(
      int origin_id, int target_id,
      RepairTimingStatistics* timing = nullptr,
      PathCatalogStatistics* statistics = nullptr);
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
  bool complete_paths_;
  bool implicit_replacement_;
  // Representation oracle controlled only by PathSamplerTestAccess.
  bool packed_catalogs_enabled_{true};
  // Construction representation oracle, independent of final catalog packing.
  bool packed_construction_enabled_{true};
  // Direct-contract oracle; no public search option or alternate heuristic.
  bool catalog_first_leg_index_enabled_{true};
  bool zero_shift_energy_bound_enabled_{true};
  [[nodiscard]] static std::optional<double> zero_shift_energy_cost_bound(
      double travel, double service, const Vehicle& vehicle, double lambda);
  const double maximum_charging_rate_wh_per_h_;
  std::optional<double> linked_triangle_slack_h_;
  friend struct PathSamplerTestAccess;
  [[nodiscard]] bool path_dominated(int origin, int target,
                                     const std::vector<int>& path) const;
  [[nodiscard]] std::optional<std::vector<int>> sample_implicit_path(
      int origin, int target, const std::vector<int>& current);
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
  // Immutable physical catalogs are shared by export, uniform sampling and
  // scoring. Large charging labels have a separate 128 MiB payload cache;
  // ordinary random ReplacePath never constructs those labels.
  LruCache<AnchorPair, std::shared_ptr<const PathCatalog>, AnchorPairHash>
      catalog_cache_{std::numeric_limits<std::size_t>::max()};
  std::size_t catalog_cache_bytes_{};
  std::size_t catalog_cache_budget_{1024 * 1024 * 1024};
  CatalogDiskCache catalog_disk_cache_;
  LruCache<AnchorPair, std::shared_ptr<const PathOptions>, AnchorPairHash>
      option_cache_{std::numeric_limits<std::size_t>::max()};
  std::size_t option_cache_bytes_{};
  // Retain the historical numerical-bound gate: zero disables replacement
  // bounds. No summary storage is allocated from this value.
  std::size_t replacement_summary_budget_{32 * 1024 * 1024};
  ReplacementSummaryStatistics replacement_summary_statistics_;
  bool bulk_replacement_bound_enabled_{true};
  bool shortfall_replacement_bound_enabled_{true};
  LruCache<VirtualCompletionCacheKey,
           std::shared_ptr<const VirtualInsertionCompletion>,
           VirtualCompletionCacheKeyHash>
      virtual_completion_cache_{8'192};
  std::array<VirtualOptionScanWorkspace, 2> virtual_option_workspaces_;
  std::vector<VirtualOptionCandidate> virtual_tolerance_band_;
};

struct InsertionCandidateStatistics {
  // Joint-branch attempted work, including work discarded by whole-S fallback.
  std::uint64_t joint_repair_calls{};
  // Local rank reuse: physical_insertion_candidates remains logical positions.
  std::uint64_t joint_rank_cache_hits{};
  std::uint64_t joint_rank_cache_misses{};
  std::uint64_t joint_rank_fresh_descriptors{};
  std::uint64_t joint_rank_reused_descriptors{};
  std::uint64_t joint_rank_score_calls{};
  std::uint64_t joint_removed_station_tokens{};
  std::uint64_t joint_customer_insertion_attempts{};
  std::uint64_t joint_station_insertion_attempts{};
  std::uint64_t joint_blocked_token_scans{};
  std::uint64_t joint_deadlock_rollbacks{};
  std::uint64_t joint_invariant_rollbacks{};
  std::uint64_t joint_max_pending_tokens{};
  std::uint64_t joint_max_removed_station_tokens{};

  std::uint64_t station_shells_retained{};
  std::uint64_t station_shells_activated{};
  std::uint64_t station_shells_discarded{};
  std::uint64_t station_conservation_checks{};
  std::uint64_t station_conservation_rollbacks{};
  std::uint64_t station_conservation_noops{};
  std::uint64_t station_conservation_effective{};
  std::uint64_t physical_insertion_decisions{};
  std::uint64_t physical_insertion_candidates{};
  std::uint64_t physical_route_profiles_built{};
  std::uint64_t repair_memo_hits{};
  std::uint64_t repair_memo_misses{};
  std::uint64_t repair_memo_epoch_changes{};
  std::uint64_t repair_memo_invalid_entries{};
  std::uint64_t repair_memo_evictions{};
  std::uint64_t repair_memo_byte_bypasses{};
  std::uint64_t repair_memo_candidates_saved{};
  std::uint64_t repair_memo_profiles_saved{};

  std::uint64_t partial_route_events{};
  std::uint64_t randomized_partial_routes{};
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
                         std::size_t granular_threshold = 100,
                         bool rebuild_station_paths = false,
                         bool preserve_station_visits = false);
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
  bool rebuild_station_paths_;
  bool preserve_station_visits_;
  std::unique_ptr<GranularNeighborhood> granular_;
  InsertionCandidateStatistics statistics_;
};

// Each pending station is one original visit, not a unique physical station ID.
struct StationOccurrenceToken {
  int station_id{};
  std::size_t source_route{};
  std::size_t source_visit{};
};

struct DestroyedSolution {
  Solution remaining;
  std::vector<int> removed_customers;
  // Dormant physical CS chains from customer-empty routes. Not fleet members.
  std::vector<std::vector<int>> station_shells;
  std::vector<StationOccurrenceToken> removed_stations;
};

class Perturbation {
 public:
  Perturbation(const Instance& instance, PlanFactory& factory,
               PathSampler& path_sampler, std::mt19937_64& random,
               SearchProfile profile = SearchProfile::Paper,
               std::size_t promising_arc_count = 40,
               std::size_t granular_threshold = 100,
               PartialRouteSubset partial_route_subset = PartialRouteSubset::Prefix,
               bool preserve_station_visits = false);
  ~Perturbation();
  [[nodiscard]] Solution apply(const Solution& solution,
                               double penalty_lambda);
  [[nodiscard]] std::size_t last_removal_count() const noexcept {
    return last_removal_count_;
  }
  [[nodiscard]] const InsertionCandidateStatistics& statistics() const
      noexcept {
    return statistics_;
  }

 private:
  friend struct PerturbationTestAccess;
  struct RepairDecisionMemo;

  [[nodiscard]] DestroyedSolution destroy(
      const Solution& solution, DestroyOperator operation,
      std::size_t removal_count, double penalty_lambda);
  [[nodiscard]] Solution repair(DestroyedSolution destroyed,
                                RepairOperator operation,
                                double penalty_lambda,
                                std::vector<InsertionPathHint>* hints =
                                    nullptr,
                                RepairDecisionMemo* memo = nullptr);
  [[nodiscard]] DestroyedSolution destroy_joint_occurrences(
      const Solution& solution, std::vector<int> removed_customers);
  [[nodiscard]] Solution repair_joint_occurrences(
      const Solution& original, DestroyedSolution destroyed,
      RepairOperator operation, double penalty_lambda);
  std::string last_joint_failure_;
  std::uint64_t joint_station_insertions_{};
  std::uint64_t joint_customer_insertions_{};
  std::uint64_t joint_blocked_token_scans_{};

  [[nodiscard]] Solution enforce_station_conservation(
      const Solution& original, Solution repaired,
      const std::vector<std::vector<int>>& expected_station_chains);
  [[nodiscard]] Solution repair_preserving_stations(
      DestroyedSolution destroyed, RepairOperator operation,
      double penalty_lambda, RepairDecisionMemo* memo = nullptr);

  const Instance& instance_;
  PlanFactory& factory_;
  PathSampler& path_sampler_;
  std::mt19937_64& random_;
  SearchProfile profile_;
  PartialRouteSubset partial_route_subset_;
  bool preserve_station_visits_{};
  std::size_t last_removal_count_{};
  // Isolated full-gate candidate: enable repair-only one-step completion.
  bool single_step_repair_completion_{true};
  std::unique_ptr<GranularNeighborhood> granular_;
  InsertionCandidateStatistics statistics_;
  std::unique_ptr<RouteSequenceScorer> sequence_scorer_;
  std::unique_ptr<RepairDecisionMemo> repair_decision_memo_;
};

struct VndStatistics {
  std::uint64_t two_opt_reverse_label_slots_saved{};
  std::uint64_t physical_profile_repeat_checks_saved{};
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
  std::uint64_t physical_summary_hits{};
  std::uint64_t physical_summary_misses{};
  std::uint64_t physical_summary_cache_groups{};
  std::uint64_t physical_summary_cache_bytes{};
  std::uint64_t physical_summary_cache_peak_bytes{};
  std::uint64_t physical_summary_cache_evictions{};
  std::uint64_t accepted_moves{};
  // VND winner attempts only; initialization and perturbation are not counted.
  std::uint64_t linked_completion_attempts{};
  std::uint64_t linked_completion_changed_attempts{};
  std::uint64_t replace_path_calls{};
  std::uint64_t replace_path_accepted{};
  // Non-adaptive VND proposals, including null implicit proposals.
  std::uint64_t replace_path_proposals{};
  std::uint64_t replace_path_null_proposals{};
  std::uint64_t replace_path_nonimproving_proposals{};
  // Experimental fixed-budget schedule; zero for legacy calls without kappa.
  std::uint64_t cooperative_vnd_calls{};
  std::uint64_t cooperative_proposal_budget{};
  std::uint64_t cooperative_entry_proposals{};
  std::uint64_t rvnd_neighborhood_selections{};
  std::uint64_t rvnd_improvement_resets{};
  std::uint64_t adaptive_customer_advances{};
  std::uint64_t adaptive_replace_path_selections{};
  std::uint64_t adaptive_inserted_replace_path_failures{};
  std::uint64_t adaptive_final_replace_path_calls{};
  std::uint64_t adaptive_final_replace_path_gaps_scanned{};
  std::uint64_t adaptive_final_replace_path_accepted{};
};

class VariableNeighborhoodDescent {
  friend struct VndTestAccess;

 public:
  VariableNeighborhoodDescent(const Instance& instance, PlanFactory& factory,
                              PathSampler& path_sampler,
                              std::mt19937_64& random,
                              std::size_t promising_arc_count = 40,
                              SearchProfile profile = SearchProfile::Paper,
                              std::size_t replace_path_attempts = 1,
                              bool project_opening_station_visits = false,
                              bool preserve_swap_station_visits = false);

  // Supplying the preceding perturbation's actual removal count selects the
  // experimental cooperative RP schedule (Paper + explicit sampler only).
  [[nodiscard]] Solution improve(
      Solution solution, double penalty_lambda,
      std::optional<std::chrono::steady_clock::time_point> deadline =
          std::nullopt,
      std::optional<std::size_t> cooperative_removal_count = std::nullopt);
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
    bool physical_reverse_two_opt{};
    bool physical_forward_two_opt{};
    std::array<std::size_t, 2> two_opt_cuts{};
    bool physical_relocate{};
    std::size_t relocate_source_route{}, relocate_target_route{};
    std::size_t relocate_source_position{}, relocate_target_slot{};
  };

  [[nodiscard]] static bool described_move_less(
      const MoveDescriptor& left, const MoveDescriptor& right,
      const std::vector<std::vector<int>>& current);
  void order_relocate_moves(const std::vector<std::vector<int>>& current);
  void order_two_opt_moves(const std::vector<std::vector<int>>& current);

  struct MoveMemoryKey {
    Neighborhood neighborhood{};
    // Move descriptors store absolute route slots.  A content version can
    // reappear in another slot after inter-route moves, so both parts are
    // required to identify a support whose cached descriptor is applicable.
    std::array<std::size_t, 2> route_indices{};
    std::array<std::uint64_t, 2> route_versions{};
    std::size_t route_count{};
    bool opens_route{};
    bool physical_customer_mode{};
    bool operator==(const MoveMemoryKey&) const = default;
  };

  struct MoveMemoryKeyHash {
    std::size_t operator()(const MoveMemoryKey& key) const noexcept {
      std::size_t seed = static_cast<std::size_t>(key.neighborhood);
      hash_combine(seed, key.opens_route);
      hash_combine(seed, key.route_count);
      hash_combine(seed, key.physical_customer_mode);
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

  // Physical inputs, not scores or winners. Absolute route slots are absent:
  // the two input identities follow the support's local slot order.
  struct PhysicalSupportKey {
    Neighborhood neighborhood{};
    std::array<std::uint64_t, 2> routes{};
    bool opens_route{};
    bool physical_customer_mode{};
    bool operator==(const PhysicalSupportKey&) const = default;
  };
  struct PhysicalSupportKeyHash {
    std::size_t operator()(const PhysicalSupportKey& key) const noexcept {
      std::size_t seed = static_cast<std::size_t>(key.neighborhood);
      hash_combine(seed, key.routes[0]);
      hash_combine(seed, key.routes[1]);
      hash_combine(seed, key.opens_route);
      hash_combine(seed, key.physical_customer_mode);
      return seed;
    }
  };
  // The customer-move cache is consumed only by generalized_cost. Keep its
  // three physical inputs, never a lambda-dependent score or feasibility proof.
  struct PhysicalScoreSummary {
    double raw_cost_h{};
    double energy_shortfall_wh{};
    double duration_excess_h{};

    PhysicalScoreSummary() = default;
    explicit PhysicalScoreSummary(const RouteCostSummary& summary) noexcept
        : raw_cost_h(summary.raw_cost_h),
          energy_shortfall_wh(summary.energy_shortfall_wh),
          duration_excess_h(summary.duration_excess_h) {}
    [[nodiscard]] double generalized_cost(const Instance& instance,
                                          double penalty_lambda) const;
  };
  struct PhysicalSummaryRow {
    std::array<PhysicalScoreSummary, 2> summaries{};
    unsigned char valid{};
  };
  struct PhysicalSummaryTable {
    std::vector<PhysicalSummaryRow> rows;
  };

  struct SupportWorkingMemory {
    std::uint64_t epoch{};
    SupportState state{SupportState::Unknown};
    bool active{};
    MoveMemoryKey key{};
    MoveMemoryEntry entry{};
    // Borrowed only in this epoch; best_customer_neighbor owns the tables
    // until the pass ends and support_memory clears this before another epoch.
    PhysicalSummaryTable* physical_summaries{};
  };

  [[nodiscard]] std::size_t select_active_neighborhood(unsigned int active);
  [[nodiscard]] std::optional<Solution> best_neighbor(
      const Solution& solution, Neighborhood neighborhood,
      double penalty_lambda,
      std::optional<std::chrono::steady_clock::time_point> deadline,
      bool cooperative_physical_customers = false);
  [[nodiscard]] std::optional<Solution> best_customer_neighbor(
      const Solution& solution, Neighborhood neighborhood,
      double penalty_lambda,
      std::optional<std::chrono::steady_clock::time_point> deadline,
      bool cooperative_physical_customers = false);
  [[nodiscard]] std::optional<Solution> replace_path_neighbor(
      const Solution& solution, double penalty_lambda,
      std::optional<std::chrono::steady_clock::time_point> deadline);
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
  std::size_t replace_path_attempts_;
  bool project_opening_station_visits_;
  // Immutable per VND/cache owner; the two physical scoring domains never mix.
  const bool preserve_swap_station_visits_;
  struct SwapPhysicalProfile {
    std::vector<int> visits;
    std::vector<std::size_t> customer_slots;
    std::unique_ptr<PhysicalRouteScorer> scorer;
  };
  std::vector<SwapPhysicalProfile> swap_profiles_;
  LruCache<MoveMemoryKey, MoveMemoryEntry, MoveMemoryKeyHash> move_memory_{
      100'000};
  struct RouteStateVersion {
    std::uint64_t id{};
    // A weak owner prevents address-reuse alias without retaining SOC buffers.
    std::weak_ptr<const RouteEvaluation> evaluation;
  };
  LruCache<std::vector<int>, RouteStateVersion, IntVectorHash>
      route_versions_{8'192};
  std::uint64_t next_route_version_{1};
  bool move_memory_usable_{true};
  bool move_memory_reusable_{};
  std::optional<double> move_memory_lambda_;
  std::size_t move_memory_route_count_{};
  // A separate bounded identity registry survives ILS calls. Eviction assigns
  // a fresh identity on reinsertion, so it can only lose physical cache hits.
  LruCache<std::vector<int>, std::uint64_t, IntVectorHash>
      physical_route_ids_{8'192};
  std::uint64_t next_physical_route_id_{1};
  LruCache<PhysicalSupportKey, std::shared_ptr<PhysicalSummaryTable>,
           PhysicalSupportKeyHash> physical_summaries_{65'536};
  std::size_t physical_summary_bytes_{};
  std::size_t physical_summary_budget_{128 * 1024 * 1024};
  std::unique_ptr<RouteSequenceScorer> sequence_scorer_;
  std::vector<SupportWorkingMemory> support_working_memory_;
  std::uint64_t support_epoch_{};
  std::vector<std::pair<int, int>> swap_customer_pairs_;
  // Rebuilt each customer-neighborhood pass; preserves open route supports.
  std::vector<unsigned char> customer_open_support_flags_;
  std::vector<unsigned char> swap_order_flags_;
  std::vector<std::tuple<int, std::size_t, std::size_t>> relocate_moves_;
  std::vector<
      std::tuple<std::size_t, std::size_t, std::size_t, std::size_t, bool>>
      two_opt_moves_;
  // Per-pass integer ordering workspace; no scoring or decisions survive it.
  std::vector<std::size_t> order_gap_offsets_;
  std::vector<unsigned char> relocate_order_flags_;
  std::vector<unsigned char> two_opt_order_flags_;
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
