#pragma once

#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <type_traits>
#include <unordered_map>
#include <variant>

#include "ils_sp/core.hpp"
#include "ils_sp/lru_cache.hpp"

namespace ils_sp {

struct ChargingSequenceLeg {
  int station_id{};
  double energy_from_previous_station_wh{};
};

// Associative resource label for a physical visit subsequence.  Customers
// between two charging opportunities are collapsed into one energy segment;
// station order is retained because nonlinear charging decisions depend on it.
struct ChargingSequenceLabel {
  bool initialized{};
  int origin_id{};
  int target_id{};
  std::vector<ChargingSequenceLeg> station_legs;
  double tail_energy_wh{};
  double travel_time_h{};
  double service_time_h{};
};

struct RouteCostSummary {
  bool feasible{};
  double raw_cost_h{};
  double travel_time_h{};
  double charging_time_h{};
  double service_time_h{};
  double route_duration_h{};
  double energy_shortfall_wh{};
  double duration_excess_h{};
};

[[nodiscard]] double generalized_cost(const RouteCostSummary& summary,
                                      const Instance& instance,
                                      double penalty_lambda);

// A physical lower bound on Algorithm 3's result for the same complete input.
// It is neither a completed evaluation nor a feasibility/optimality result.
struct NonlinearCostLowerBound {
  double travel_time_h{};
  double energy_shortfall_wh{};
  double duration_excess_h{};

  [[nodiscard]] double generalized_cost(const Instance& instance,
                                        double penalty_lambda) const;
};

// A half-open customer interval in one incumbent route.  Reversed slices are
// traversed from end - 1 down to begin.  Customer-neighborhood descriptors use
// a constant number of these slices, so they can be scored without first
// materializing and rescanning a complete candidate customer sequence.
struct IncumbentRouteSlice {
  std::size_t route_index{};
  std::size_t begin{};
  std::size_t end{};
  bool reversed{};
  bool operator==(const IncumbentRouteSlice&) const = default;
};

class NonlinearSequenceEvaluator {
 public:
  explicit NonlinearSequenceEvaluator(const Instance& instance)
      : instance_(instance) {}

  [[nodiscard]] ChargingSequenceLabel edge_label(
      int origin_id, std::span<const int> station_ids, int target_id) const;
  [[nodiscard]] ChargingSequenceLabel route_label(const Route& route) const;
  void append(ChargingSequenceLabel& left,
              const ChargingSequenceLabel& right) const;
  [[nodiscard]] ChargingSequenceLabel concatenate(
      ChargingSequenceLabel left,
      const ChargingSequenceLabel& right) const;
  [[nodiscard]] RouteCostSummary evaluate(
      const ChargingSequenceLabel& route) const;
  [[nodiscard]] RouteCostSummary evaluate(
      std::span<const ChargingSequenceLabel* const> components) const;
  [[nodiscard]] double generalized_cost(
      const ChargingSequenceLabel& route, double penalty_lambda) const;
  [[nodiscard]] double generalized_cost(
      std::span<const ChargingSequenceLabel* const> components,
      double penalty_lambda) const;

 private:
  friend class PathSampler;
  friend struct NonlinearSequenceEvaluatorTestAccess;

  using BoundedEvaluation =
      std::variant<RouteCostSummary, NonlinearCostLowerBound>;
  struct BoundedLimits {
    double penalty_lambda{};
    double incumbent_limit{};
    double minimum_limit{};
    bool all_charging_rates_finite{};
  };

  // The caller certifies the complete physical label, immutable Instance data
  // and runtime FP scope. Disable builds with a different arithmetic kernel.
  [[nodiscard]] static constexpr bool bounded_evaluation_compatible() noexcept {
#if defined(__x86_64__) && defined(__SSE2__) && !defined(__FMA__) && \
    !defined(__FMA4__) && !defined(__FAST_MATH__) && \
    (!defined(__FINITE_MATH_ONLY__) || __FINITE_MATH_ONLY__ == 0) && \
    (!defined(__FLT_EVAL_METHOD__) || __FLT_EVAL_METHOD__ == 0)
    return sizeof(double) == 8 && std::numeric_limits<double>::is_iec559 &&
           std::numeric_limits<double>::digits == 53;
#else
    return false;
#endif
  }

  [[nodiscard]] BoundedEvaluation evaluate_bounded(
      const ChargingSequenceLabel& route, double penalty_lambda,
      double incumbent_limit, double minimum_limit,
      bool all_charging_rates_finite = false) const;
  [[nodiscard]] const ChargingCurve& station_curve(int station_id) const;
  template <bool Bounded>
  [[nodiscard]] std::conditional_t<Bounded, BoundedEvaluation, RouteCostSummary>
  evaluate_data(
      double travel_time_h, double service_time_h,
      std::span<const ChargingSequenceLeg> station_legs,
      double tail_energy_wh,
      std::conditional_t<Bounded, BoundedLimits, std::monostate> limits = {}) const;

  const Instance& instance_;
  // Only immutable identities, bound lazily on the first actual CS evaluation.
  // Each evaluator belongs to one serial worker and one immutable Instance.
  mutable std::vector<const ChargingCurve*> station_curves_by_id_;
  mutable bool station_curves_ready_{};
};

class RouteEvaluator {
 public:
  explicit RouteEvaluator(const Instance& instance,
                          std::size_t cache_capacity = 200'000);

  [[nodiscard]] std::shared_ptr<const RouteEvaluation> evaluate(
      const Route& route);
  [[nodiscard]] std::shared_ptr<const RouteEvaluation> evaluate_schedule(
      const Route& route,
      const std::vector<std::optional<double>>& station_departure_wh_by_visit);

  [[nodiscard]] std::size_t cache_size() const noexcept;
  [[nodiscard]] std::size_t cache_hits() const noexcept;
  [[nodiscard]] std::size_t cache_misses() const noexcept;
  [[nodiscard]] std::size_t cache_evictions() const noexcept;
  [[nodiscard]] const Instance& instance() const noexcept { return instance_; }

 private:
  void validate_route(const Route& route) const;
  [[nodiscard]] RouteEvaluation evaluate_uncached(const Route& route) const;

  const Instance& instance_;
  LruCache<std::vector<int>, std::shared_ptr<const RouteEvaluation>,
           IntVectorHash>
      cache_;
};

struct AnchorPair {
  int origin{};
  int target{};
  bool operator==(const AnchorPair&) const = default;
};

struct AnchorPairHash {
  std::size_t operator()(const AnchorPair& pair) const noexcept {
    std::size_t seed = 0;
    hash_combine(seed, pair.origin);
    hash_combine(seed, pair.target);
    return seed;
  }
};

class PlanFactory {
 public:
  using GapMap =
      std::unordered_map<AnchorPair, std::vector<int>, AnchorPairHash>;

  PlanFactory(const Instance& instance, RouteEvaluator& evaluator,
              std::size_t exact_cache_capacity = 20'000,
              bool allow_inherited_station_revisits = false);

  [[nodiscard]] Plan make_plan(
      std::vector<int> customer_ids,
      std::vector<std::vector<int>> gaps);
  // Preserve physical station occurrences, including repeats created when
  // customer removal joins two formerly distinct charging paths.
  [[nodiscard]] Plan plan_visits(std::vector<int> visits);
  [[nodiscard]] Plan plan_sequence(const std::vector<int>& customer_ids,
                                   const Solution* reference = nullptr);
  [[nodiscard]] Plan plan_sequence(const std::vector<int>& customer_ids,
                                   const GapMap& reference_gaps);
  [[nodiscard]] std::optional<Plan> exact_plan(
      const std::vector<int>& customer_ids);
  void publish_exact(Plan plan);

  [[nodiscard]] GapMap gap_map(const Solution& reference) const;
  [[nodiscard]] std::size_t exact_cache_size() const noexcept {
    return exact_cache_.size();
  }
  [[nodiscard]] std::size_t exact_cache_hits() const noexcept {
    return exact_cache_.hits();
  }
  [[nodiscard]] std::size_t exact_cache_misses() const noexcept {
    return exact_cache_.misses();
  }
  [[nodiscard]] std::size_t exact_cache_evictions() const noexcept {
    return exact_cache_.evictions();
  }

 private:
  void validate_node_sequence(std::span<const int> ids, NodeKind kind,
                              const char* error);
  [[nodiscard]] Plan assemble_plan(std::vector<int> customer_ids,
                                  std::vector<std::vector<int>> gaps,
                                  Route route);

  const Instance& instance_;
  RouteEvaluator& evaluator_;
  LruCache<std::vector<int>, Plan, IntVectorHash> exact_cache_;
  bool allow_inherited_station_revisits_{};
  std::vector<std::uint64_t> validation_marks_;
  std::uint64_t validation_epoch_{};
};

// Fixed physical-route labels for the explicit customer-only LNS edit mode.
// Slots and removal positions index Route::visits, not customer-only gaps.
// Candidate scoring never changes, deduplicates, or completes station paths.
class PhysicalRouteScorer {
 public:
  PhysicalRouteScorer(const Instance& instance, const Plan& plan,
                      bool support_swaps = false);
  [[nodiscard]] RouteCostSummary evaluate_replacement(
      std::size_t physical_customer_position, int customer_id) const;
  [[nodiscard]] RouteCostSummary evaluate_swap(
      std::size_t first_physical_position,
      std::size_t second_physical_position) const;

  // Remove one customer and reinsert it at an ORIGINAL-route physical slot.
  // All other visits stay in order. No route-deletion shortcut applies here.
  [[nodiscard]] RouteCostSummary evaluate_relocate(
      std::size_t physical_customer_position, std::size_t original_physical_slot) const;
  [[nodiscard]] RouteCostSummary evaluate_insertion(
      std::size_t physical_slot, int customer_id) const;
  [[nodiscard]] RouteCostSummary evaluate_removal(
      std::size_t physical_customer_position) const;

 private:
  const Instance& instance_;
  NonlinearSequenceEvaluator evaluator_;
  // Includes the implicit departure and return depot at the two ends.
  std::vector<int> nodes_;
  std::unordered_set<int> customers_;
  std::vector<ChargingSequenceLabel> prefixes_;
  std::vector<ChargingSequenceLabel> suffixes_;
  // Optional immutable range labels for the unchanged interior of an intra
  // Swap. Existing LNS insertion/removal callers do not build this tree.
  std::size_t swap_range_base_{};
  std::vector<ChargingSequenceLabel> swap_ranges_;
};

// Vidal-style incumbent-subsequence preprocessing adapted to the Montoya
// nonlinear charging evaluator.  Consecutive customer fragments (and their
// reversal) are indexed in segment trees.  A candidate route is evaluated by
// concatenating O(log n) stored labels per changed fragment and then replaying
// only its ordered charging opportunities.
class RouteSequenceScorer {
 public:
  RouteSequenceScorer(const Instance& instance, const Solution& reference,
                      const PlanFactory::GapMap& reference_gaps);
  ~RouteSequenceScorer();
  RouteSequenceScorer(RouteSequenceScorer&&) noexcept;
  RouteSequenceScorer& operator=(RouteSequenceScorer&&) noexcept;
  RouteSequenceScorer(const RouteSequenceScorer&) = delete;
  RouteSequenceScorer& operator=(const RouteSequenceScorer&) = delete;

  // Rebind the incumbent while retaining preprocessed labels for physical
  // routes that have not changed.
  void reset(const Solution& reference,
             const PlanFactory::GapMap& reference_gaps);

  [[nodiscard]] RouteCostSummary evaluate(
      const std::vector<int>& customer_sequence) const;
  // Score insertion into one incumbent route without copying and rescanning
  // its complete customer sequence.  The label append order is identical to
  // evaluate(materialized_sequence), including inherited directed gaps.
  [[nodiscard]] RouteCostSummary evaluate_insertion(
      std::size_t route_index, std::size_t position, int customer_id) const;
  [[nodiscard]] double generalized_cost(
      const std::vector<int>& customer_sequence,
      double penalty_lambda) const;
  [[nodiscard]] RouteCostSummary evaluate(
      std::span<const IncumbentRouteSlice> slices) const;
  [[nodiscard]] double generalized_cost(
      std::span<const IncumbentRouteSlice> slices,
      double penalty_lambda) const;

 private:
  struct Implementation;
  std::unique_ptr<Implementation> implementation_;
};

}  // namespace ils_sp
