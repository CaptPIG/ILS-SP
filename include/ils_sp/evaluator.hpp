#pragma once

#include <optional>
#include <span>
#include <unordered_map>

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
  [[nodiscard]] RouteCostSummary evaluate_data(
      double travel_time_h, double service_time_h,
      std::span<const ChargingSequenceLeg> station_legs,
      double tail_energy_wh) const;

  const Instance& instance_;
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
              std::size_t exact_cache_capacity = 20'000);

  [[nodiscard]] Plan make_plan(
      std::vector<int> customer_ids,
      std::vector<std::vector<int>> gaps);
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
  const Instance& instance_;
  RouteEvaluator& evaluator_;
  LruCache<std::vector<int>, Plan, IntVectorHash> exact_cache_;
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
