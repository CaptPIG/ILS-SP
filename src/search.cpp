#include "ils_sp/search.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>
#include <set>
#include <stdexcept>
#include <tuple>

namespace ils_sp {
namespace {

constexpr double kLinkedBoundSafety = 100.0 * kCostTolerance;

template <typename Container>
std::size_t random_index(const Container& values, std::mt19937_64& random) {
  if (values.empty()) {
    throw std::invalid_argument("cannot select from an empty collection");
  }
  return std::uniform_int_distribution<std::size_t>(0, values.size() - 1)(
      random);
}

[[nodiscard]] bool deadline_reached(
    std::optional<std::chrono::steady_clock::time_point> deadline) {
  return deadline.has_value() && std::chrono::steady_clock::now() >= *deadline;
}

[[nodiscard]] double elapsed_seconds(
    std::chrono::steady_clock::time_point started) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                       started)
      .count();
}

struct InsertionCompletionOutcome {
  Plan completed_plan;
  double absolute_cost{};
  bool changed{};
  bool infeasible_to_feasible{};
};

struct InsertionMove {
  int customer_id{};
  std::size_t route_index{};
  std::size_t position{};
  double incumbent_route_cost{};
  double direct_delta{};
  double linked_lower_bound_delta{};
  double delta{};
  bool evaluated{};
  std::shared_ptr<const InsertionCompletionOutcome> completion;
  std::shared_ptr<const VirtualInsertionCompletion> virtual_completion;
  bool infeasible_to_feasible{};
};

struct RejectedInsertionDescriptor {
  int customer_id{};
  std::size_t route_index{};
  std::size_t position{};
};

struct InsertionCompletionKey {
  std::vector<int> source_physical_visits;
  int customer_id{};
  std::size_t position{};
  bool opens_route{};
  bool operator==(const InsertionCompletionKey&) const = default;
};

struct InsertionCompletionKeyHash {
  std::size_t operator()(const InsertionCompletionKey& key) const noexcept {
    std::size_t seed = IntVectorHash{}(key.source_physical_visits);
    hash_combine(seed, key.customer_id);
    hash_combine(seed, key.position);
    hash_combine(seed, key.opens_route);
    return seed;
  }
};

using InsertionCompletionMemo =
    LruCache<InsertionCompletionKey,
             std::shared_ptr<const InsertionCompletionOutcome>,
             InsertionCompletionKeyHash>;

constexpr std::size_t kInsertionCompletionMemoCapacity = 8'192;

using VirtualInsertionCompletionMemo =
    LruCache<InsertionCompletionKey,
             std::shared_ptr<const VirtualInsertionCompletion>,
             InsertionCompletionKeyHash>;

constexpr std::size_t kVirtualInsertionCompletionMemoCapacity = 8'192;

[[nodiscard]] bool insertion_better(const InsertionMove& left,
                                    const InsertionMove& right) {
  return std::tie(left.delta, left.customer_id, left.route_index,
                  left.position) <
         std::tie(right.delta, right.customer_id, right.route_index,
                  right.position);
}

[[nodiscard]] bool direct_insertion_better(const InsertionMove& left,
                                           const InsertionMove& right) {
  return std::tie(left.direct_delta, left.customer_id, left.route_index,
                  left.position) <
         std::tie(right.direct_delta, right.customer_id, right.route_index,
                  right.position);
}

[[nodiscard]] bool same_insertion(const InsertionMove& left,
                                  const InsertionMove& right) {
  return std::tie(left.customer_id, left.route_index, left.position) ==
         std::tie(right.customer_id, right.route_index, right.position);
}

[[nodiscard]] bool completes_construction_winner(SearchProfile profile) {
  return profile == SearchProfile::LinkedWinner ||
         profile == SearchProfile::LinkedInsertionCandidates;
}

[[nodiscard]] bool completes_vnd_winner(SearchProfile profile) {
  return profile == SearchProfile::LinkedWinner ||
         profile == SearchProfile::LinkedInsertionCandidates;
}

[[nodiscard]] bool completes_insertion_candidates(SearchProfile profile) {
  return profile == SearchProfile::LinkedInsertionCandidates;
}

[[nodiscard]] bool uses_granular_topology(SearchProfile profile) {
  return profile == SearchProfile::GranularAdaptive;
}

[[nodiscard]] bool uses_adaptive_vnd(SearchProfile profile) {
  return profile == SearchProfile::GranularAdaptive;
}

[[nodiscard]] double maximum_charging_rate(const Instance& instance) {
  double result = 0.0;
  for (const auto& [type, curve] : instance.charging_curves()) {
    (void)type;
    for (std::size_t index = 1; index < curve.points().size(); ++index) {
      const ChargingPoint& left = curve.points()[index - 1];
      const ChargingPoint& right = curve.points()[index];
      const double energy = right.energy_wh - left.energy_wh;
      const double time = right.time_h - left.time_h;
      if (energy <= 0.0) continue;
      if (time <= 1e-15) return std::numeric_limits<double>::infinity();
      const long double rate =
          static_cast<long double>(energy) / static_cast<long double>(time);
      const double rounded_rate = static_cast<double>(rate);
      result = std::max(
          result,
          std::nextafter(rounded_rate,
                         std::numeric_limits<double>::infinity()));
    }
  }
  return result;
}

// Optimistic relaxation for a linked insertion candidate.  It keeps the
// direct/inherited travel and service, allows charging anywhere at the global
// fastest marginal rate, and ignores capacity/reachability constraints.  It
// therefore cannot exceed the cost returned by repair_new_gaps().
[[nodiscard]] double linked_route_lower_bound(
    const RouteCostSummary& summary, const Instance& instance,
    double penalty_lambda, double maximum_rate_wh_per_h) {
  const long double travel = summary.travel_time_h;
  const long double service = summary.service_time_h;
  const long double speed = instance.vehicle().speed_kmph;
  const long double consumption =
      instance.vehicle().consumption_wh_per_km;
  const long double capacity = instance.vehicle().battery_capacity_wh;
  const long double maximum_duration =
      instance.vehicle().max_route_duration_h;
  const long double lambda = penalty_lambda;
  const long double required_energy = std::max(
      0.0L, travel * speed * consumption - capacity);

  const auto relaxed_cost = [&](long double charged_energy) {
    long double charging_time = 0.0L;
    if (charged_energy > 0.0L) {
      if (maximum_rate_wh_per_h <= 0.0) {
        return std::numeric_limits<long double>::infinity();
      }
      if (std::isfinite(maximum_rate_wh_per_h)) {
        charging_time = charged_energy / maximum_rate_wh_per_h;
      }
    }
    const long double shortfall =
        std::max(0.0L, required_energy - charged_energy) / consumption;
    const long double duration_excess = std::max(
        0.0L, travel + service + charging_time - maximum_duration);
    return travel + charging_time + lambda * (shortfall + duration_excess);
  };

  long double lower_bound = relaxed_cost(0.0L);
  lower_bound = std::min(lower_bound, relaxed_cost(required_energy));
  if (std::isfinite(maximum_rate_wh_per_h) &&
      maximum_rate_wh_per_h > 0.0) {
    const long double duration_kink = std::clamp(
        static_cast<long double>(maximum_rate_wh_per_h) *
            (maximum_duration - travel - service),
        0.0L, required_energy);
    lower_bound = std::min(lower_bound, relaxed_cost(duration_kink));
  }
  return std::nextafter(static_cast<double>(lower_bound),
                        -std::numeric_limits<double>::infinity());
}

[[nodiscard]] std::vector<int> insertion_sequence(
    const Solution& current, const InsertionMove& move) {
  if (move.route_index == current.plans.size()) {
    return {move.customer_id};
  }
  std::vector<int> sequence = current.plans[move.route_index].customer_ids;
  sequence.insert(sequence.begin() +
                      static_cast<std::ptrdiff_t>(move.position),
                  move.customer_id);
  return sequence;
}

[[nodiscard]] double sequence_cost(
    const std::vector<int>& sequence, const RouteSequenceScorer& scorer,
    double penalty_lambda) {
  return scorer.generalized_cost(sequence, penalty_lambda);
}

[[nodiscard]] std::vector<InsertionMove> insertion_moves(
    const Solution& current, int customer_id, const Instance& instance,
    const RouteSequenceScorer& scorer, double penalty_lambda,
    bool calculate_linked_lower_bound, double maximum_rate_wh_per_h,
    bool rank_by_linked_lower_bound = false,
    const GranularNeighborhood* granular = nullptr,
    std::uint64_t* granular_rejections = nullptr,
    std::vector<RejectedInsertionDescriptor>* rejected = nullptr) {
  std::vector<InsertionMove> result;
  result.reserve(instance.customer_ids().size() + current.plans.size() + 1);
  for (std::size_t route_index = 0; route_index < current.plans.size();
       ++route_index) {
    const Plan& old_plan = current.plans[route_index];
    const double old_cost =
        generalized_cost(*old_plan.evaluation, instance, penalty_lambda);
    for (std::size_t position = 0; position <= old_plan.customer_ids.size();
         ++position) {
      if (granular != nullptr) {
        if (!granular->allows_insertion(old_plan.customer_ids, position,
                                        customer_id)) {
          if (granular_rejections != nullptr) ++*granular_rejections;
          if (rejected != nullptr) {
            rejected->push_back(RejectedInsertionDescriptor{
                .customer_id = customer_id,
                .route_index = route_index,
                .position = position});
          }
          continue;
        }
      }
      const RouteCostSummary summary =
          scorer.evaluate_insertion(route_index, position, customer_id);
      const double delta =
          generalized_cost(summary, instance, penalty_lambda) - old_cost;
      double lower_bound = delta;
      if (calculate_linked_lower_bound) {
        lower_bound = std::nextafter(
            linked_route_lower_bound(summary, instance, penalty_lambda,
                                     maximum_rate_wh_per_h) -
                old_cost,
            -std::numeric_limits<double>::infinity());
      }
      result.push_back(InsertionMove{.customer_id = customer_id,
                                     .route_index = route_index,
                                     .position = position,
                                     .incumbent_route_cost = old_cost,
                                     .direct_delta = delta,
                                     .linked_lower_bound_delta = lower_bound,
                                     .delta = rank_by_linked_lower_bound
                                                  ? lower_bound
                                                  : delta,
                                     .evaluated = false,
                                     .completion = nullptr,
                                     .virtual_completion = nullptr,
                                     .infeasible_to_feasible = false});
    }
  }
  const std::vector<int> singleton{customer_id};
  double singleton_cost{};
  double singleton_lower_bound{};
  if (calculate_linked_lower_bound) {
    const RouteCostSummary singleton_summary = scorer.evaluate(singleton);
    singleton_cost =
        generalized_cost(singleton_summary, instance, penalty_lambda);
    singleton_lower_bound = linked_route_lower_bound(
        singleton_summary, instance, penalty_lambda,
        maximum_rate_wh_per_h);
  } else {
    singleton_cost = sequence_cost(singleton, scorer, penalty_lambda);
    singleton_lower_bound = singleton_cost;
  }
  // Opening a route is the completeness fallback for a repair state in which
  // none of the remaining customers is currently adjacent to this customer's
  // directed Gamma set.  It has no customer-customer arc to test, so it is
  // retained explicitly rather than silently disabling route creation.
  result.push_back(InsertionMove{
      .customer_id = customer_id,
      .route_index = current.plans.size(),
      .position = 0,
      .incumbent_route_cost = 0.0,
      .direct_delta = singleton_cost,
      .linked_lower_bound_delta = singleton_lower_bound,
      .delta = rank_by_linked_lower_bound ? singleton_lower_bound
                                          : singleton_cost,
      .evaluated = false,
      .completion = nullptr,
      .virtual_completion = nullptr,
      .infeasible_to_feasible = false});

  return result;
}

[[nodiscard]] InsertionMove rejected_insertion_move(
    const Solution& current, const RejectedInsertionDescriptor& descriptor,
    const Instance& instance, const RouteSequenceScorer& scorer,
    double penalty_lambda, double maximum_rate_wh_per_h) {
  if (descriptor.route_index >= current.plans.size()) {
    throw std::out_of_range("granular escape has invalid route");
  }
  const Plan& old_plan = current.plans[descriptor.route_index];
  const double old_cost =
      generalized_cost(*old_plan.evaluation, instance, penalty_lambda);
  const RouteCostSummary summary = scorer.evaluate_insertion(
      descriptor.route_index, descriptor.position, descriptor.customer_id);
  const double direct_cost =
      generalized_cost(summary, instance, penalty_lambda);
  const double lower_bound = std::nextafter(
      linked_route_lower_bound(summary, instance, penalty_lambda,
                               maximum_rate_wh_per_h) -
          old_cost,
      -std::numeric_limits<double>::infinity());
  return InsertionMove{.customer_id = descriptor.customer_id,
                       .route_index = descriptor.route_index,
                       .position = descriptor.position,
                       .incumbent_route_cost = old_cost,
                       .direct_delta = direct_cost - old_cost,
                       .linked_lower_bound_delta = lower_bound,
                       .delta = direct_cost - old_cost,
                       .evaluated = false,
                       .completion = nullptr,
                       .virtual_completion = nullptr,
                       .infeasible_to_feasible = false};
}

void append_granular_escape(
    std::vector<InsertionMove>& candidates,
    const std::vector<RejectedInsertionDescriptor>& rejected,
    const Solution& current, const Instance& instance,
    const RouteSequenceScorer& scorer, double penalty_lambda,
    double maximum_rate_wh_per_h, std::mt19937_64& random,
    InsertionCandidateStatistics& statistics) {
  statistics.granular_escape_candidates_seen += rejected.size();
  if (rejected.empty()) return;
  const RejectedInsertionDescriptor& selected =
      rejected[random_index(rejected, random)];
  candidates.push_back(rejected_insertion_move(
      current, selected, instance, scorer, penalty_lambda,
      maximum_rate_wh_per_h));
  ++statistics.granular_escape_candidates_selected;
}

void evaluate_virtual_insertion(
    InsertionMove& move, const Solution& current, PathSampler& path_sampler,
    double penalty_lambda, InsertionCandidateStatistics& statistics,
    VirtualInsertionCompletionMemo& completion_memo,
    RepairTimingStatistics* timing) {
  if (move.virtual_completion != nullptr) return;
  InsertionCompletionKey memo_key{
      .source_physical_visits =
          move.route_index == current.plans.size()
              ? std::vector<int>{}
              : current.plans[move.route_index].route.visits,
      .customer_id = move.customer_id,
      .position = move.position,
      .opens_route = move.route_index == current.plans.size()};
  std::shared_ptr<const VirtualInsertionCompletion> outcome;
  bool scored_options = false;
  if (auto cached = completion_memo.get(memo_key); cached.has_value()) {
    ++statistics.virtual_completion_memo_hits;
    outcome = std::move(*cached);
  } else {
    ++statistics.virtual_completion_memo_misses;
    const Plan* source = move.route_index == current.plans.size()
                             ? nullptr
                             : &current.plans[move.route_index];
    outcome = std::make_shared<const VirtualInsertionCompletion>(
        path_sampler.virtual_complete_insertion(
            source, move.position, move.customer_id, penalty_lambda, timing));
    scored_options = true;
    completion_memo.put(std::move(memo_key), outcome);
  }
  move.delta = outcome->generalized_cost_h - move.incumbent_route_cost;
  if (move.linked_lower_bound_delta > move.delta + kLinkedBoundSafety) {
    throw std::logic_error(
        "virtual insertion lower bound exceeded its completed cost");
  }
  move.infeasible_to_feasible = outcome->infeasible_to_feasible;
  move.virtual_completion = std::move(outcome);
  move.evaluated = true;
  ++statistics.virtual_candidates_evaluated;
  if (scored_options) {
    statistics.virtual_options_scored +=
        move.virtual_completion->options_scored;
    statistics.virtual_options_pruned_by_lower_bound +=
        move.virtual_completion->options_pruned_by_lower_bound;
  }
  if (move.virtual_completion->changed) {
    ++statistics.virtual_candidates_changed;
  }
  if (move.infeasible_to_feasible) {
    ++statistics.virtual_infeasible_to_feasible;
  }
}

void evaluate_linked_insertion(
    InsertionMove& move, const Solution& current, const Instance& instance,
    PlanFactory& factory, PathSampler& path_sampler,
    const PlanFactory::GapMap& reference_gaps, double penalty_lambda,
    InsertionCandidateStatistics& statistics,
    InsertionCompletionMemo& completion_memo,
    bool repair_affected_paths = false) {
  if (move.evaluated) return;
  InsertionCompletionKey memo_key{
      .source_physical_visits =
          move.route_index == current.plans.size()
              ? std::vector<int>{}
              : current.plans[move.route_index].route.visits,
      .customer_id = move.customer_id,
      .position = move.position,
      .opens_route = move.route_index == current.plans.size()};
  std::shared_ptr<const InsertionCompletionOutcome> outcome;
  if (auto cached = completion_memo.get(memo_key); cached.has_value()) {
    ++statistics.completion_memo_hits;
    outcome = std::move(*cached);
  } else {
    ++statistics.completion_memo_misses;
    Plan direct = factory.plan_sequence(insertion_sequence(current, move),
                                        reference_gaps);
    const bool direct_feasible = direct.evaluation->feasible;
    const auto direct_gaps = direct.gaps;
    Plan completed = std::move(direct);
    if (repair_affected_paths) {
      completed = path_sampler.repair_affected_paths(
          std::move(completed), reference_gaps, penalty_lambda, std::nullopt,
          &statistics.affected_path_completion);
    } else if (!(completed.evaluation->feasible &&
                 completed.evaluation->charging_time_h <= kTimeToleranceH)) {
      // A feasible zero-charge direct endpoint is already shortest for every
      // introduced gap.  Any station detour has nonnegative travel and
      // charging time.
      completed = path_sampler.repair_new_gaps(
          std::move(completed), reference_gaps, penalty_lambda);
    }
    const double absolute_cost =
        generalized_cost(*completed.evaluation, instance, penalty_lambda);
    const bool changed = completed.gaps != direct_gaps;
    const bool infeasible_to_feasible =
        !direct_feasible && completed.evaluation->feasible;
    outcome = std::make_shared<const InsertionCompletionOutcome>(
        InsertionCompletionOutcome{
            .completed_plan = std::move(completed),
            .absolute_cost = absolute_cost,
            .changed = changed,
            .infeasible_to_feasible = infeasible_to_feasible});
    completion_memo.put(std::move(memo_key), outcome);
  }
  move.delta = outcome->absolute_cost - move.incumbent_route_cost;
  if (!repair_affected_paths && move.linked_lower_bound_delta >
      move.delta + kLinkedBoundSafety) {
    throw std::logic_error("linked insertion lower bound exceeded its cost");
  }
  move.infeasible_to_feasible = outcome->infeasible_to_feasible;
  ++statistics.candidates_evaluated;
  if (outcome->changed) ++statistics.candidates_changed;
  if (move.infeasible_to_feasible) ++statistics.infeasible_to_feasible;
  move.evaluated = true;
  move.completion = std::move(outcome);
}

[[nodiscard]] bool insertion_bound_better(
    const std::vector<InsertionMove>& moves, std::size_t left,
    std::size_t right) {
  const auto left_key =
      std::tie(moves[left].linked_lower_bound_delta, moves[left].customer_id,
               moves[left].route_index, moves[left].position);
  const auto right_key =
      std::tie(moves[right].linked_lower_bound_delta, moves[right].customer_id,
               moves[right].route_index, moves[right].position);
  return left_key < right_key || (left_key == right_key && left < right);
}

// A min-priority view of the immutable linked lower-bound keys.  Constructing
// the heap is linear; only the prefix that survives the exact bound is popped.
class InsertionBoundQueue {
  struct Comparator {
    const std::vector<InsertionMove>* moves{};

    bool operator()(std::size_t left, std::size_t right) const {
      return insertion_bound_better(*moves, right, left);
    }
  };

 public:
  explicit InsertionBoundQueue(const std::vector<InsertionMove>& moves)
      : moves_(moves), heap_(moves.size()) {
    std::iota(heap_.begin(), heap_.end(), 0);
    std::make_heap(heap_.begin(), heap_.end(), comparator());
  }

  [[nodiscard]] bool empty() const noexcept { return heap_.empty(); }

  [[nodiscard]] std::size_t pop() {
    std::pop_heap(heap_.begin(), heap_.end(), comparator());
    const std::size_t result = heap_.back();
    heap_.pop_back();
    return result;
  }

 private:
  [[nodiscard]] Comparator comparator() const {
    return Comparator{&moves_};
  }

  const std::vector<InsertionMove>& moves_;
  std::vector<std::size_t> heap_;
};

struct BestTwoInsertionIndices {
  std::optional<std::size_t> best;
  std::optional<std::size_t> second;
};

template <typename Better>
void consider_best_two(BestTwoInsertionIndices& result, std::size_t index,
                       Better better) {
  if (!result.best.has_value() || better(index, *result.best)) {
    result.second = result.best;
    result.best = index;
  } else if (!result.second.has_value() || better(index, *result.second)) {
    result.second = index;
  }
}

struct DirectCustomerInsertionRanking {
  BestTwoInsertionIndices indices;
  InsertionMove best;
  double regret{};
};

[[nodiscard]] DirectCustomerInsertionRanking direct_customer_ranking(
    const std::vector<InsertionMove>& moves) {
  if (moves.empty()) throw std::logic_error("customer has no insertion");
  BestTwoInsertionIndices indices;
  for (std::size_t index = 0; index < moves.size(); ++index) {
    consider_best_two(indices, index, [&](std::size_t left,
                                          std::size_t right) {
      return direct_insertion_better(moves[left], moves[right]);
    });
  }
  const std::size_t best = *indices.best;
  const double regret =
      indices.second.has_value()
          ? moves[*indices.second].direct_delta - moves[best].direct_delta
          : std::numeric_limits<double>::infinity();
  return DirectCustomerInsertionRanking{
      .indices = indices, .best = moves[best], .regret = regret};
}

[[nodiscard]] InsertionMove choose_linked_greedy(
    std::vector<InsertionMove> moves, const Solution& current,
    const Instance& instance, PlanFactory& factory, PathSampler& path_sampler,
    const PlanFactory::GapMap& reference_gaps, double penalty_lambda,
    InsertionCandidateStatistics& statistics,
    InsertionCompletionMemo& completion_memo) {
  if (moves.empty()) throw std::logic_error("no insertion candidate exists");
  ++statistics.decisions;
  ++statistics.greedy_decisions;
  statistics.candidates_generated += moves.size();

  const auto direct_best = std::min_element(
      moves.begin(), moves.end(), direct_insertion_better);
  const std::size_t direct_index =
      static_cast<std::size_t>(std::distance(moves.begin(), direct_best));
  const InsertionMove direct_identity = moves[direct_index];
  evaluate_linked_insertion(moves[direct_index], current, instance, factory,
                            path_sampler, reference_gaps, penalty_lambda,
                            statistics, completion_memo);
  const double completed_direct_delta = moves[direct_index].delta;
  std::size_t best_index = direct_index;

  InsertionBoundQueue bound_queue(moves);
  while (!bound_queue.empty()) {
    const std::size_t index = bound_queue.pop();
    if (moves[index].evaluated) continue;
    if (moves[index].linked_lower_bound_delta >
        moves[best_index].delta + kLinkedBoundSafety) {
      break;
    }
    evaluate_linked_insertion(moves[index], current, instance, factory,
                              path_sampler, reference_gaps, penalty_lambda,
                              statistics, completion_memo);
    if (insertion_better(moves[index], moves[best_index])) {
      moves[best_index].completion.reset();
      best_index = index;
    } else {
      moves[index].completion.reset();
    }
  }
  statistics.candidates_pruned_by_bound += std::count_if(
      moves.begin(), moves.end(), [](const InsertionMove& move) {
        return !move.evaluated;
      });
  if (!same_insertion(moves[best_index], direct_identity)) {
    ++statistics.winner_flips;
    if (moves[best_index].delta <
        completed_direct_delta - kCostTolerance) {
      ++statistics.strict_winner_flips;
    }
  }
  if (moves[best_index].infeasible_to_feasible) {
    ++statistics.selected_infeasible_to_feasible;
  }
  return std::move(moves[best_index]);
}

struct CustomerInsertionRanking {
  InsertionMove direct_best;
  double direct_regret{};
  double completed_direct_best_delta{};
  InsertionMove linked_best;
  double linked_regret{};
};

[[nodiscard]] std::optional<CustomerInsertionRanking> rank_linked_customer(
    std::vector<InsertionMove> moves,
    const DirectCustomerInsertionRanking& direct, const Solution& current,
    const Instance& instance, PlanFactory& factory, PathSampler& path_sampler,
    const PlanFactory::GapMap& reference_gaps, double penalty_lambda,
    InsertionCandidateStatistics& statistics,
    InsertionCompletionMemo& completion_memo,
    std::optional<double> incumbent_linked_regret,
    bool require_complete_ranking) {
  statistics.candidates_generated += moves.size();
  const std::size_t direct_best_index = *direct.indices.best;
  InsertionMove direct_best = direct.best;

  BestTwoInsertionIndices exact;
  const auto evaluate_and_rank = [&](std::size_t index) {
    evaluate_linked_insertion(moves[index], current, instance, factory,
                              path_sampler, reference_gaps, penalty_lambda,
                              statistics, completion_memo);
    const std::optional<std::size_t> previous_best = exact.best;
    consider_best_two(exact, index, [&](std::size_t left,
                                        std::size_t right) {
      return insertion_better(moves[left], moves[right]);
    });
    if (exact.best != previous_best) {
      if (previous_best.has_value()) {
        moves[*previous_best].completion.reset();
      }
    } else {
      moves[index].completion.reset();
    }
  };

  evaluate_and_rank(direct_best_index);
  if (direct.indices.second.has_value()) {
    evaluate_and_rank(*direct.indices.second);
  }
  const double completed_direct_best_delta = moves[direct_best_index].delta;

  InsertionBoundQueue bound_queue(moves);
  while (!bound_queue.empty()) {
    const std::size_t index = bound_queue.pop();
    if (moves[index].evaluated) continue;
    const double second_upper =
        exact.second.has_value() ? moves[*exact.second].delta
                                 : std::numeric_limits<double>::infinity();
    if (moves[index].linked_lower_bound_delta >
        second_upper + kLinkedBoundSafety) {
      break;
    }
    if (!require_complete_ranking && incumbent_linked_regret.has_value() &&
        exact.best.has_value() && exact.second.has_value()) {
      const double exact_best = moves[*exact.best].delta;
      const double exact_second = moves[*exact.second].delta;
      const double next_lower_bound =
          moves[index].linked_lower_bound_delta;
      if (std::isfinite(*incumbent_linked_regret) &&
          std::isfinite(exact_best) && std::isfinite(exact_second) &&
          std::isfinite(next_lower_bound)) {
        const double outward_next_lower_bound = std::nextafter(
            next_lower_bound, -std::numeric_limits<double>::infinity());
        const long double best_lower_bound = std::min(
            static_cast<long double>(exact_best),
            static_cast<long double>(outward_next_lower_bound) -
                static_cast<long double>(kLinkedBoundSafety));
        const long double upper_long = std::max(
            0.0L, static_cast<long double>(exact_second) - best_lower_bound);
        if (std::isfinite(upper_long) &&
            upper_long <=
                static_cast<long double>(
                    std::numeric_limits<double>::max())) {
          const double regret_upper_bound = std::nextafter(
              static_cast<double>(upper_long),
              std::numeric_limits<double>::infinity());
          if (*incumbent_linked_regret > regret_upper_bound &&
              *incumbent_linked_regret - regret_upper_bound >
                  kCostTolerance) {
            ++statistics.customers_pruned_by_regret;
            statistics.candidates_pruned_by_regret += std::count_if(
                moves.begin(), moves.end(), [](const InsertionMove& move) {
                  return !move.evaluated;
                });
            return std::nullopt;
          }
        }
      }
    }
    evaluate_and_rank(index);
  }
  statistics.candidates_pruned_by_bound += std::count_if(
      moves.begin(), moves.end(), [](const InsertionMove& move) {
        return !move.evaluated;
      });
  const double linked_regret =
      exact.second.has_value()
          ? moves[*exact.second].delta - moves[*exact.best].delta
          : std::numeric_limits<double>::infinity();
  return CustomerInsertionRanking{
      .direct_best = std::move(direct_best),
      .direct_regret = direct.regret,
      .completed_direct_best_delta = completed_direct_best_delta,
      .linked_best = std::move(moves[*exact.best]),
      .linked_regret = linked_regret};
}

[[nodiscard]] bool regret_better(double candidate_regret,
                                 const InsertionMove& candidate,
                                 double incumbent_regret,
                                 const InsertionMove& incumbent) {
  return candidate_regret > incumbent_regret + kCostTolerance ||
         (std::abs(candidate_regret - incumbent_regret) <= kCostTolerance &&
          insertion_better(candidate, incumbent));
}

[[nodiscard]] InsertionMove choose_linked_two_regret(
    const Solution& current, const std::vector<int>& unassigned,
    const Instance& instance, PlanFactory& factory, PathSampler& path_sampler,
    const RouteSequenceScorer& scorer,
    const PlanFactory::GapMap& reference_gaps, double penalty_lambda,
    double maximum_rate_wh_per_h,
    InsertionCandidateStatistics& statistics,
    InsertionCompletionMemo& completion_memo) {
  ++statistics.decisions;
  ++statistics.two_regret_decisions;
  struct PreparedCustomerRanking {
    std::vector<InsertionMove> moves;
    DirectCustomerInsertionRanking direct;
  };
  std::vector<PreparedCustomerRanking> prepared;
  prepared.reserve(unassigned.size());
  for (const int customer_id : unassigned) {
    std::vector<InsertionMove> moves = insertion_moves(
        current, customer_id, instance, scorer, penalty_lambda, true,
        maximum_rate_wh_per_h);
    DirectCustomerInsertionRanking direct =
        direct_customer_ranking(moves);
    prepared.push_back(PreparedCustomerRanking{
        .moves = std::move(moves), .direct = std::move(direct)});
  }
  std::size_t direct_winner_index = 0;
  for (std::size_t index = 1; index < prepared.size(); ++index) {
    if (regret_better(prepared[index].direct.regret,
                      prepared[index].direct.best,
                      prepared[direct_winner_index].direct.regret,
                      prepared[direct_winner_index].direct.best)) {
      direct_winner_index = index;
    }
  }
  const DirectCustomerInsertionRanking direct_winner =
      prepared[direct_winner_index].direct;
  std::optional<CustomerInsertionRanking> linked_winner;
  double direct_winner_linked_regret{};
  for (std::size_t index = 0; index < prepared.size(); ++index) {
    std::optional<CustomerInsertionRanking> linked = rank_linked_customer(
        std::move(prepared[index].moves), prepared[index].direct, current,
        instance, factory, path_sampler, reference_gaps, penalty_lambda,
        statistics, completion_memo,
        linked_winner.has_value()
            ? std::optional<double>{linked_winner->linked_regret}
            : std::nullopt,
        index == direct_winner_index);
    if (index == direct_winner_index) {
      if (!linked.has_value()) {
        throw std::logic_error(
            "direct 2-Regret winner was incompletely ranked");
      }
      direct_winner_linked_regret = linked->linked_regret;
    }
    if (linked.has_value() &&
        (!linked_winner.has_value() ||
         regret_better(linked->linked_regret, linked->linked_best,
                       linked_winner->linked_regret,
                       linked_winner->linked_best))) {
      linked_winner = std::move(linked);
    }
  }
  if (!linked_winner.has_value()) {
    throw std::logic_error("2-Regret insertion did not select a customer");
  }

  if (!same_insertion(linked_winner->linked_best, direct_winner.best)) {
    ++statistics.winner_flips;
    bool strict = false;
    if (linked_winner->linked_best.customer_id !=
        direct_winner.best.customer_id) {
      strict = linked_winner->linked_regret >
               direct_winner_linked_regret + kCostTolerance;
    } else {
      strict = linked_winner->linked_best.delta <
               linked_winner->completed_direct_best_delta -
                   kCostTolerance;
    }
    if (strict) ++statistics.strict_winner_flips;
  }
  if (linked_winner->linked_best.infeasible_to_feasible) {
    ++statistics.selected_infeasible_to_feasible;
  }
  return std::move(linked_winner->linked_best);
}

[[nodiscard]] InsertionMove choose_lazy_path_greedy(
    std::vector<InsertionMove> moves, const Solution& current,
    PathSampler& path_sampler, double penalty_lambda,
    InsertionCandidateStatistics& statistics,
    VirtualInsertionCompletionMemo& completion_memo,
    RepairTimingStatistics* timing) {
  if (moves.empty()) throw std::logic_error("no insertion candidate exists");
  ++statistics.decisions;
  ++statistics.greedy_decisions;
  statistics.candidates_generated += moves.size();
  statistics.candidates_evaluated += moves.size();
  const auto direct_ranking_started = std::chrono::steady_clock::now();
  const std::size_t direct_best = static_cast<std::size_t>(std::distance(
      moves.begin(),
      std::min_element(moves.begin(), moves.end(), direct_insertion_better)));
  if (timing != nullptr) {
    timing->granular_direct_scoring_s +=
        elapsed_seconds(direct_ranking_started);
  }
  const InsertionMove direct_identity = moves[direct_best];
  const auto bound_started = std::chrono::steady_clock::now();
  InsertionBoundQueue bound_queue(moves);
  if (timing != nullptr) {
    timing->frontier_selection_s += elapsed_seconds(bound_started);
  }
  const auto completion_started = std::chrono::steady_clock::now();
  evaluate_virtual_insertion(moves[direct_best], current, path_sampler,
                             penalty_lambda, statistics, completion_memo,
                             timing);
  const double completed_direct_delta = moves[direct_best].delta;
  std::size_t best = direct_best;
  while (!bound_queue.empty()) {
    const std::size_t index = bound_queue.pop();
    if (moves[index].evaluated) continue;
    if (moves[index].linked_lower_bound_delta >
        moves[best].delta + kLinkedBoundSafety) {
      break;
    }
    evaluate_virtual_insertion(moves[index], current, path_sampler,
                               penalty_lambda, statistics, completion_memo,
                               timing);
    if (insertion_better(moves[index], moves[best])) {
      moves[best].virtual_completion.reset();
      best = index;
    } else {
      moves[index].virtual_completion.reset();
    }
  }
  if (timing != nullptr) {
    timing->virtual_completion_s += elapsed_seconds(completion_started);
  }
  statistics.virtual_candidates_pruned_by_bound += std::count_if(
      moves.begin(), moves.end(),
      [](const InsertionMove& move) { return !move.evaluated; });
  if (!same_insertion(moves[best], direct_identity)) {
    ++statistics.winner_flips;
    ++statistics.virtual_winner_flips;
    if (moves[best].delta < completed_direct_delta - kCostTolerance) {
      ++statistics.strict_winner_flips;
      ++statistics.virtual_strict_winner_flips;
    }
  }
  if (moves[best].infeasible_to_feasible) {
    ++statistics.selected_infeasible_to_feasible;
  }
  return std::move(moves[best]);
}

struct VirtualCustomerRanking {
  InsertionMove direct_best;
  double direct_regret{};
  double completed_direct_best_delta{};
  InsertionMove virtual_best;
  double virtual_regret{};
};

[[nodiscard]] std::optional<VirtualCustomerRanking> rank_virtual_customer(
    std::vector<InsertionMove> moves,
    const DirectCustomerInsertionRanking& direct, const Solution& current,
    PathSampler& path_sampler, double penalty_lambda,
    InsertionCandidateStatistics& statistics,
    VirtualInsertionCompletionMemo& completion_memo,
    std::optional<double> incumbent_virtual_regret,
    bool require_complete_ranking, RepairTimingStatistics* timing) {
  const std::size_t direct_best_index = *direct.indices.best;
  BestTwoInsertionIndices exact;
  const auto evaluate_and_rank = [&](std::size_t index) {
    evaluate_virtual_insertion(moves[index], current, path_sampler,
                               penalty_lambda, statistics, completion_memo,
                               timing);
    const std::optional<std::size_t> previous_best = exact.best;
    consider_best_two(exact, index, [&](std::size_t left,
                                        std::size_t right) {
      return insertion_better(moves[left], moves[right]);
    });
    if (exact.best != previous_best) {
      if (previous_best.has_value()) {
        moves[*previous_best].virtual_completion.reset();
      }
    } else {
      moves[index].virtual_completion.reset();
    }
  };

  evaluate_and_rank(direct_best_index);
  if (direct.indices.second.has_value()) {
    evaluate_and_rank(*direct.indices.second);
  }
  const double completed_direct_best_delta = moves[direct_best_index].delta;
  const auto bound_started = std::chrono::steady_clock::now();
  InsertionBoundQueue bound_queue(moves);
  if (timing != nullptr) {
    timing->frontier_selection_s += elapsed_seconds(bound_started);
  }
  while (!bound_queue.empty()) {
    const std::size_t index = bound_queue.pop();
    if (moves[index].evaluated) continue;
    const double second_upper =
        exact.second.has_value() ? moves[*exact.second].delta
                                 : std::numeric_limits<double>::infinity();
    if (moves[index].linked_lower_bound_delta >
        second_upper + kLinkedBoundSafety) {
      break;
    }
    if (!require_complete_ranking && incumbent_virtual_regret.has_value() &&
        exact.best.has_value() && exact.second.has_value()) {
      const double exact_best = moves[*exact.best].delta;
      const double exact_second = moves[*exact.second].delta;
      const double next_lower_bound = moves[index].linked_lower_bound_delta;
      if (std::isfinite(*incumbent_virtual_regret) &&
          std::isfinite(exact_best) && std::isfinite(exact_second) &&
          std::isfinite(next_lower_bound)) {
        const double outward_next_lower_bound = std::nextafter(
            next_lower_bound, -std::numeric_limits<double>::infinity());
        const long double best_lower_bound = std::min(
            static_cast<long double>(exact_best),
            static_cast<long double>(outward_next_lower_bound) -
                static_cast<long double>(kLinkedBoundSafety));
        const long double upper_long = std::max(
            0.0L,
            static_cast<long double>(exact_second) - best_lower_bound);
        if (std::isfinite(upper_long) &&
            upper_long <=
                static_cast<long double>(
                    std::numeric_limits<double>::max())) {
          const double regret_upper_bound = std::nextafter(
              static_cast<double>(upper_long),
              std::numeric_limits<double>::infinity());
          if (*incumbent_virtual_regret > regret_upper_bound &&
              *incumbent_virtual_regret - regret_upper_bound >
                  kCostTolerance) {
            ++statistics.virtual_customers_pruned_by_regret;
            statistics.virtual_candidates_pruned_by_regret += std::count_if(
                moves.begin(), moves.end(),
                [](const InsertionMove& move) { return !move.evaluated; });
            return std::nullopt;
          }
        }
      }
    }
    evaluate_and_rank(index);
  }
  statistics.virtual_candidates_pruned_by_bound += std::count_if(
      moves.begin(), moves.end(),
      [](const InsertionMove& move) { return !move.evaluated; });
  const double virtual_regret =
      exact.second.has_value()
          ? moves[*exact.second].delta - moves[*exact.best].delta
          : std::numeric_limits<double>::infinity();
  return VirtualCustomerRanking{
      .direct_best = direct.best,
      .direct_regret = direct.regret,
      .completed_direct_best_delta = completed_direct_best_delta,
      .virtual_best = std::move(moves[*exact.best]),
      .virtual_regret = virtual_regret};
}

[[nodiscard]] InsertionMove choose_lazy_path_two_regret(
    std::vector<InsertionMove> moves, const Solution& current,
    const std::vector<int>& unassigned, PathSampler& path_sampler,
    double penalty_lambda,
    InsertionCandidateStatistics& statistics,
    VirtualInsertionCompletionMemo& completion_memo,
    RepairTimingStatistics* timing) {
  if (moves.empty()) throw std::logic_error("no insertion candidate exists");
  ++statistics.decisions;
  ++statistics.two_regret_decisions;
  statistics.candidates_generated += moves.size();
  statistics.candidates_evaluated += moves.size();
  const auto direct_ranking_started = std::chrono::steady_clock::now();
  struct PreparedCustomerRanking {
    std::vector<InsertionMove> moves;
    DirectCustomerInsertionRanking direct;
  };
  std::unordered_map<int, std::size_t> customer_index;
  customer_index.reserve(unassigned.size());
  std::vector<PreparedCustomerRanking> prepared(unassigned.size());
  for (std::size_t index = 0; index < unassigned.size(); ++index) {
    customer_index.emplace(unassigned[index], index);
  }
  for (InsertionMove& move : moves) {
    prepared[customer_index.at(move.customer_id)].moves.push_back(
        std::move(move));
  }
  for (PreparedCustomerRanking& ranking : prepared) {
    ranking.direct = direct_customer_ranking(ranking.moves);
  }
  std::size_t direct_winner_index = 0;
  for (std::size_t index = 1; index < prepared.size(); ++index) {
    if (regret_better(prepared[index].direct.regret,
                      prepared[index].direct.best,
                      prepared[direct_winner_index].direct.regret,
                      prepared[direct_winner_index].direct.best)) {
      direct_winner_index = index;
    }
  }
  const DirectCustomerInsertionRanking direct_winner =
      prepared[direct_winner_index].direct;
  if (timing != nullptr) {
    timing->granular_direct_scoring_s +=
        elapsed_seconds(direct_ranking_started);
  }
  const auto completion_started = std::chrono::steady_clock::now();
  std::optional<VirtualCustomerRanking> virtual_winner;
  double direct_winner_virtual_regret{};
  for (std::size_t index = 0; index < prepared.size(); ++index) {
    std::optional<VirtualCustomerRanking> ranked = rank_virtual_customer(
        std::move(prepared[index].moves), prepared[index].direct, current,
        path_sampler, penalty_lambda, statistics, completion_memo,
        virtual_winner.has_value()
            ? std::optional<double>{virtual_winner->virtual_regret}
            : std::nullopt,
        index == direct_winner_index, timing);
    if (index == direct_winner_index) {
      if (!ranked.has_value()) {
        throw std::logic_error(
            "direct 2-Regret winner was incompletely virtual-ranked");
      }
      direct_winner_virtual_regret = ranked->virtual_regret;
    }
    if (ranked.has_value() &&
        (!virtual_winner.has_value() ||
         regret_better(ranked->virtual_regret, ranked->virtual_best,
                       virtual_winner->virtual_regret,
                       virtual_winner->virtual_best))) {
      virtual_winner = std::move(ranked);
    }
  }
  if (timing != nullptr) {
    timing->virtual_completion_s += elapsed_seconds(completion_started);
  }
  if (!virtual_winner.has_value()) {
    throw std::logic_error("2-Regret virtual insertion selected no customer");
  }
  if (!same_insertion(virtual_winner->virtual_best, direct_winner.best)) {
    ++statistics.winner_flips;
    ++statistics.virtual_winner_flips;
    const bool strict =
        virtual_winner->virtual_best.customer_id !=
                direct_winner.best.customer_id
            ? virtual_winner->virtual_regret >
                  direct_winner_virtual_regret + kCostTolerance
            : virtual_winner->virtual_best.delta <
                  virtual_winner->completed_direct_best_delta -
                      kCostTolerance;
    if (strict) {
      ++statistics.strict_winner_flips;
      ++statistics.virtual_strict_winner_flips;
    }
  }
  if (virtual_winner->virtual_best.infeasible_to_feasible) {
    ++statistics.selected_infeasible_to_feasible;
  }
  return std::move(virtual_winner->virtual_best);
}

void apply_insertion(Solution& current, std::size_t route_index, Plan plan) {
  if (route_index == current.plans.size()) {
    current.plans.push_back(std::move(plan));
  } else {
    current.plans[route_index] = std::move(plan);
  }
}

void update_insertion_hints(std::vector<InsertionPathHint>& hints,
                            const Solution& current,
                            const InsertionMove& move,
                            const Instance& instance) {
  const auto erase_pair = [&](const AnchorPair& pair) {
    std::erase_if(hints, [&](const InsertionPathHint& hint) {
      return hint.gap == pair;
    });
  };
  if (move.route_index < current.plans.size()) {
    const auto& customers = current.plans[move.route_index].customer_ids;
    const int previous = move.position == 0
                             ? instance.depot().id
                             : customers[move.position - 1];
    const int next = move.position == customers.size()
                         ? instance.depot().id
                         : customers[move.position];
    erase_pair(AnchorPair{previous, next});
  }
  if (move.virtual_completion == nullptr) return;
  for (std::size_t side = 0; side < 2; ++side) {
    erase_pair(move.virtual_completion->gaps[side]);
    if (!move.virtual_completion->paths[side].empty()) {
      hints.push_back(InsertionPathHint{
          .gap = move.virtual_completion->gaps[side],
          .stations = move.virtual_completion->paths[side]});
    }
  }
}

[[nodiscard]] std::vector<int> all_customers(const Solution& solution) {
  std::vector<int> result;
  for (const Plan& plan : solution.plans) {
    result.insert(result.end(), plan.customer_ids.begin(), plan.customer_ids.end());
  }
  return result;
}

void accumulate_dirty_direct_anchor_pairs(
    const PlanFactory::GapMap& previous_gaps, const Solution& current,
    const Instance& instance,
    std::unordered_set<AnchorPair, AnchorPairHash>& dirty) {
  for (const Plan& plan : current.plans) {
    std::vector<int> anchors;
    anchors.reserve(plan.customer_ids.size() + 2);
    anchors.push_back(instance.depot().id);
    anchors.insert(anchors.end(), plan.customer_ids.begin(),
                   plan.customer_ids.end());
    anchors.push_back(instance.depot().id);
    for (std::size_t gap = 0; gap < plan.gaps.size(); ++gap) {
      if (!plan.gaps[gap].empty()) continue;
      const AnchorPair pair{anchors[gap], anchors[gap + 1]};
      const auto previous = previous_gaps.find(pair);
      if (previous == previous_gaps.end() ||
          previous->second != plan.gaps[gap]) {
        dirty.insert(pair);
      }
    }
  }
}

[[nodiscard]] std::vector<std::vector<int>> remove_customers(
    const Solution& solution, const std::vector<int>& removed) {
  const std::unordered_set<int> removed_set(removed.begin(), removed.end());
  std::vector<std::vector<int>> result;
  for (const Plan& plan : solution.plans) {
    std::vector<int> remaining;
    for (const int customer_id : plan.customer_ids) {
      if (!removed_set.contains(customer_id)) {
        remaining.push_back(customer_id);
      }
    }
    if (!remaining.empty()) {
      result.push_back(std::move(remaining));
    }
  }
  return result;
}

}  // namespace

GranularNeighborhood::GranularNeighborhood(const Instance& instance,
                                           std::size_t maximum_neighbors) {
  if (maximum_neighbors == 0) {
    throw std::invalid_argument("promising arc count must be positive");
  }
  const std::size_t count = instance.customer_ids().size();
  customer_rank_by_id_.reserve(count);
  for (std::size_t rank = 0; rank < count; ++rank) {
    customer_rank_by_id_.emplace(instance.customer_ids()[rank], rank);
  }
  promising_matrix_.assign(count * count, 0);
  if (count <= 1) {
    gamma_size_ = 0;
    return;
  }
  gamma_size_ = std::min(count - 1, maximum_neighbors);
  for (const int customer_id : instance.customer_ids()) {
    std::vector<int> nearest;
    nearest.reserve(count - 1);
    for (const int other_id : instance.customer_ids()) {
      if (other_id != customer_id) nearest.push_back(other_id);
    }
    std::sort(nearest.begin(), nearest.end(), [&](int left, int right) {
      return std::tuple{instance.distance_km(customer_id, left), left} <
             std::tuple{instance.distance_km(customer_id, right), right};
    });
    nearest.resize(gamma_size_);
    for (const int neighbor : nearest) {
      promising_edges_.insert(edge_key(customer_id, neighbor));
      promising_matrix_[customer_rank(customer_id) * count +
                         customer_rank(neighbor)] = 1;
    }
    directed_neighbors_.emplace(customer_id, std::move(nearest));
  }
  for (const int customer_id : instance.customer_ids()) {
    promising_neighbors_.try_emplace(customer_id);
  }
  for (const auto& [customer_id, neighbors] : directed_neighbors_) {
    for (const int neighbor : neighbors) {
      promising_neighbors_.at(customer_id).push_back(neighbor);
      promising_neighbors_.at(neighbor).push_back(customer_id);
    }
  }
  for (auto& [customer_id, neighbors] : promising_neighbors_) {
    std::sort(neighbors.begin(), neighbors.end());
    neighbors.erase(std::unique(neighbors.begin(), neighbors.end()),
                    neighbors.end());
    for (const int neighbor : neighbors) {
      if (customer_id < neighbor) {
        promising_pairs_.emplace_back(customer_id, neighbor);
      }
    }
  }
  std::sort(promising_pairs_.begin(), promising_pairs_.end());
}

const std::vector<int>& GranularNeighborhood::neighbors_of(
    int customer_id) const {
  return directed_neighbors_.at(customer_id);
}

const std::vector<int>& GranularNeighborhood::promising_neighbors_of(
    int customer_id) const {
  return promising_neighbors_.at(customer_id);
}

std::uint64_t GranularNeighborhood::edge_key(int left_id, int right_id) {
  return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(left_id))
          << 32U) |
         static_cast<std::uint64_t>(static_cast<std::uint32_t>(right_id));
}

std::size_t GranularNeighborhood::customer_rank(int customer_id) const {
  return customer_rank_by_id_.at(customer_id);
}

bool GranularNeighborhood::promising_by_rank(std::size_t left_rank,
                                             std::size_t right_rank) const {
  const std::size_t count = customer_rank_by_id_.size();
  if (left_rank >= count || right_rank >= count) {
    throw std::out_of_range("promising arc customer rank is out of range");
  }
  return promising_matrix_[left_rank * count + right_rank] != 0;
}

bool GranularNeighborhood::promising(int left_id, int right_id) const {
  const auto left = customer_rank_by_id_.find(left_id);
  const auto right = customer_rank_by_id_.find(right_id);
  return left != customer_rank_by_id_.end() &&
         right != customer_rank_by_id_.end() &&
         promising_by_rank(left->second, right->second);
}

bool GranularNeighborhood::allows_insertion(std::span<const int> route,
                                            std::size_t position,
                                            int customer_id) const {
  if (position > route.size()) {
    throw std::out_of_range("granular insertion position is out of range");
  }
  return (position > 0 && promising(route[position - 1], customer_id)) ||
         (position < route.size() && promising(customer_id, route[position]));
}

GranularNeighborhood::EdgeSet GranularNeighborhood::edge_set(
    const std::vector<std::vector<int>>& sequences) const {
  EdgeSet result;
  for (const auto& sequence : sequences) {
    for (std::size_t index = 1; index < sequence.size(); ++index) {
      result.insert(edge_key(sequence[index - 1], sequence[index]));
    }
  }
  return result;
}

bool GranularNeighborhood::allows(
    const std::vector<std::vector<int>>& current,
    const std::vector<std::vector<int>>& candidate) const {
  return allows(edge_set(current), candidate);
}

bool GranularNeighborhood::allows(
    const EdgeSet& current_edges,
    const std::vector<std::vector<int>>& candidate) const {
  if (promising_edges_.empty()) return true;
  for (const auto& sequence : candidate) {
    for (std::size_t index = 1; index < sequence.size(); ++index) {
      const std::uint64_t edge = edge_key(sequence[index - 1], sequence[index]);
      if (current_edges.contains(edge)) continue;
      if (promising_edges_.contains(edge)) return true;
    }
  }
  // Vidal's restriction requires every tested move to introduce at least one
  // promising customer arc; depot-only changes have no such arc.
  return false;
}

PathSampler::PathSampler(const Instance& instance, PlanFactory& factory,
                         std::mt19937_64& random)
    : instance_(instance), factory_(factory), random_(random) {
  build_curve_order();
  build_reachable_station_index();
  build_shortest_station_paths();
}

void PathSampler::build_curve_order() {
  std::vector<std::string> curve_types;
  curve_types.reserve(instance_.charging_curves().size());
  for (const auto& [station_type, curve] : instance_.charging_curves()) {
    (void)curve;
    curve_types.push_back(station_type);
  }
  std::sort(curve_types.begin(), curve_types.end());
  for (std::size_t left = 0; left < curve_types.size(); ++left) {
    for (std::size_t right = left + 1; right < curve_types.size(); ++right) {
      const ChargingCurve& left_curve =
          instance_.charging_curves().at(curve_types[left]);
      const ChargingCurve& right_curve =
          instance_.charging_curves().at(curve_types[right]);
      if (!charging_curve_no_slower(left_curve, right_curve) &&
          !charging_curve_no_slower(right_curve, left_curve)) {
        throw std::invalid_argument(
            "paper path dominance requires comparable charging curves");
      }
    }
  }
  std::sort(curve_types.begin(), curve_types.end(),
            [&](const std::string& left, const std::string& right) {
              const bool left_no_slower = charging_curve_no_slower(
                  instance_.charging_curves().at(left),
                  instance_.charging_curves().at(right));
              const bool right_no_slower = charging_curve_no_slower(
                  instance_.charging_curves().at(right),
                  instance_.charging_curves().at(left));
              if (left_no_slower != right_no_slower) return left_no_slower;
              return left < right;
            });
  for (std::size_t rank = 0; rank < curve_types.size(); ++rank) {
    curve_rank_by_type_.emplace(curve_types[rank], rank);
  }
  station_curve_ranks_.reserve(instance_.station_ids().size());
  for (const int station_id : instance_.station_ids()) {
    station_curve_ranks_.push_back(
        curve_rank_by_type_.at(instance_.node(station_id).station_type));
  }
}

void PathSampler::build_reachable_station_index() {
  const auto& stations = instance_.station_ids();
  const double capacity = instance_.vehicle().battery_capacity_wh;
  for (const Node& anchor : instance_.nodes()) {
    std::vector<ReachableStation>& reachable_from =
        stations_reachable_from_anchor_[anchor.id];
    std::vector<ReachableStation>& reaching =
        stations_reaching_anchor_[anchor.id];
    reachable_from.reserve(stations.size());
    reaching.reserve(stations.size());
    for (std::size_t station_index = 0; station_index < stations.size();
         ++station_index) {
      const double distance_from =
          instance_.distance_km(anchor.id, stations[station_index]);
      if (distance_from * instance_.vehicle().consumption_wh_per_km <=
          capacity + kEnergyToleranceWh) {
        reachable_from.push_back(ReachableStation{
            .station_index = station_index, .distance_km = distance_from});
      }
      const double distance_to =
          instance_.distance_km(stations[station_index], anchor.id);
      if (distance_to * instance_.vehicle().consumption_wh_per_km <=
          capacity + kEnergyToleranceWh) {
        reaching.push_back(ReachableStation{
            .station_index = station_index, .distance_km = distance_to});
      }
    }
  }
}

void PathSampler::build_shortest_station_paths() {
  const auto& stations = instance_.station_ids();
  const std::size_t count = stations.size();
  const std::size_t curve_count = curve_rank_by_type_.size();
  const double infinity = std::numeric_limits<double>::infinity();
  std::vector<std::vector<double>> edges(
      count, std::vector<double>(count, infinity));
  const double capacity = instance_.vehicle().battery_capacity_wh;
  for (std::size_t left = 0; left < count; ++left) {
    edges[left][left] = 0.0;
    for (std::size_t right = 0; right < count; ++right) {
      if (left == right) continue;
      const double distance =
          instance_.distance_km(stations[left], stations[right]);
      if (distance * instance_.vehicle().consumption_wh_per_km <=
          capacity + kEnergyToleranceWh) {
        edges[left][right] = distance;
      }
    }
  }
  shortest_station_paths_.assign(
      count, std::vector<std::vector<StationPath>>(count));
  struct QueueState {
    double distance{};
    std::size_t station{};
    std::size_t fastest_curve_rank{};
    std::size_t slowest_curve_rank{};
    std::vector<int> path;
  };
  const auto farther = [](const QueueState& left, const QueueState& right) {
    if (left.distance != right.distance) {
      return left.distance > right.distance;
    }
    if (left.path != right.path) return left.path > right.path;
    return std::tie(left.station, left.fastest_curve_rank,
                    left.slowest_curve_rank) >
           std::tie(right.station, right.fastest_curve_rank,
                    right.slowest_curve_rank);
  };
  for (std::size_t source = 0; source < count; ++source) {
    const auto state_index = [&](std::size_t station,
                                 std::size_t fastest_curve_rank,
                                 std::size_t slowest_curve_rank) {
      return (station * curve_count + fastest_curve_rank) * curve_count +
             slowest_curve_rank;
    };
    const std::size_t state_count = count * curve_count * curve_count;
    std::vector<double> distances(state_count, infinity);
    std::vector<std::vector<int>> paths(state_count);
    std::priority_queue<QueueState, std::vector<QueueState>,
                        decltype(farther)>
        queue(farther);
    const std::size_t source_rank = station_curve_ranks_[source];
    const std::size_t source_state =
        state_index(source, source_rank, source_rank);
    distances[source_state] = 0.0;
    paths[source_state] = {stations[source]};
    queue.push(QueueState{.distance = 0.0,
                          .station = source,
                          .fastest_curve_rank = source_rank,
                          .slowest_curve_rank = source_rank,
                          .path = paths[source_state]});
    while (!queue.empty()) {
      QueueState current = queue.top();
      queue.pop();
      const std::size_t current_state = state_index(
          current.station, current.fastest_curve_rank,
          current.slowest_curve_rank);
      if (current.distance != distances[current_state] ||
          current.path != paths[current_state]) {
        continue;
      }
      for (std::size_t next = 0; next < count; ++next) {
        if (next == current.station ||
            !std::isfinite(edges[current.station][next]) ||
            std::find(current.path.begin(), current.path.end(),
                      stations[next]) != current.path.end()) {
          continue;
        }
        const double candidate_distance =
            current.distance + edges[current.station][next];
        std::vector<int> candidate_path = current.path;
        candidate_path.push_back(stations[next]);
        const std::size_t candidate_fastest = std::min(
            current.fastest_curve_rank, station_curve_ranks_[next]);
        const std::size_t candidate_slowest = std::max(
            current.slowest_curve_rank, station_curve_ranks_[next]);
        const std::size_t candidate_state =
            state_index(next, candidate_fastest, candidate_slowest);
        if (candidate_distance < distances[candidate_state] ||
            (candidate_distance == distances[candidate_state] &&
             (paths[candidate_state].empty() ||
              candidate_path < paths[candidate_state]))) {
          distances[candidate_state] = candidate_distance;
          paths[candidate_state] = candidate_path;
          queue.push(QueueState{.distance = candidate_distance,
                                .station = next,
                                .fastest_curve_rank = candidate_fastest,
                                .slowest_curve_rank = candidate_slowest,
                                .path = std::move(candidate_path)});
        }
      }
    }
    for (std::size_t target = 0; target < count; ++target) {
      for (std::size_t fastest = 0; fastest < curve_count; ++fastest) {
        for (std::size_t slowest = fastest; slowest < curve_count;
             ++slowest) {
          const std::size_t state = state_index(target, fastest, slowest);
          if (paths[state].empty()) continue;
          shortest_station_paths_[source][target].push_back(StationPath{
              .distance = distances[state],
              .stations = std::move(paths[state]),
              .fastest_curve_rank = fastest,
              .slowest_curve_rank = slowest});
        }
      }
      std::sort(shortest_station_paths_[source][target].begin(),
                shortest_station_paths_[source][target].end(),
                [](const StationPath& left, const StationPath& right) {
                  return std::tie(left.fastest_curve_rank,
                                  left.slowest_curve_rank, left.distance,
                                  left.stations) <
                         std::tie(right.fastest_curve_rank,
                                  right.slowest_curve_rank, right.distance,
                                  right.stations);
                });
    }
  }
}

std::vector<std::vector<int>> PathSampler::alternatives(
    int origin_id, int target_id, const std::vector<int>& current_path) {
  const auto options = path_options(origin_id, target_id);
  std::vector<std::vector<int>> result;
  result.reserve(options->size());
  for (const auto& option : *options) {
    if (option.stations != current_path) result.push_back(option.stations);
  }
  std::sort(result.begin(), result.end());
  return result;
}

std::shared_ptr<const PathSampler::PathOptions> PathSampler::path_options(
    int origin_id, int target_id, RepairTimingStatistics* timing) {
  const AnchorPair key{origin_id, target_id};
  if (auto cached = option_cache_.get(key); cached.has_value()) {
    if (timing != nullptr) ++timing->path_option_cache_hits;
    return *cached;
  }
  if (timing != nullptr) ++timing->path_option_cache_misses;
  if (option_keys_built_.insert(key).second) {
    ++option_first_builds_;
  } else {
    ++option_rebuilds_;
  }
  const auto build_started = std::chrono::steady_clock::now();
  struct CandidatePath {
    const StationPath* station_path{};
    double first_leg{};
    double last_leg{};
    double total_distance{};
  };
  std::vector<CandidatePath> candidates;
  const auto& reachable_first =
      stations_reachable_from_anchor_.at(origin_id);
  const auto& reachable_last = stations_reaching_anchor_.at(target_id);
  candidates.reserve(reachable_first.size() * reachable_last.size());
  for (const ReachableStation& first : reachable_first) {
    for (const ReachableStation& last : reachable_last) {
      for (const StationPath& station_path :
           shortest_station_paths_[first.station_index][last.station_index]) {
        candidates.push_back(CandidatePath{
            .station_path = &station_path,
            .first_leg = first.distance_km,
            .last_leg = last.distance_km,
            .total_distance = first.distance_km + station_path.distance +
                              last.distance_km});
      }
    }
  }
  option_candidates_generated_ += candidates.size();
  if (timing != nullptr) {
    timing->path_option_candidates_generated += candidates.size();
  }

  std::sort(candidates.begin(), candidates.end(), [](const CandidatePath& left,
                                                      const CandidatePath& right) {
    return std::tie(left.first_leg, left.last_leg, left.total_distance,
                    left.station_path->fastest_curve_rank,
                    left.station_path->slowest_curve_rank,
                    left.station_path->stations) <
           std::tie(right.first_leg, right.last_leg, right.total_distance,
                    right.station_path->fastest_curve_rank,
                    right.station_path->slowest_curve_rank,
                    right.station_path->stations);
  });
  const auto dominates = [](const CandidatePath& left,
                            const CandidatePath& right) {
    // TRC Section 4.3 requires even the slowest curve on the dominating
    // path to be no slower than the fastest curve on the dominated path.
    return left.station_path->slowest_curve_rank <=
               right.station_path->fastest_curve_rank &&
           left.first_leg <= right.first_leg + kCostTolerance &&
           left.last_leg <= right.last_leg + kCostTolerance &&
           left.total_distance <= right.total_distance + kCostTolerance;
  };
  struct RankedCandidate {
    CandidatePath candidate;
    bool active{true};
  };
  const std::size_t curve_count = curve_rank_by_type_.size();
  std::vector<RankedCandidate> ranked;
  ranked.reserve(candidates.size());
  std::vector<std::vector<std::size_t>> by_fastest(curve_count);
  std::vector<std::vector<std::size_t>> by_slowest(curve_count);
  for (CandidatePath& candidate : candidates) {
    bool is_dominated = false;
    for (std::size_t slowest = 0;
         slowest <= candidate.station_path->fastest_curve_rank;
         ++slowest) {
      for (const std::size_t index : by_slowest[slowest]) {
        if (ranked[index].active &&
            dominates(ranked[index].candidate, candidate)) {
          is_dominated = true;
          break;
        }
      }
      if (is_dominated) break;
    }
    if (is_dominated) {
      continue;
    }
    for (std::size_t fastest = candidate.station_path->slowest_curve_rank;
         fastest < curve_count; ++fastest) {
      for (const std::size_t index : by_fastest[fastest]) {
        if (ranked[index].active &&
            dominates(candidate, ranked[index].candidate)) {
          ranked[index].active = false;
        }
      }
    }
    const std::size_t index = ranked.size();
    const std::size_t fastest = candidate.station_path->fastest_curve_rank;
    const std::size_t slowest = candidate.station_path->slowest_curve_rank;
    ranked.push_back(
        RankedCandidate{.candidate = std::move(candidate), .active = true});
    by_fastest[fastest].push_back(index);
    by_slowest[slowest].push_back(index);
  }
  std::vector<CandidatePath> nondominated;
  nondominated.reserve(ranked.size());
  for (RankedCandidate& candidate : ranked) {
    if (candidate.active) {
      nondominated.push_back(std::move(candidate.candidate));
    }
  }
  option_nondominated_generated_ += nondominated.size();
  if (timing != nullptr) {
    timing->path_option_nondominated += nondominated.size();
  }

  auto result = std::make_shared<PathOptions>();
  result->reserve(nondominated.size() + 1);
  NonlinearSequenceEvaluator sequence_evaluator(instance_);
  result->push_back(PathOption{
      .stations = {},
      .total_distance_km = instance_.distance_km(origin_id, target_id),
      .label = sequence_evaluator.edge_label(origin_id, {}, target_id)});
  for (const CandidatePath& candidate : nondominated) {
    std::vector<int> path_stations = candidate.station_path->stations;
    ChargingSequenceLabel label = sequence_evaluator.edge_label(
        origin_id, path_stations, target_id);
    result->push_back(PathOption{
        .stations = std::move(path_stations),
        .total_distance_km = candidate.total_distance,
        .label = std::move(label)});
  }
  std::sort(result->begin(), result->end(),
            [](const PathOption& left, const PathOption& right) {
              return std::tuple{left.total_distance_km,
                                left.stations.size(), left.stations} <
                     std::tuple{right.total_distance_km,
                                right.stations.size(), right.stations};
            });
  const std::size_t evictions_before = option_cache_.evictions();
  option_cache_.put(key, result);
  if (timing != nullptr) {
    timing->path_option_cache_evictions +=
        option_cache_.evictions() - evictions_before;
    timing->path_option_cache_miss_build_s += elapsed_seconds(build_started);
  }
  return result;
}

std::vector<std::vector<int>> PathSampler::paths_between(int origin_id,
                                                         int target_id) {
  const auto options = path_options(origin_id, target_id);
  std::vector<std::vector<int>> paths;
  paths.reserve(options->size());
  for (const PathOption& option : *options) {
    paths.push_back(option.stations);
  }
  std::sort(paths.begin(), paths.end());
  return paths;
}

VirtualInsertionCompletion PathSampler::virtual_complete_insertion(
    const Plan* source_plan, std::size_t position, int customer_id,
    double penalty_lambda, RepairTimingStatistics* timing) {
  const std::size_t source_size =
      source_plan == nullptr ? 0 : source_plan->customer_ids.size();
  if (position > source_size) {
    throw std::out_of_range("virtual insertion has invalid position");
  }
  VirtualCompletionCacheKey completion_key{
      .source_physical_visits =
          source_plan == nullptr ? std::vector<int>{}
                                 : source_plan->route.visits,
      .customer_id = customer_id,
      .position = position,
      .penalty_lambda_bits = std::bit_cast<std::uint64_t>(penalty_lambda),
      .semantics_version = 1,
      .opens_route = source_plan == nullptr};
  if (auto cached = virtual_completion_cache_.get(completion_key);
      cached.has_value()) {
    if (timing != nullptr) ++timing->virtual_completion_cross_cache_hits;
    return **cached;
  }
  if (timing != nullptr) ++timing->virtual_completion_cross_cache_misses;

  std::vector<int> customers;
  if (source_plan != nullptr) customers = source_plan->customer_ids;
  customers.insert(customers.begin() + static_cast<std::ptrdiff_t>(position),
                   customer_id);
  std::vector<int> anchors;
  anchors.reserve(customers.size() + 2);
  anchors.push_back(instance_.depot().id);
  anchors.insert(anchors.end(), customers.begin(), customers.end());
  anchors.push_back(instance_.depot().id);

  const std::array<std::size_t, 2> inserted_gaps{position, position + 1};
  std::vector<std::vector<int>> gaps(customers.size() + 1);
  if (source_plan != nullptr) {
    if (source_plan->gaps.size() != source_size + 1) {
      throw std::logic_error("source plan has invalid gap count");
    }
    for (std::size_t gap = 0; gap < gaps.size(); ++gap) {
      if (gap == inserted_gaps[0] || gap == inserted_gaps[1]) continue;
      gaps[gap] = source_plan->gaps[gap < position ? gap : gap - 1];
    }
  }

  NonlinearSequenceEvaluator evaluator(instance_);
  std::vector<ChargingSequenceLabel> edges;
  edges.reserve(gaps.size());
  for (std::size_t gap = 0; gap < gaps.size(); ++gap) {
    edges.push_back(
        evaluator.edge_label(anchors[gap], gaps[gap], anchors[gap + 1]));
  }
  const auto evaluate_edges = [&]() {
    ChargingSequenceLabel route;
    for (const ChargingSequenceLabel& edge : edges) {
      evaluator.append(route, edge);
    }
    return evaluator.evaluate(route);
  };

  const RouteCostSummary direct_summary = evaluate_edges();
  RouteCostSummary incumbent_summary = direct_summary;
  double incumbent_cost =
      generalized_cost(incumbent_summary, instance_, penalty_lambda);
  const std::array<std::shared_ptr<const PathOptions>, 2> options{
      path_options(anchors[inserted_gaps[0]],
                   anchors[inserted_gaps[0] + 1], timing),
      path_options(anchors[inserted_gaps[1]],
                   anchors[inserted_gaps[1] + 1], timing)};
  std::array<std::vector<int>, 2> selected_paths;
  std::uint64_t options_scored = 0;
  std::uint64_t options_pruned_by_lower_bound = 0;

  // The insertion changes exactly two adjacent gaps.  Preserve the historical
  // concatenation tree at their boundaries, but build the unchanged outer
  // prefix and suffix only once.  Each fixed-point pass then updates just the
  // prefix/suffix context affected by the two current edge choices.
  const std::size_t left_gap = inserted_gaps[0];
  const std::size_t right_gap = inserted_gaps[1];
  ChargingSequenceLabel outer_prefix;
  for (std::size_t gap = 0; gap < left_gap; ++gap) {
    outer_prefix = evaluator.concatenate(outer_prefix, edges[gap]);
  }
  ChargingSequenceLabel outer_suffix;
  for (std::size_t gap = edges.size(); gap > right_gap + 1; --gap) {
    outer_suffix =
        evaluator.concatenate(edges[gap - 1], outer_suffix);
  }

  const auto option_scan_started = std::chrono::steady_clock::now();
  while (true) {
    const ChargingSequenceLabel right_prefix =
        evaluator.concatenate(outer_prefix, edges[left_gap]);
    const ChargingSequenceLabel left_suffix =
        evaluator.concatenate(edges[right_gap], outer_suffix);
    const std::array<const ChargingSequenceLabel*, 2> scoring_prefixes{
        &outer_prefix, &right_prefix};
    const std::array<const ChargingSequenceLabel*, 2> scoring_suffixes{
        &left_suffix, &outer_suffix};

    double minimum_cost = std::numeric_limits<double>::infinity();
    virtual_tolerance_band_.clear();
    for (std::size_t side = 0; side < inserted_gaps.size(); ++side) {
      const PathOptions& side_options = *options[side];
      VirtualOptionScanWorkspace& workspace = virtual_option_workspaces_[side];
      workspace.lower_bounds.resize(side_options.size());
      workspace.suffix_minimum_lower_bound.resize(side_options.size() + 1);
      workspace.suffix_candidate_count.resize(side_options.size() + 1);
      workspace.suffix_minimum_lower_bound[side_options.size()] =
          std::numeric_limits<double>::infinity();
      workspace.suffix_candidate_count[side_options.size()] = 0;
      for (std::size_t option_index = side_options.size(); option_index-- > 0;) {
        const PathOption& option = side_options[option_index];
        workspace.suffix_minimum_lower_bound[option_index] =
            workspace.suffix_minimum_lower_bound[option_index + 1];
        workspace.suffix_candidate_count[option_index] =
            workspace.suffix_candidate_count[option_index + 1];
        if (option.stations == selected_paths[side]) continue;
        const double travel_lower_bound =
            scoring_prefixes[side]->travel_time_h +
            option.label.travel_time_h +
            scoring_suffixes[side]->travel_time_h;
        const double service_lower_bound =
            scoring_prefixes[side]->service_time_h +
            option.label.service_time_h +
            scoring_suffixes[side]->service_time_h;
        const double generalized_lower_bound =
            travel_lower_bound +
            penalty_lambda *
                std::max(0.0, travel_lower_bound + service_lower_bound -
                                  instance_.vehicle().max_route_duration_h);
        workspace.lower_bounds[option_index] = generalized_lower_bound;
        ++workspace.suffix_candidate_count[option_index];
        // A non-finite value cannot certify that any remaining option is
        // prunable.  The ordinary per-option test below retains its historical
        // IEEE-754 behavior for such a value.
        const double certified_lower_bound =
            std::isfinite(generalized_lower_bound)
                ? generalized_lower_bound
                : -std::numeric_limits<double>::infinity();
        workspace.suffix_minimum_lower_bound[option_index] = std::min(
            certified_lower_bound,
            workspace.suffix_minimum_lower_bound[option_index + 1]);
      }
      for (std::size_t option_index = 0;
           option_index < side_options.size(); ++option_index) {
        const PathOption& option = side_options[option_index];
        if (option.stations == selected_paths[side]) continue;
        if (workspace.suffix_minimum_lower_bound[option_index] >=
            incumbent_cost - kCostTolerance) {
          options_scored += workspace.suffix_candidate_count[option_index];
          options_pruned_by_lower_bound +=
              workspace.suffix_candidate_count[option_index];
          break;
        }
        ++options_scored;
        const double generalized_lower_bound =
            workspace.lower_bounds[option_index];
        if (generalized_lower_bound >=
            incumbent_cost - kCostTolerance) {
          ++options_pruned_by_lower_bound;
          continue;
        }
        ChargingSequenceLabel candidate = *scoring_prefixes[side];
        evaluator.append(candidate, option.label);
        evaluator.append(candidate, *scoring_suffixes[side]);
        const double cost = evaluator.generalized_cost(candidate,
                                                       penalty_lambda);
        if (cost >= incumbent_cost - kCostTolerance) continue;
        if (cost < minimum_cost) {
          minimum_cost = cost;
          std::erase_if(virtual_tolerance_band_,
                        [&](const VirtualOptionCandidate& item) {
            return item.cost > minimum_cost + kCostTolerance;
          });
        }
        if (cost <= minimum_cost + kCostTolerance) {
          virtual_tolerance_band_.push_back(VirtualOptionCandidate{
              .cost = cost, .side = side, .option_index = option_index});
        }
      }
    }
    if (virtual_tolerance_band_.empty()) break;

    const VirtualOptionCandidate* best = nullptr;
    std::vector<int> best_visits;
    for (const VirtualOptionCandidate& candidate : virtual_tolerance_band_) {
      const PathOption& option =
          (*options[candidate.side])[candidate.option_index];
      auto candidate_gaps = gaps;
      candidate_gaps[inserted_gaps[candidate.side]] =
          option.stations;
      std::vector<int> visits =
          expand_route(customers, candidate_gaps).visits;
      if (best == nullptr || visits < best_visits) {
        best = &candidate;
        best_visits = std::move(visits);
      }
    }
    const PathOption& best_option =
        (*options[best->side])[best->option_index];
    selected_paths[best->side] = best_option.stations;
    gaps[inserted_gaps[best->side]] = best_option.stations;
    edges[inserted_gaps[best->side]] = best_option.label;
    incumbent_summary = evaluate_edges();
    incumbent_cost =
        generalized_cost(incumbent_summary, instance_, penalty_lambda);
  }
  if (timing != nullptr) {
    timing->nonlinear_option_scanning_s +=
        elapsed_seconds(option_scan_started);
  }

  VirtualInsertionCompletion result{
      .summary = incumbent_summary,
      .generalized_cost_h = incumbent_cost,
      .gaps = {AnchorPair{anchors[inserted_gaps[0]],
                          anchors[inserted_gaps[0] + 1]},
               AnchorPair{anchors[inserted_gaps[1]],
                          anchors[inserted_gaps[1] + 1]}},
      .paths = std::move(selected_paths),
      .options_scored = options_scored,
      .options_pruned_by_lower_bound = options_pruned_by_lower_bound,
      .changed = !gaps[inserted_gaps[0]].empty() ||
                 !gaps[inserted_gaps[1]].empty(),
      .infeasible_to_feasible =
          !direct_summary.feasible && incumbent_summary.feasible};
  const std::size_t evictions_before = virtual_completion_cache_.evictions();
  virtual_completion_cache_.put(
      std::move(completion_key),
      std::make_shared<const VirtualInsertionCompletion>(result));
  if (timing != nullptr) {
    timing->virtual_completion_cross_cache_evictions +=
        virtual_completion_cache_.evictions() - evictions_before;
  }
  return result;
}

HintedPathResult PathSampler::apply_insertion_hints(
    Solution solution, std::span<const InsertionPathHint> hints,
    double penalty_lambda,
    std::optional<std::chrono::steady_clock::time_point> deadline) {
  HintedPathResult result{.solution = std::move(solution)};
  for (const InsertionPathHint& hint : hints) {
    if (deadline_reached(deadline)) break;
    bool found = false;
    for (Plan& plan : result.solution.plans) {
      if (plan.exact_charging) continue;
      std::vector<int> anchors;
      anchors.reserve(plan.customer_ids.size() + 2);
      anchors.push_back(instance_.depot().id);
      anchors.insert(anchors.end(), plan.customer_ids.begin(),
                     plan.customer_ids.end());
      anchors.push_back(instance_.depot().id);
      for (std::size_t gap = 0; gap < plan.gaps.size(); ++gap) {
        if (AnchorPair{anchors[gap], anchors[gap + 1]} != hint.gap) continue;
        found = true;
        if (plan.gaps[gap] == hint.stations) break;
        ++result.attempts;
        auto improved = improving_atomic_replacement(
            plan, gap, hint.stations, penalty_lambda, deadline);
        if (improved.has_value()) {
          plan = std::move(*improved);
          ++result.accepted;
        }
        break;
      }
      if (found) break;
    }
  }
  return result;
}

std::optional<Solution> PathSampler::random_replacement(
    const Solution& solution) {
  std::size_t gap_count = 0;
  for (const Plan& plan : solution.plans) gap_count += plan.gaps.size();
  if (gap_count == 0) return std::nullopt;
  std::size_t selected_ordinal =
      std::uniform_int_distribution<std::size_t>(0, gap_count - 1)(random_);
  std::size_t selected_route = 0;
  while (selected_ordinal >= solution.plans[selected_route].gaps.size()) {
    selected_ordinal -= solution.plans[selected_route].gaps.size();
    ++selected_route;
  }
  const std::size_t selected_gap = selected_ordinal;
  Solution candidate = solution;
  Plan& plan = candidate.plans[selected_route];
  std::vector<int> anchors;
  anchors.reserve(plan.customer_ids.size() + 2);
  anchors.push_back(instance_.depot().id);
  anchors.insert(anchors.end(), plan.customer_ids.begin(),
                 plan.customer_ids.end());
  anchors.push_back(instance_.depot().id);
  auto options = alternatives(anchors[selected_gap], anchors[selected_gap + 1],
                              plan.gaps[selected_gap]);
  if (options.empty()) return std::nullopt;
  auto gaps = plan.gaps;
  gaps[selected_gap] = options[random_index(options, random_)];
  plan = factory_.make_plan(plan.customer_ids, std::move(gaps));
  return candidate;
}

std::optional<Plan> PathSampler::improving_atomic_replacement(
    const Plan& plan, std::size_t gap, const std::vector<int>& stations,
    double penalty_lambda,
    std::optional<std::chrono::steady_clock::time_point> deadline) {
  if (plan.exact_charging || gap >= plan.gaps.size() ||
      stations == plan.gaps[gap] || deadline_reached(deadline)) {
    return std::nullopt;
  }
  std::vector<int> anchors;
  anchors.reserve(plan.customer_ids.size() + 2);
  anchors.push_back(instance_.depot().id);
  anchors.insert(anchors.end(), plan.customer_ids.begin(),
                 plan.customer_ids.end());
  anchors.push_back(instance_.depot().id);

  NonlinearSequenceEvaluator sequence_evaluator(instance_);
  std::vector<ChargingSequenceLabel> edges;
  edges.reserve(plan.gaps.size());
  for (std::size_t index = 0; index < plan.gaps.size(); ++index) {
    if (deadline_reached(deadline)) return std::nullopt;
    edges.push_back(sequence_evaluator.edge_label(
        anchors[index], plan.gaps[index], anchors[index + 1]));
  }
  ChargingSequenceLabel prefix;
  for (std::size_t index = 0; index < gap; ++index) {
    sequence_evaluator.append(prefix, edges[index]);
  }
  ChargingSequenceLabel suffix;
  for (std::size_t index = plan.gaps.size(); index-- > gap + 1;) {
    suffix = sequence_evaluator.concatenate(edges[index], suffix);
  }
  const ChargingSequenceLabel replacement = sequence_evaluator.edge_label(
      anchors[gap], stations, anchors[gap + 1]);
  ChargingSequenceLabel candidate = prefix;
  sequence_evaluator.append(candidate, replacement);
  sequence_evaluator.append(candidate, suffix);
  const double incumbent_fast_cost =
      sequence_evaluator.generalized_cost(
          sequence_evaluator.route_label(plan.route), penalty_lambda);
  const double candidate_fast_cost =
      sequence_evaluator.generalized_cost(candidate, penalty_lambda);
  if (candidate_fast_cost >= incumbent_fast_cost - kCostTolerance ||
      deadline_reached(deadline)) {
    return std::nullopt;
  }

  auto gaps = plan.gaps;
  gaps[gap] = stations;
  Plan materialized = factory_.make_plan(plan.customer_ids, std::move(gaps));
  if (deadline_reached(deadline) ||
      generalized_cost(*materialized.evaluation, instance_, penalty_lambda) >=
          generalized_cost(*plan.evaluation, instance_, penalty_lambda) -
              kCostTolerance) {
    return std::nullopt;
  }
  return materialized;
}

LimitedRandomPathResult PathSampler::limited_random_replacement(
    Solution solution, double penalty_lambda,
    std::span<const AnchorPair> eligible_anchor_pairs,
    std::optional<std::chrono::steady_clock::time_point> deadline) {
  LimitedRandomPathResult result{.solution = std::move(solution)};
  if (eligible_anchor_pairs.empty() || deadline_reached(deadline)) {
    return result;
  }
  const std::unordered_set<AnchorPair, AnchorPairHash> eligible(
      eligible_anchor_pairs.begin(), eligible_anchor_pairs.end());
  struct EligibleGap {
    std::size_t route{};
    std::size_t gap{};
    std::shared_ptr<const PathOptions> options;
  };
  std::vector<EligibleGap> gaps;
  for (std::size_t route = 0; route < result.solution.plans.size(); ++route) {
    const Plan& plan = result.solution.plans[route];
    if (plan.exact_charging) continue;
    std::vector<int> anchors;
    anchors.reserve(plan.customer_ids.size() + 2);
    anchors.push_back(instance_.depot().id);
    anchors.insert(anchors.end(), plan.customer_ids.begin(),
                   plan.customer_ids.end());
    anchors.push_back(instance_.depot().id);
    for (std::size_t gap = 0; gap < plan.gaps.size(); ++gap) {
      if (!eligible.contains(AnchorPair{anchors[gap], anchors[gap + 1]}))
        continue;
      const auto options = path_options(anchors[gap], anchors[gap + 1]);
      if (std::none_of(options->begin(), options->end(),
                       [&](const PathOption& option) {
                         return option.stations != plan.gaps[gap];
                       })) {
        continue;
      }
      gaps.push_back(EligibleGap{.route = route,
                                 .gap = gap,
                                 .options = options});
    }
  }
  result.eligible_gaps = gaps.size();
  if (gaps.empty()) return result;
  result.failure_limit = std::min<std::size_t>(
      gaps.size(), std::clamp<std::size_t>(
                       static_cast<std::size_t>(std::ceil(
                           std::sqrt(static_cast<double>(gaps.size())))),
                       4, 16));

  struct AttemptedAlternative {
    std::size_t gap{};
    std::vector<int> stations;
  };
  std::vector<std::vector<AttemptedAlternative>> attempted_by_route(
      result.solution.plans.size());
  std::size_t consecutive_failures = 0;
  while (consecutive_failures < result.failure_limit &&
         !deadline_reached(deadline)) {
    std::vector<std::size_t> selectable_gaps;
    selectable_gaps.reserve(gaps.size());
    for (std::size_t gap_index = 0; gap_index < gaps.size(); ++gap_index) {
      const EligibleGap& candidate_gap = gaps[gap_index];
      const Plan& plan = result.solution.plans[candidate_gap.route];
      const auto& attempted = attempted_by_route[candidate_gap.route];
      if (std::any_of(candidate_gap.options->begin(),
                      candidate_gap.options->end(),
                      [&](const PathOption& option) {
                        if (option.stations == plan.gaps[candidate_gap.gap])
                          return false;
                        return std::none_of(
                            attempted.begin(), attempted.end(),
                            [&](const AttemptedAlternative& previous) {
                              return previous.gap == candidate_gap.gap &&
                                     previous.stations == option.stations;
                            });
                      })) {
        selectable_gaps.push_back(gap_index);
      }
    }
    if (selectable_gaps.empty()) break;

    const EligibleGap& selected_gap =
        gaps[selectable_gaps[random_index(selectable_gaps, random_)]];
    const Plan& incumbent = result.solution.plans[selected_gap.route];
    const auto& attempted = attempted_by_route[selected_gap.route];
    std::vector<const PathOption*> untried;
    for (const PathOption& option : *selected_gap.options) {
      if (option.stations == incumbent.gaps[selected_gap.gap]) continue;
      if (std::any_of(attempted.begin(), attempted.end(),
                      [&](const AttemptedAlternative& previous) {
                        return previous.gap == selected_gap.gap &&
                               previous.stations == option.stations;
                      })) {
        continue;
      }
      untried.push_back(&option);
    }
    if (untried.empty()) continue;
    const PathOption& selected_option =
        *untried[random_index(untried, random_)];
    attempted_by_route[selected_gap.route].push_back(AttemptedAlternative{
        .gap = selected_gap.gap, .stations = selected_option.stations});
    ++result.attempts;
    auto improved = improving_atomic_replacement(
        incumbent, selected_gap.gap, selected_option.stations,
        penalty_lambda, deadline);
    if (improved.has_value()) {
      result.solution.plans[selected_gap.route] = std::move(*improved);
      attempted_by_route[selected_gap.route].clear();
      ++result.accepted;
      consecutive_failures = 0;
    } else {
      ++result.failed_attempts;
      ++consecutive_failures;
    }
  }
  return result;
}

QuasiExhaustivePathResult PathSampler::quasi_exhaustive_replacement(
    Solution solution, double penalty_lambda,
    std::span<const AnchorPair> eligible_anchor_pairs,
    std::optional<std::chrono::steady_clock::time_point> deadline) {
  QuasiExhaustivePathResult result{.solution = std::move(solution)};
  if (eligible_anchor_pairs.empty() || deadline_reached(deadline)) {
    return result;
  }
  const std::unordered_set<AnchorPair, AnchorPairHash> eligible(
      eligible_anchor_pairs.begin(), eligible_anchor_pairs.end());
  std::vector<std::pair<std::size_t, std::size_t>> gaps;
  for (std::size_t route = 0; route < result.solution.plans.size(); ++route) {
    const Plan& plan = result.solution.plans[route];
    if (plan.exact_charging) continue;
    std::vector<int> anchors;
    anchors.reserve(plan.customer_ids.size() + 2);
    anchors.push_back(instance_.depot().id);
    anchors.insert(anchors.end(), plan.customer_ids.begin(),
                   plan.customer_ids.end());
    anchors.push_back(instance_.depot().id);
    for (std::size_t gap = 0; gap < plan.gaps.size(); ++gap) {
      if (!eligible.contains(AnchorPair{anchors[gap], anchors[gap + 1]}))
        continue;
      const auto options = path_options(anchors[gap], anchors[gap + 1]);
      if (std::any_of(options->begin(), options->end(),
                      [&](const PathOption& option) {
                        return option.stations != plan.gaps[gap];
                      })) {
        gaps.emplace_back(route, gap);
      }
    }
  }
  std::shuffle(gaps.begin(), gaps.end(), random_);
  for (const auto [route, gap] : gaps) {
    if (deadline_reached(deadline)) break;
    ++result.gaps_scanned;
    const std::array<std::size_t, 1> eligible_gap{gap};
    auto improved = best_replacement_in_gaps(
        result.solution.plans[route], penalty_lambda, eligible_gap,
        std::span<const ExcludedPathOption>{}, deadline, nullptr);
    if (!improved.has_value()) continue;
    const double incumbent_cost = generalized_cost(
        *result.solution.plans[route].evaluation, instance_, penalty_lambda);
    const double improved_cost = generalized_cost(
        *improved->evaluation, instance_, penalty_lambda);
    if (improved_cost >= incumbent_cost - kCostTolerance) continue;
    result.solution.plans[route] = std::move(*improved);
    ++result.accepted;
  }
  return result;
}

std::optional<Plan> PathSampler::best_replacement(const Plan& plan,
                                                  double penalty_lambda) {
  std::vector<std::size_t> eligible(plan.gaps.size());
  std::iota(eligible.begin(), eligible.end(), 0);
  return best_replacement_in_gaps(
      plan, penalty_lambda, eligible,
      std::span<const ExcludedPathOption>{}, std::nullopt, nullptr);
}

std::optional<Plan> PathSampler::best_replacement_in_gaps(
    const Plan& plan, double penalty_lambda,
    std::span<const std::size_t> eligible_gaps,
    std::span<const ExcludedPathOption> excluded,
    std::optional<std::chrono::steady_clock::time_point> deadline,
    PathCompletionStatistics* statistics) {
  if (plan.exact_charging || deadline_reached(deadline)) return std::nullopt;
  if (statistics != nullptr) ++statistics->search_passes;
  const auto setup_started = std::chrono::steady_clock::now();
  std::vector<int> anchors;
  anchors.reserve(plan.customer_ids.size() + 2);
  anchors.push_back(instance_.depot().id);
  anchors.insert(anchors.end(), plan.customer_ids.begin(),
                 plan.customer_ids.end());
  anchors.push_back(instance_.depot().id);
  NonlinearSequenceEvaluator sequence_evaluator(instance_);
  std::vector<ChargingSequenceLabel> edges;
  edges.reserve(plan.gaps.size());
  for (std::size_t gap_index = 0; gap_index < plan.gaps.size(); ++gap_index) {
    if (deadline_reached(deadline)) return std::nullopt;
    edges.push_back(sequence_evaluator.edge_label(
        anchors[gap_index], plan.gaps[gap_index],
        anchors[gap_index + 1]));
  }
  for (const std::size_t gap_index : eligible_gaps) {
    if (gap_index >= plan.gaps.size()) {
      throw std::out_of_range("eligible path-repair gap is invalid");
    }
  }
  std::vector<ChargingSequenceLabel> prefixes(plan.gaps.size() + 1);
  for (std::size_t gap_index = 0; gap_index < plan.gaps.size(); ++gap_index) {
    if (deadline_reached(deadline)) return std::nullopt;
    prefixes[gap_index + 1] = sequence_evaluator.concatenate(
        prefixes[gap_index], edges[gap_index]);
  }
  std::vector<ChargingSequenceLabel> suffixes(plan.gaps.size() + 1);
  for (std::size_t gap_index = plan.gaps.size(); gap_index-- > 0;) {
    if (deadline_reached(deadline)) return std::nullopt;
    suffixes[gap_index] = sequence_evaluator.concatenate(
        edges[gap_index], suffixes[gap_index + 1]);
  }
  const double incumbent_cost =
      generalized_cost(*plan.evaluation, instance_, penalty_lambda);
  struct ReplacementCandidate {
    double cost{};
    std::size_t gap{};
    std::vector<int> stations;
  };
  double minimum_cost = std::numeric_limits<double>::infinity();
  std::vector<ReplacementCandidate> tolerance_band;
  if (statistics != nullptr) {
    statistics->setup_s += elapsed_seconds(setup_started);
  }
  for (const std::size_t gap_index : eligible_gaps) {
    if (deadline_reached(deadline)) break;
    const auto catalog_started = std::chrono::steady_clock::now();
    const std::size_t cache_hits_before = option_cache_.hits();
    const std::size_t cache_misses_before = option_cache_.misses();
    const auto options =
        path_options(anchors[gap_index], anchors[gap_index + 1]);
    if (statistics != nullptr) {
      statistics->option_catalog_s += elapsed_seconds(catalog_started);
      statistics->option_cache_hits +=
          option_cache_.hits() - cache_hits_before;
      statistics->option_cache_misses +=
          option_cache_.misses() - cache_misses_before;
    }
    if (deadline_reached(deadline)) break;
    const auto scoring_started = std::chrono::steady_clock::now();
    for (const PathOption& option : *options) {
      if (deadline_reached(deadline)) break;
      if (option.stations == plan.gaps[gap_index]) continue;
      if (std::any_of(excluded.begin(), excluded.end(),
                      [&](const ExcludedPathOption& item) {
                        return item.gap == gap_index &&
                               item.stations == option.stations;
      })) {
        continue;
      }
      if (statistics != nullptr) ++statistics->options_scored;
      ChargingSequenceLabel candidate = prefixes[gap_index];
      sequence_evaluator.append(candidate, option.label);
      sequence_evaluator.append(candidate, suffixes[gap_index + 1]);
      const double candidate_cost =
          sequence_evaluator.generalized_cost(candidate, penalty_lambda);
      if (candidate_cost >= incumbent_cost - kCostTolerance) continue;
      if (candidate_cost < minimum_cost) {
        minimum_cost = candidate_cost;
        std::erase_if(tolerance_band, [&](const ReplacementCandidate& item) {
          return item.cost > minimum_cost + kCostTolerance;
        });
      }
      if (candidate_cost <= minimum_cost + kCostTolerance) {
        tolerance_band.push_back(ReplacementCandidate{
            .cost = candidate_cost,
            .gap = gap_index,
            .stations = option.stations});
      }
    }
    if (statistics != nullptr) {
      statistics->fast_scoring_s += elapsed_seconds(scoring_started);
    }
  }
  if (deadline_reached(deadline)) return std::nullopt;
  if (tolerance_band.empty()) return std::nullopt;

  const auto materialization_started = std::chrono::steady_clock::now();
  const ReplacementCandidate* best = nullptr;
  std::vector<int> best_visits;
  for (const ReplacementCandidate& candidate : tolerance_band) {
    if (deadline_reached(deadline)) return std::nullopt;
    auto candidate_gaps = plan.gaps;
    candidate_gaps[candidate.gap] = candidate.stations;
    std::vector<int> candidate_visits =
        expand_route(plan.customer_ids, candidate_gaps).visits;
    if (best == nullptr || candidate_visits < best_visits) {
      best = &candidate;
      best_visits = std::move(candidate_visits);
    }
  }
  auto best_gaps = plan.gaps;
  best_gaps[best->gap] = best->stations;
  Plan result = factory_.make_plan(plan.customer_ids, std::move(best_gaps));
  if (statistics != nullptr) {
    ++statistics->full_materializations;
    statistics->full_materialization_s +=
        elapsed_seconds(materialization_started);
  }
  if (deadline_reached(deadline)) return std::nullopt;
  return result;
}

Plan PathSampler::improve_path_subset(
    Plan plan, double penalty_lambda,
    std::span<const std::size_t> eligible_gaps,
    std::optional<std::chrono::steady_clock::time_point> deadline,
    PathCompletionStatistics* statistics) {
  while (!eligible_gaps.empty() && !deadline_reached(deadline)) {
    std::vector<ExcludedPathOption> excluded;
    bool accepted = false;
    while (!deadline_reached(deadline)) {
      auto improved = best_replacement_in_gaps(
          plan, penalty_lambda, eligible_gaps, excluded, deadline,
          statistics);
      if (!improved.has_value()) return plan;
      const auto acceptance_started = std::chrono::steady_clock::now();
      const double incumbent_cost =
          generalized_cost(*plan.evaluation, instance_, penalty_lambda);
      const double improved_cost =
          generalized_cost(*improved->evaluation, instance_, penalty_lambda);
      if (improved_cost < incumbent_cost - kCostTolerance) {
        plan = std::move(*improved);
        accepted = true;
        if (statistics != nullptr) {
          ++statistics->accepted_steps;
          statistics->acceptance_s += elapsed_seconds(acceptance_started);
        }
        break;
      }

      // The fast nonlinear label and the full physical evaluator use the same
      // semantics but can differ at a floating-point tolerance boundary.  Do
      // not terminate the solver or hide a second valid option: exclude only
      // the rejected descriptor and resume the same deterministic ranking.
      std::optional<std::size_t> changed_gap;
      for (std::size_t gap = 0; gap < plan.gaps.size(); ++gap) {
        if (plan.gaps[gap] == improved->gaps[gap]) continue;
        if (changed_gap.has_value()) return plan;
        changed_gap = gap;
      }
      if (!changed_gap.has_value()) return plan;
      excluded.push_back(ExcludedPathOption{
          .gap = *changed_gap,
          .stations = improved->gaps[*changed_gap]});
      if (statistics != nullptr) {
        statistics->acceptance_s += elapsed_seconds(acceptance_started);
      }
    }
    if (!accepted) break;
  }
  return plan;
}

Plan PathSampler::repair_new_gaps(
    Plan plan, const PlanFactory::GapMap& inherited_gaps,
    double penalty_lambda,
    std::optional<std::chrono::steady_clock::time_point> deadline) {
  if (plan.exact_charging || deadline_reached(deadline)) return plan;
  std::vector<int> anchors;
  anchors.reserve(plan.customer_ids.size() + 2);
  anchors.push_back(instance_.depot().id);
  anchors.insert(anchors.end(), plan.customer_ids.begin(),
                 plan.customer_ids.end());
  anchors.push_back(instance_.depot().id);

  std::vector<std::size_t> introduced;
  introduced.reserve(plan.gaps.size());
  for (std::size_t gap = 0; gap < plan.gaps.size(); ++gap) {
    if (!inherited_gaps.contains(AnchorPair{anchors[gap], anchors[gap + 1]})) {
      introduced.push_back(gap);
    }
  }
  return improve_path_subset(std::move(plan), penalty_lambda, introduced,
                             deadline, nullptr);
}

Plan PathSampler::repair_affected_paths(
    Plan plan, const PlanFactory::GapMap& inherited_gaps,
    double penalty_lambda,
    std::optional<std::chrono::steady_clock::time_point> deadline,
    PathCompletionStatistics* statistics) {
  if (plan.exact_charging || deadline_reached(deadline)) return plan;
  const auto total_started = std::chrono::steady_clock::now();
  const auto eligibility_started = total_started;
  if (statistics != nullptr) ++statistics->calls;
  std::vector<int> anchors;
  anchors.reserve(plan.customer_ids.size() + 2);
  anchors.push_back(instance_.depot().id);
  anchors.insert(anchors.end(), plan.customer_ids.begin(),
                 plan.customer_ids.end());
  anchors.push_back(instance_.depot().id);

  std::vector<bool> selected(plan.gaps.size(), false);
  for (std::size_t gap = 0; gap < plan.gaps.size(); ++gap) {
    const bool introduced =
        !inherited_gaps.contains(AnchorPair{anchors[gap], anchors[gap + 1]});
    if (introduced || !plan.gaps[gap].empty()) selected[gap] = true;
    if (introduced && gap > 0) selected[gap - 1] = true;
    if (introduced && gap + 1 < plan.gaps.size()) selected[gap + 1] = true;
  }
  std::vector<std::size_t> eligible;
  eligible.reserve(plan.gaps.size());
  for (std::size_t gap = 0; gap < selected.size(); ++gap) {
    if (selected[gap]) eligible.push_back(gap);
  }
  if (statistics != nullptr) {
    statistics->eligible_gaps += eligible.size();
    statistics->eligibility_s += elapsed_seconds(eligibility_started);
  }
  Plan result = improve_path_subset(std::move(plan), penalty_lambda, eligible,
                                    deadline, statistics);
  if (statistics != nullptr) {
    statistics->total_s += elapsed_seconds(total_started);
  }
  return result;
}

InitialSolutionBuilder::InitialSolutionBuilder(const Instance& instance,
                                               PlanFactory& factory,
                                               PathSampler& path_sampler,
                                               std::mt19937_64& random,
                                               SearchProfile profile,
                                               std::size_t promising_arc_count,
                                               std::size_t granular_threshold)
    : instance_(instance),
      factory_(factory),
      path_sampler_(path_sampler),
      random_(random),
      profile_(profile),
      granular_(profile == SearchProfile::GranularAdaptive &&
                        instance.customer_ids().size() > granular_threshold
                    ? std::make_unique<GranularNeighborhood>(
                          instance, promising_arc_count)
                    : nullptr),
      maximum_charging_rate_wh_per_h_(
          completes_insertion_candidates(profile) ||
                  uses_granular_topology(profile)
              ? maximum_charging_rate(instance)
              : 0.0) {}

Solution InitialSolutionBuilder::build(double penalty_lambda) {
  std::vector<int> unassigned = instance_.customer_ids();
  const std::size_t first_index = random_index(unassigned, random_);
  const int first = unassigned[first_index];
  unassigned.erase(unassigned.begin() +
                   static_cast<std::ptrdiff_t>(first_index));
  Plan first_plan = factory_.plan_sequence({first});
  if (completes_construction_winner(profile_)) {
    first_plan = path_sampler_.repair_new_gaps(
        std::move(first_plan), PlanFactory::GapMap{}, penalty_lambda);
  }
  Solution current{{std::move(first_plan)}};
  InsertionCompletionMemo completion_memo{
      kInsertionCompletionMemoCapacity};
  VirtualInsertionCompletionMemo virtual_completion_memo{
      kVirtualInsertionCompletionMemoCapacity};
  std::unique_ptr<RouteSequenceScorer> sequence_scorer;
  while (!unassigned.empty()) {
    const PlanFactory::GapMap reference_gaps = factory_.gap_map(current);
    if (sequence_scorer == nullptr) {
      sequence_scorer = std::make_unique<RouteSequenceScorer>(
          instance_, current, reference_gaps);
    } else {
      sequence_scorer->reset(current, reference_gaps);
    }
    const RouteSequenceScorer& scorer = *sequence_scorer;
    std::optional<InsertionMove> best;
    const bool complete_all_candidates =
        completes_insertion_candidates(profile_);
    if (uses_granular_topology(profile_)) {
      std::vector<InsertionMove> candidates;
      std::vector<RejectedInsertionDescriptor> rejected;
      for (const int customer_id : unassigned) {
        auto moves = insertion_moves(
            current, customer_id, instance_, scorer, penalty_lambda, true,
            maximum_charging_rate_wh_per_h_, false, granular_.get(),
            &statistics_.candidates_rejected_by_granular,
            granular_ != nullptr ? &rejected : nullptr);
        candidates.insert(candidates.end(),
                          std::make_move_iterator(moves.begin()),
                          std::make_move_iterator(moves.end()));
      }
      if (granular_ != nullptr) {
        append_granular_escape(
            candidates, rejected, current, instance_, scorer,
            penalty_lambda, maximum_charging_rate_wh_per_h_, random_,
            statistics_);
      }
      best = choose_lazy_path_greedy(
          std::move(candidates), current, path_sampler_, penalty_lambda,
          statistics_, virtual_completion_memo, nullptr);
    } else if (complete_all_candidates) {
      std::vector<InsertionMove> candidates;
      for (const int customer_id : unassigned) {
        auto moves = insertion_moves(
            current, customer_id, instance_, scorer, penalty_lambda,
            true, maximum_charging_rate_wh_per_h_, false, granular_.get(),
            &statistics_.candidates_rejected_by_granular);
        candidates.insert(candidates.end(),
                          std::make_move_iterator(moves.begin()),
                          std::make_move_iterator(moves.end()));
      }
      best = choose_linked_greedy(
          std::move(candidates), current, instance_, factory_, path_sampler_,
          reference_gaps, penalty_lambda, statistics_, completion_memo);
    } else {
      for (const int customer_id : unassigned) {
        for (InsertionMove move : insertion_moves(
                 current, customer_id, instance_, scorer, penalty_lambda,
                 false, 0.0, false, granular_.get(),
                 &statistics_.candidates_rejected_by_granular)) {
          if (!best.has_value() || insertion_better(move, *best)) {
            best = std::move(move);
          }
        }
      }
    }
    const int inserted = best->customer_id;
    // An insertion changes customer order only.  The two gaps incident to the
    // inserted customer are absent from the current partial solution and are
    // therefore materialized as direct paths by plan_sequence().
    Plan selected =
        best->completion && !uses_granular_topology(profile_)
            ? best->completion->completed_plan
            : factory_.plan_sequence(insertion_sequence(current, *best),
                                     reference_gaps);
    if (profile_ == SearchProfile::LinkedWinner) {
      selected = path_sampler_.repair_new_gaps(
          std::move(selected), reference_gaps, penalty_lambda);
    }
    apply_insertion(current, best->route_index, std::move(selected));
    std::erase(unassigned, inserted);
  }
  current.validate_partition(instance_);
  return current;
}

Perturbation::Perturbation(const Instance& instance, PlanFactory& factory,
                           PathSampler& path_sampler,
                           std::mt19937_64& random,
                           SearchProfile profile,
                           std::size_t promising_arc_count,
                           std::size_t granular_threshold)
    : instance_(instance),
      factory_(factory),
      path_sampler_(path_sampler),
      random_(random),
      profile_(profile),
      granular_(uses_granular_topology(profile) &&
                        instance.customer_ids().size() >
                            granular_threshold
                    ? std::make_unique<GranularNeighborhood>(
                          instance, promising_arc_count)
                    : nullptr),
      maximum_charging_rate_wh_per_h_(
          completes_insertion_candidates(profile) ||
                  uses_granular_topology(profile)
              ? maximum_charging_rate(instance)
              : 0.0) {}

DestroyedSolution Perturbation::destroy(const Solution& solution,
                                        DestroyOperator operation,
                                        std::size_t removal_count,
                                        double penalty_lambda) {
  std::vector<int> customers = all_customers(solution);
  removal_count = std::min(removal_count, customers.size());
  std::vector<int> removed;
  removed.reserve(removal_count);

  if (operation == DestroyOperator::RandomRemoval) {
    std::shuffle(customers.begin(), customers.end(), random_);
    removed.assign(customers.begin(), customers.begin() +
                                          static_cast<std::ptrdiff_t>(removal_count));
  } else if (operation == DestroyOperator::RandomRoute) {
    std::vector<std::size_t> route_order(solution.plans.size());
    std::iota(route_order.begin(), route_order.end(), 0);
    std::shuffle(route_order.begin(), route_order.end(), random_);
    for (const std::size_t route_index : route_order) {
      const auto& route_customers = solution.plans[route_index].customer_ids;
      removed.insert(removed.end(), route_customers.begin(),
                     route_customers.end());
      if (removed.size() >= removal_count) break;
    }
  } else if (operation == DestroyOperator::ClosestRemoval) {
    const int seed = customers[random_index(customers, random_)];
    std::sort(customers.begin(), customers.end(), [&](int left, int right) {
      return std::tuple{instance_.distance_km(seed, left), left} <
             std::tuple{instance_.distance_km(seed, right), right};
    });
    removed.assign(customers.begin(), customers.begin() +
                                          static_cast<std::ptrdiff_t>(removal_count));
  } else {
    struct TargetContribution {
      double value{};
      int customer_id{};
      std::size_t route_index{};
    };
    std::vector<TargetContribution> contribution;
    const auto reference_gaps = factory_.gap_map(solution);
    if (sequence_scorer_ == nullptr) {
      sequence_scorer_ = std::make_unique<RouteSequenceScorer>(
          instance_, solution, reference_gaps);
    } else {
      sequence_scorer_->reset(solution, reference_gaps);
    }
    const RouteSequenceScorer& scorer = *sequence_scorer_;
    for (std::size_t route_index = 0; route_index < solution.plans.size();
         ++route_index) {
      const Plan& plan = solution.plans[route_index];
      const double current_cost =
          generalized_cost(*plan.evaluation, instance_, penalty_lambda);
      for (const int customer_id : plan.customer_ids) {
        std::vector<int> reduced = plan.customer_ids;
        std::erase(reduced, customer_id);
        const double reduced_cost =
            reduced.empty()
                ? 0.0
                : sequence_cost(reduced, scorer, penalty_lambda);
        contribution.push_back(TargetContribution{
            .value = current_cost - reduced_cost,
            .customer_id = customer_id,
            .route_index = route_index});
      }
    }
    std::sort(contribution.begin(), contribution.end(),
              [](const auto& left, const auto& right) {
                if (std::abs(left.value - right.value) > kCostTolerance)
                  return left.value > right.value;
                return std::tie(left.customer_id, left.route_index) <
                       std::tie(right.customer_id, right.route_index);
              });
    const int target = contribution.front().customer_id;
    const std::size_t target_route = contribution.front().route_index;
    removed.push_back(target);
    struct RankedRoute {
      double distance_to_target{};
      std::size_t route_index{};
    };
    std::vector<RankedRoute> ranked;
    for (std::size_t route_index = 0; route_index < solution.plans.size();
         ++route_index) {
      if (route_index == target_route) continue;
      const Plan& plan = solution.plans[route_index];
      double route_distance = std::numeric_limits<double>::infinity();
      for (const int customer_id : plan.customer_ids) {
        route_distance = std::min(
            route_distance, instance_.distance_km(target, customer_id));
      }
      ranked.push_back(RankedRoute{.distance_to_target = route_distance,
                                   .route_index = route_index});
    }
    std::sort(ranked.begin(), ranked.end(), [](const auto& left,
                                               const auto& right) {
      return std::tie(left.distance_to_target, left.route_index) <
             std::tie(right.distance_to_target, right.route_index);
    });
    for (const RankedRoute& candidate : ranked) {
      const auto& route_customers =
          solution.plans[candidate.route_index].customer_ids;
      removed.insert(removed.end(), route_customers.begin(),
                     route_customers.end());
      if (removed.size() >= removal_count) break;
    }
    // A single-route solution has no other route to ruin.  In that boundary
    // case, remove the remainder of the target route to satisfy the quota.
    if (removed.size() < removal_count) {
      for (const int customer_id :
           solution.plans[target_route].customer_ids) {
        if (customer_id != target) removed.push_back(customer_id);
      }
    }
  }

  auto sequences = remove_customers(solution, removed);
  Solution remaining =
      build_solution_from_sequences(sequences, factory_, &solution);
  return DestroyedSolution{.remaining = std::move(remaining),
                           .removed_customers = std::move(removed)};
}

Solution Perturbation::repair(DestroyedSolution destroyed,
                              RepairOperator operation,
                              double penalty_lambda,
                              std::vector<InsertionPathHint>* hints) {
  Solution current = std::move(destroyed.remaining);
  std::vector<int> unassigned = std::move(destroyed.removed_customers);
  InsertionCompletionMemo completion_memo{
      kInsertionCompletionMemoCapacity};
  VirtualInsertionCompletionMemo virtual_completion_memo{
      kVirtualInsertionCompletionMemoCapacity};
  while (!unassigned.empty()) {
    // Only paths that still exist in the current partial solution are
    // inherited.  Every insertion creates two direct gaps, even when the
    // customer previously occupied the same position before destruction.
    const PlanFactory::GapMap reference_gaps = factory_.gap_map(current);
    if (sequence_scorer_ == nullptr) {
      sequence_scorer_ = std::make_unique<RouteSequenceScorer>(
          instance_, current, reference_gaps);
    } else {
      sequence_scorer_->reset(current, reference_gaps);
    }
    const RouteSequenceScorer& scorer = *sequence_scorer_;
    std::optional<InsertionMove> chosen;
    const bool lazy_path_repair = uses_granular_topology(profile_);
    if (lazy_path_repair) {
      const auto direct_scoring_started = std::chrono::steady_clock::now();
      std::vector<InsertionMove> candidates;
      std::vector<RejectedInsertionDescriptor> rejected;
      for (const int customer_id : unassigned) {
        auto moves = insertion_moves(
            current, customer_id, instance_, scorer, penalty_lambda, true,
            maximum_charging_rate_wh_per_h_, false, granular_.get(),
            &statistics_.candidates_rejected_by_granular,
            granular_ != nullptr ? &rejected : nullptr);
        candidates.insert(candidates.end(),
                          std::make_move_iterator(moves.begin()),
                          std::make_move_iterator(moves.end()));
      }
      if (granular_ != nullptr) {
        append_granular_escape(
            candidates, rejected, current, instance_, scorer,
            penalty_lambda, maximum_charging_rate_wh_per_h_, random_,
            statistics_);
      }
      statistics_.repair_timing.granular_direct_scoring_s +=
          elapsed_seconds(direct_scoring_started);
      if (operation == RepairOperator::GreedyInsertion) {
        chosen = choose_lazy_path_greedy(
            std::move(candidates), current, path_sampler_, penalty_lambda,
            statistics_, virtual_completion_memo,
            &statistics_.repair_timing);
      } else {
        chosen = choose_lazy_path_two_regret(
            std::move(candidates), current, unassigned, path_sampler_,
            penalty_lambda, statistics_, virtual_completion_memo,
            &statistics_.repair_timing);
      }
    } else if (completes_insertion_candidates(profile_) &&
        operation == RepairOperator::GreedyInsertion) {
      std::vector<InsertionMove> candidates;
      for (const int customer_id : unassigned) {
        auto moves = insertion_moves(
            current, customer_id, instance_, scorer, penalty_lambda,
            true, maximum_charging_rate_wh_per_h_, false, granular_.get(),
            &statistics_.candidates_rejected_by_granular);
        candidates.insert(candidates.end(),
                          std::make_move_iterator(moves.begin()),
                          std::make_move_iterator(moves.end()));
      }
      chosen = choose_linked_greedy(
          std::move(candidates), current, instance_, factory_, path_sampler_,
          reference_gaps, penalty_lambda, statistics_, completion_memo);
    } else if (completes_insertion_candidates(profile_)) {
      chosen = choose_linked_two_regret(
          current, unassigned, instance_, factory_, path_sampler_, scorer,
          reference_gaps, penalty_lambda, maximum_charging_rate_wh_per_h_,
          statistics_, completion_memo);
    } else if (operation == RepairOperator::GreedyInsertion) {
      for (const int customer_id : unassigned) {
        auto moves = insertion_moves(
            current, customer_id, instance_, scorer, penalty_lambda,
            false, 0.0, false,
            granular_.get(),
            &statistics_.candidates_rejected_by_granular);
        for (InsertionMove move : std::move(moves)) {
          if (!chosen.has_value() || insertion_better(move, *chosen)) {
            chosen = std::move(move);
          }
        }
      }
    } else {
      double best_regret = -1.0;
      for (const int customer_id : unassigned) {
        auto moves = insertion_moves(current, customer_id, instance_, scorer,
                                     penalty_lambda, false, 0.0, false,
                                     granular_.get(),
                                     &statistics_.candidates_rejected_by_granular);
        std::sort(moves.begin(), moves.end(), insertion_better);
        const double regret = moves.size() == 1
                                  ? std::numeric_limits<double>::infinity()
                                  : moves[1].delta - moves[0].delta;
        if (!chosen.has_value() || regret > best_regret + kCostTolerance ||
            (std::abs(regret - best_regret) <= kCostTolerance &&
             insertion_better(moves[0], *chosen))) {
          best_regret = regret;
          chosen = std::move(moves[0]);
        }
      }
    }
    const auto materialization_started = std::chrono::steady_clock::now();
    const int inserted = chosen->customer_id;
    if (hints != nullptr && lazy_path_repair) {
      update_insertion_hints(*hints, current, *chosen, instance_);
    }
    Plan selected =
        chosen->completion && !lazy_path_repair
            ? chosen->completion->completed_plan
            : factory_.plan_sequence(insertion_sequence(current, *chosen),
                                     reference_gaps);
    if (profile_ == SearchProfile::LinkedWinner) {
      selected = path_sampler_.repair_new_gaps(
          std::move(selected), reference_gaps, penalty_lambda);
    }
    apply_insertion(current, chosen->route_index, std::move(selected));
    std::erase(unassigned, inserted);
    statistics_.repair_timing.winner_materialization_s +=
        elapsed_seconds(materialization_started);
  }
  current.validate_partition(instance_);
  return current;
}

Solution Perturbation::apply(const Solution& solution, double penalty_lambda) {
  const std::size_t customers = instance_.customer_ids().size();
  const std::size_t lower = std::min<std::size_t>(customers, 3);
  const std::size_t upper = std::max(
      lower, static_cast<std::size_t>(std::floor(std::sqrt(customers))));
  const std::size_t removal_count =
      std::uniform_int_distribution<std::size_t>(lower, upper)(random_);
  const auto destroy_operation = static_cast<DestroyOperator>(
      std::uniform_int_distribution<int>(0, 3)(random_));
  const auto repair_operation = static_cast<RepairOperator>(
      std::uniform_int_distribution<int>(0, 1)(random_));
  std::vector<InsertionPathHint> hints;
  const auto destroy_started = std::chrono::steady_clock::now();
  DestroyedSolution destroyed =
      destroy(solution, destroy_operation, removal_count, penalty_lambda);
  statistics_.repair_timing.destroy_s += elapsed_seconds(destroy_started);
  Solution repaired = repair(
      std::move(destroyed), repair_operation, penalty_lambda,
      uses_granular_topology(profile_) ? &hints : nullptr);
  if (!uses_granular_topology(profile_) || hints.empty()) return repaired;
  const auto post_hint_started = std::chrono::steady_clock::now();
  HintedPathResult hinted = path_sampler_.apply_insertion_hints(
      std::move(repaired), hints, penalty_lambda);
  statistics_.repair_timing.post_hint_s += elapsed_seconds(post_hint_started);
  statistics_.post_repair_hint_attempts += hinted.attempts;
  statistics_.post_repair_hint_accepted += hinted.accepted;
  return std::move(hinted.solution);
}

VariableNeighborhoodDescent::VariableNeighborhoodDescent(
    const Instance& instance, PlanFactory& factory, PathSampler& path_sampler,
    std::mt19937_64& random, std::size_t promising_arc_count,
    SearchProfile profile)
    : instance_(instance),
      factory_(factory),
      granular_(instance, promising_arc_count),
      path_sampler_(path_sampler),
      random_(random),
      profile_(profile) {}

bool VariableNeighborhoodDescent::expired(
    std::optional<std::chrono::steady_clock::time_point> deadline) const {
  return deadline.has_value() && std::chrono::steady_clock::now() >= *deadline;
}

Solution VariableNeighborhoodDescent::improve(
    Solution solution, double penalty_lambda,
    std::optional<std::chrono::steady_clock::time_point> deadline) {
  move_memory_.clear();
  route_versions_.clear();
  next_route_version_ = 1;

  const auto accept = [&](std::optional<Solution> candidate) {
    if (!candidate.has_value() ||
        generalized_cost(*candidate, instance_, penalty_lambda) >=
            generalized_cost(solution, instance_, penalty_lambda) -
                kCostTolerance) {
      return false;
    }
    ++statistics_.accepted_moves;
    const std::size_t previous_route_count = solution.plans.size();
    solution = std::move(*candidate);
    if (solution.plans.size() != previous_route_count) {
      move_memory_.clear();
    }
    return true;
  };

  if (uses_adaptive_vnd(profile_)) {
    constexpr std::array customer_order{
        Neighborhood::Swap, Neighborhood::Relocate,
        Neighborhood::TwoOptStar};
    for (const Neighborhood neighborhood : customer_order) {
      std::unordered_set<AnchorPair, AnchorPairHash> dirty_pairs;
      while (!expired(deadline)) {
        auto best = best_neighbor(solution, neighborhood, penalty_lambda,
                                  deadline);
        if (expired(deadline)) break;
        const PlanFactory::GapMap previous_gaps = factory_.gap_map(solution);
        if (!accept(std::move(best))) break;
        accumulate_dirty_direct_anchor_pairs(previous_gaps, solution,
                                             instance_, dirty_pairs);
      }
      if (expired(deadline)) break;
      if (std::bernoulli_distribution(0.5)(random_)) {
        ++statistics_.adaptive_replace_path_selections;
        const std::vector<AnchorPair> stage_dirty(dirty_pairs.begin(),
                                                  dirty_pairs.end());
        solution = inserted_replace_path(std::move(solution), penalty_lambda,
                                         stage_dirty, deadline);
      }
      ++statistics_.adaptive_customer_advances;
    }
    if (!expired(deadline)) {
      solution = final_replace_path(std::move(solution), penalty_lambda,
                                    deadline);
    }
    solution.validate_partition(instance_);
    return solution;
  }

  std::size_t neighborhood_index = 0;
  while (neighborhood_index < order_.size() && !expired(deadline)) {
    auto best = best_neighbor(solution, order_[neighborhood_index],
                              penalty_lambda, deadline);
    // A linked completion may have reached the absolute deadline while
    // constructing a temporary endpoint.  Do not publish any partial
    // completion after the boundary.
    if (expired(deadline)) break;
    // Customer moves are ranked by the same incremental direct-gap scorer in
    // every profile.  LinkedWinner may deterministically improve only the
    // selected winner's newly-created gaps before this full acceptance check.
    if (accept(std::move(best))) {
      neighborhood_index = 0;
    } else {
      ++neighborhood_index;
    }
  }
  solution.validate_partition(instance_);
  return solution;
}

std::optional<Solution> VariableNeighborhoodDescent::best_neighbor(
    const Solution& solution, Neighborhood neighborhood, double penalty_lambda,
    std::optional<std::chrono::steady_clock::time_point> deadline) {
  if (neighborhood == Neighborhood::ReplacePath) {
    return replace_path_neighbor(solution, penalty_lambda);
  }
  return best_customer_neighbor(solution, neighborhood, penalty_lambda,
                                deadline);
}

std::optional<Solution> VariableNeighborhoodDescent::best_customer_neighbor(
    const Solution& solution, Neighborhood neighborhood, double penalty_lambda,
    std::optional<std::chrono::steady_clock::time_point> deadline) {
  ++statistics_.customer_neighborhood_calls;
  const auto current = solution.customer_sequences();
  const auto reference_gaps = factory_.gap_map(solution);
  if (sequence_scorer_ == nullptr) {
    sequence_scorer_ = std::make_unique<RouteSequenceScorer>(
        instance_, solution, reference_gaps);
  } else {
    sequence_scorer_->reset(solution, reference_gaps);
  }
  const RouteSequenceScorer& scorer = *sequence_scorer_;
  auto prepared_plan = [&](const std::vector<int>& sequence) {
    auto gaps = direct_gaps(sequence.size());
    std::vector<int> anchors;
    anchors.reserve(sequence.size() + 2);
    anchors.push_back(instance_.depot().id);
    anchors.insert(anchors.end(), sequence.begin(), sequence.end());
    anchors.push_back(instance_.depot().id);
    for (std::size_t gap = 0; gap + 1 < anchors.size(); ++gap) {
      const auto known =
          reference_gaps.find(AnchorPair{anchors[gap], anchors[gap + 1]});
      if (known != reference_gaps.end()) gaps[gap] = known->second;
    }
    return factory_.make_plan(sequence, std::move(gaps));
  };
  std::vector<double> current_route_costs;
  current_route_costs.reserve(solution.plans.size());
  for (const Plan& plan : solution.plans) {
    current_route_costs.push_back(
        generalized_cost(*plan.evaluation, instance_, penalty_lambda));
  }

  std::vector<std::uint64_t> current_route_versions;
  current_route_versions.reserve(solution.plans.size());
  for (const Plan& plan : solution.plans) {
    auto [version, inserted] = route_versions_.try_emplace(
        plan.route.visits, next_route_version_);
    if (inserted) ++next_route_version_;
    current_route_versions.push_back(version->second);
  }

  struct CustomerLocation {
    std::size_t route{};
    std::size_t position{};
    std::size_t customer_slot{};
  };
  const std::size_t customer_count = instance_.customer_ids().size();
  std::vector<CustomerLocation> location(customer_count);
  std::vector<std::vector<std::size_t>> current_ranks;
  current_ranks.reserve(current.size());
  for (std::size_t route = 0; route < current.size(); ++route) {
    std::vector<std::size_t> route_ranks;
    route_ranks.reserve(current[route].size());
    for (std::size_t position = 0; position < current[route].size();
         ++position) {
      const std::size_t rank =
          granular_.customer_rank(current[route][position]);
      location[rank] = CustomerLocation{route, position, rank};
      route_ranks.push_back(rank);
    }
    current_ranks.push_back(std::move(route_ranks));
  }
  const auto customer_location = [&](int customer_id) -> CustomerLocation {
    return location[granular_.customer_rank(customer_id)];
  };
  const auto is_current_arc = [&](std::size_t origin_rank,
                                  std::size_t target_rank) {
    const CustomerLocation& origin = location[origin_rank];
    const CustomerLocation& target = location[target_rank];
    return origin.route == target.route &&
           origin.position + 1 == target.position;
  };

  // Forward slice arcs already belong to the incumbent.  A reversed slice can
  // introduce promising arcs internally, so retain a prefix count for those
  // arcs and make the final granular guard constant in the number of slices.
  std::vector<std::vector<std::size_t>> reverse_promising_prefix;
  reverse_promising_prefix.reserve(current.size());
  for (const auto& sequence : current_ranks) {
    std::vector<std::size_t> prefix(sequence.size(), 0);
    for (std::size_t position = 0; position + 1 < sequence.size();
         ++position) {
      prefix[position + 1] =
          prefix[position] +
          static_cast<std::size_t>(
              !is_current_arc(sequence[position + 1], sequence[position]) &&
              granular_.promising_by_rank(sequence[position + 1],
                                           sequence[position]));
    }
    reverse_promising_prefix.push_back(std::move(prefix));
  }

  const auto add_slice = [&](DescribedRouteChange& change,
                             std::size_t source_route, std::size_t begin,
                             std::size_t end, bool reversed = false) {
    if (begin == end) return;
    if (source_route >= current.size() || begin > end ||
        end > current[source_route].size()) {
      throw std::out_of_range("move descriptor has invalid route slice");
    }
    if (change.slice_count >= change.slices.size()) {
      throw std::logic_error("move descriptor exceeded its slice bound");
    }
    change.slices[change.slice_count++] = IncumbentRouteSlice{
        .route_index = source_route,
        .begin = begin,
        .end = end,
        .reversed = reversed};
  };

  const auto one_change = [](DescribedRouteChange change) {
    MoveDescriptor move;
    move.changes[0] = std::move(change);
    move.change_count = 1;
    return move;
  };
  const auto two_changes = [](DescribedRouteChange first,
                              DescribedRouteChange second) {
    if (first.route == second.route) {
      throw std::logic_error("move descriptor changes one route twice");
    }
    if (second.route < first.route) std::swap(first, second);
    MoveDescriptor move;
    move.changes[0] = std::move(first);
    move.changes[1] = std::move(second);
    move.change_count = 2;
    return move;
  };

  const auto materialize_route = [&](const DescribedRouteChange& change) {
    std::size_t route_customer_count = 0;
    for (std::size_t index = 0; index < change.slice_count; ++index) {
      route_customer_count += change.slices[index].end -
                              change.slices[index].begin;
    }
    std::vector<int> sequence;
    sequence.reserve(route_customer_count);
    for (std::size_t index = 0; index < change.slice_count; ++index) {
      const IncumbentRouteSlice& slice = change.slices[index];
      if (!slice.reversed) {
        sequence.insert(
            sequence.end(),
            current[slice.route_index].begin() +
                static_cast<std::ptrdiff_t>(slice.begin),
            current[slice.route_index].begin() +
                static_cast<std::ptrdiff_t>(slice.end));
      } else {
        for (std::size_t position = slice.end; position > slice.begin;
             --position) {
          sequence.push_back(current[slice.route_index][position - 1]);
        }
      }
    }
    return sequence;
  };

  using RouteChange = std::pair<std::size_t, std::vector<int>>;
  using RouteChanges = std::vector<RouteChange>;
  const auto materialize_changes = [&](const MoveDescriptor& move) {
    RouteChanges changes;
    changes.reserve(move.change_count);
    for (std::size_t index = 0; index < move.change_count; ++index) {
      changes.emplace_back(move.changes[index].route,
                           materialize_route(move.changes[index]));
    }
    return changes;
  };
  const auto move_less = [&](const MoveDescriptor& left,
                             const MoveDescriptor& right) {
    return materialize_changes(left) < materialize_changes(right);
  };

  const auto first_customer_rank = [&](const IncumbentRouteSlice& slice) {
    return slice.reversed
               ? current_ranks[slice.route_index][slice.end - 1]
               : current_ranks[slice.route_index][slice.begin];
  };
  const auto last_customer_rank = [&](const IncumbentRouteSlice& slice) {
    return slice.reversed
               ? current_ranks[slice.route_index][slice.begin]
               : current_ranks[slice.route_index][slice.end - 1];
  };
  const auto move_allowed = [&](const MoveDescriptor& move) {
    if (granular_.gamma_size() == 0) return true;
    for (std::size_t change_index = 0;
         change_index < move.change_count; ++change_index) {
      const DescribedRouteChange& change = move.changes[change_index];
      std::optional<std::size_t> previous;
      for (std::size_t slice_index = 0;
           slice_index < change.slice_count; ++slice_index) {
        const IncumbentRouteSlice& slice = change.slices[slice_index];
        if (slice.reversed && slice.end > slice.begin + 1) {
          const auto& prefix =
              reverse_promising_prefix[slice.route_index];
          if (prefix[slice.end - 1] > prefix[slice.begin]) return true;
        }
        const std::size_t first = first_customer_rank(slice);
        if (previous.has_value() &&
            !is_current_arc(*previous, first) &&
            granular_.promising_by_rank(*previous, first)) {
          return true;
        }
        previous = last_customer_rank(slice);
      }
    }
    return false;
  };

  const std::size_t route_count = current.size();
  const std::size_t support_count =
      route_count * route_count + 2 * route_count;
  if (++support_epoch_ == 0) {
    for (SupportWorkingMemory& memory : support_working_memory_) {
      memory.epoch = 0;
    }
    ++support_epoch_;
  }
  if (support_working_memory_.size() < support_count) {
    support_working_memory_.resize(support_count);
  }
  const auto support_memory = [&](std::size_t index)
      -> SupportWorkingMemory& {
    SupportWorkingMemory& memory = support_working_memory_[index];
    if (memory.epoch != support_epoch_) {
      memory.epoch = support_epoch_;
      memory.state = SupportState::Unknown;
      memory.active = false;
      memory.entry.candidate_cost =
          std::numeric_limits<double>::infinity();
      memory.entry.move.change_count = 0;
    }
    return memory;
  };

  const auto support_index = [&](std::size_t first,
                                 std::optional<std::size_t> second,
                                 bool opens_route) {
    if (first >= route_count ||
        (second.has_value() && *second >= route_count) ||
        (opens_route && second.has_value())) {
      throw std::out_of_range("move support contains an invalid route");
    }
    if (opens_route) return route_count + first;
    if (!second.has_value()) return first;
    const auto [lower, upper] = std::minmax(first, *second);
    return 2 * route_count + lower * route_count + upper;
  };

  const auto support_key = [&](std::size_t first,
                               std::optional<std::size_t> second,
                               bool opens_route) {
    if (first >= current.size() ||
        (second.has_value() && *second >= current.size())) {
      throw std::out_of_range("move support contains an invalid route");
    }
    if (second.has_value() && *second < first) std::swap(first, *second);
    MoveMemoryKey key{.neighborhood = neighborhood,
                      .route_indices = {},
                      .route_versions = {},
                      .route_count = 1,
                      .opens_route = opens_route};
    key.route_indices[0] = first;
    key.route_versions[0] = current_route_versions[first];
    if (second.has_value()) {
      key.route_indices[key.route_count] = *second;
      key.route_versions[key.route_count] = current_route_versions[*second];
      ++key.route_count;
    }
    return key;
  };

  // Vidal's route-pair move memory is useful only if it is consulted before
  // rebuilding every customer-level move on unchanged routes.  Remember the
  // first lookup for each support in this neighborhood pass, then skip all of
  // its descriptors when the cached best move is still valid.
  const auto support_is_closed = [&](std::size_t first,
                                     std::optional<std::size_t> second,
                                     bool opens_route) {
    const std::size_t index = support_index(first, second, opens_route);
    SupportWorkingMemory& memory = support_memory(index);
    if (memory.state == SupportState::Unknown) {
      MoveMemoryKey key = support_key(first, second, opens_route);
      if (auto cached = move_memory_.get(key); cached.has_value()) {
        memory.active = true;
        memory.key = std::move(key);
        memory.entry = std::move(*cached);
        ++statistics_.route_support_groups;
        ++statistics_.route_support_cache_hits;
        ++statistics_.route_support_buckets_skipped;
        memory.state = SupportState::Closed;
      } else {
        memory.state = SupportState::Open;
      }
    }
    return memory.state == SupportState::Closed;
  };

  struct MoveSupport {
    std::size_t first{};
    std::optional<std::size_t> second;
    bool opens_route{};
  };
  const auto move_support = [&](const MoveDescriptor& move) {
    if (move.change_count == 0 || move.change_count > 2 ||
        move.changes[0].route >= route_count) {
      throw std::logic_error("move descriptor has invalid route support");
    }
    MoveSupport support{.first = move.changes[0].route,
                        .second = std::nullopt,
                        .opens_route = false};
    if (move.change_count == 2) {
      const std::size_t second = move.changes[1].route;
      if (second == route_count) {
        support.opens_route = true;
      } else if (second < route_count) {
        support.second = second;
      } else {
        throw std::logic_error("move descriptor opens an invalid route");
      }
    }
    return support;
  };

  struct ReusableRelocateSource {
    std::size_t route{};
    std::size_t customer_slot{};
  };
  std::vector<std::optional<RouteCostSummary>> relocate_source_summaries(
      neighborhood == Neighborhood::Relocate ? customer_count : 0);
  auto consider_move = [&](MoveDescriptor move,
                           std::optional<ReusableRelocateSource> reusable) {
    ++statistics_.descriptors_considered;
    if (expired(deadline)) return;
    const MoveSupport support = move_support(move);
    const std::size_t support_slot = support_index(
        support.first, support.second, support.opens_route);
    SupportWorkingMemory& memory = support_memory(support_slot);
    if (memory.state == SupportState::Closed) {
      ++statistics_.descriptors_skipped_by_memory;
      return;
    }
    if (memory.state == SupportState::Unknown) {
      throw std::logic_error("move support was not prepared before scoring");
    }
    if (!memory.active) {
      memory.active = true;
      memory.key =
          support_key(support.first, support.second, support.opens_route);
      ++statistics_.route_support_groups;
    }
    if (!move_allowed(move)) {
      ++statistics_.granular_rejections;
      return;
    }

    double candidate_cost = 0.0;
    for (std::size_t index = 0; index < move.change_count; ++index) {
      const DescribedRouteChange& change = move.changes[index];
      if (change.slice_count == 0) continue;
      const auto slices = std::span<const IncumbentRouteSlice>{
          change.slices.data(), change.slice_count};
      if (reusable.has_value() && change.route == reusable->route) {
        auto& summary =
            relocate_source_summaries.at(reusable->customer_slot);
        if (!summary.has_value()) {
          summary = scorer.evaluate(slices);
          ++statistics_.slice_routes_evaluated;
          ++statistics_.relocate_source_summary_misses;
        } else {
          ++statistics_.relocate_source_summary_hits;
        }
        candidate_cost +=
            generalized_cost(*summary, instance_, penalty_lambda);
      } else {
        ++statistics_.slice_routes_evaluated;
        candidate_cost += scorer.generalized_cost(slices, penalty_lambda);
      }
    }
    ++statistics_.moves_evaluated;
    if (candidate_cost < memory.entry.candidate_cost - kCostTolerance ||
        (std::abs(candidate_cost - memory.entry.candidate_cost) <=
             kCostTolerance &&
         move_less(move, memory.entry.move))) {
      memory.entry = MoveMemoryEntry{.candidate_cost = candidate_cost,
                                     .move = std::move(move)};
    }
  };
  const auto consider = [&](MoveDescriptor move) {
    consider_move(std::move(move), std::nullopt);
  };

  if (neighborhood == Neighborhood::Swap) {
    std::vector<std::pair<int, int>>& customer_pairs = swap_customer_pairs_;
    customer_pairs.clear();
    customer_pairs.reserve(granular_.promising_pairs().size() * 4);
    const auto generate_swap_bucket = [&](std::size_t first_route,
                                          std::size_t second_route) {
      const auto add_from_route = [&](std::size_t source_route,
                                      std::size_t target_route) {
        for (std::size_t position = 0;
             position < current[source_route].size(); ++position) {
          const int candidate_customer = current[source_route][position];
          const auto add_from_adjacent = [&](int adjacent_customer) {
            for (const int promising_customer :
                 granular_.promising_neighbors_of(adjacent_customer)) {
              if (customer_location(promising_customer).route != target_route ||
                  candidate_customer == promising_customer) {
                continue;
              }
              customer_pairs.push_back(
                  std::minmax(candidate_customer, promising_customer));
            }
          };
          if (position > 0) {
            add_from_adjacent(current[source_route][position - 1]);
          }
          if (position + 1 < current[source_route].size()) {
            add_from_adjacent(current[source_route][position + 1]);
          }
        }
      };
      add_from_route(first_route, second_route);
      if (first_route != second_route) {
        add_from_route(second_route, first_route);
      }
    };
    for (std::size_t first_route = 0; first_route < route_count;
         ++first_route) {
      if (!support_is_closed(first_route, std::nullopt, false)) {
        generate_swap_bucket(first_route, first_route);
      }
      for (std::size_t second_route = first_route + 1;
           second_route < route_count; ++second_route) {
        if (!support_is_closed(first_route, second_route, false)) {
          generate_swap_bucket(first_route, second_route);
        }
      }
    }
    std::sort(customer_pairs.begin(), customer_pairs.end());
    customer_pairs.erase(
        std::unique(customer_pairs.begin(), customer_pairs.end()),
        customer_pairs.end());
    for (const auto& [left_customer, right_customer] : customer_pairs) {
      const CustomerLocation left = customer_location(left_customer);
      const CustomerLocation right = customer_location(right_customer);
      if (left.route == right.route) {
        const std::size_t low = std::min(left.position, right.position);
        const std::size_t high = std::max(left.position, right.position);
        if (high < low + 2) continue;
        DescribedRouteChange changed{.route = left.route};
        add_slice(changed, left.route, 0, low);
        add_slice(changed, left.route, high, high + 1);
        add_slice(changed, left.route, low + 1, high);
        add_slice(changed, left.route, low, low + 1);
        add_slice(changed, left.route, high + 1,
                  current[left.route].size());
        consider(one_change(std::move(changed)));
      } else {
        DescribedRouteChange changed_left{.route = left.route};
        add_slice(changed_left, left.route, 0, left.position);
        add_slice(changed_left, right.route, right.position,
                  right.position + 1);
        add_slice(changed_left, left.route, left.position + 1,
                  current[left.route].size());
        DescribedRouteChange changed_right{.route = right.route};
        add_slice(changed_right, right.route, 0, right.position);
        add_slice(changed_right, left.route, left.position,
                  left.position + 1);
        add_slice(changed_right, right.route, right.position + 1,
                  current[right.route].size());
        consider(two_changes(std::move(changed_left),
                             std::move(changed_right)));
      }
    }
  } else if (neighborhood == Neighborhood::Relocate) {
    using RelocateMove = std::tuple<int, std::size_t, std::size_t>;
    std::vector<RelocateMove>& moves = relocate_moves_;
    moves.clear();
    moves.reserve(instance_.customer_ids().size() *
                  (4 * granular_.gamma_size() + 1));
    const auto generate_relocate_bucket = [&](std::size_t source_route,
                                              std::size_t target_route) {
      for (const int customer : current[source_route]) {
        for (const int neighbor :
             granular_.promising_neighbors_of(customer)) {
          const CustomerLocation target = customer_location(neighbor);
          if (target.route != target_route) continue;
          moves.emplace_back(customer, target_route, target.position);
          moves.emplace_back(customer, target_route, target.position + 1);
        }
      }
    };
    for (std::size_t first_route = 0; first_route < route_count;
         ++first_route) {
      if (!support_is_closed(first_route, std::nullopt, false)) {
        generate_relocate_bucket(first_route, first_route);
      }
      if (current[first_route].size() > 1 &&
          !support_is_closed(first_route, std::nullopt, true)) {
        for (const int customer : current[first_route]) {
          moves.emplace_back(customer, current.size(), 0);
        }
      }
      for (std::size_t second_route = first_route + 1;
           second_route < route_count; ++second_route) {
        if (support_is_closed(first_route, second_route, false)) continue;
        generate_relocate_bucket(first_route, second_route);
        generate_relocate_bucket(second_route, first_route);
      }
    }
    std::sort(moves.begin(), moves.end());
    moves.erase(std::unique(moves.begin(), moves.end()), moves.end());
    for (const auto& [customer, target_route, original_gap] : moves) {
      const CustomerLocation source = customer_location(customer);
      if (target_route == source.route) {
        std::size_t insertion = original_gap;
        if (source.position < original_gap) --insertion;
        if (insertion == source.position) continue;
        DescribedRouteChange changed{.route = source.route};
        if (source.position < insertion) {
          add_slice(changed, source.route, 0, source.position);
          add_slice(changed, source.route, source.position + 1,
                    original_gap);
          add_slice(changed, source.route, source.position,
                    source.position + 1);
          add_slice(changed, source.route, original_gap,
                    current[source.route].size());
        } else {
          add_slice(changed, source.route, 0, insertion);
          add_slice(changed, source.route, source.position,
                    source.position + 1);
          add_slice(changed, source.route, insertion, source.position);
          add_slice(changed, source.route, source.position + 1,
                    current[source.route].size());
        }
        consider(one_change(std::move(changed)));
      } else if (target_route == current.size()) {
        DescribedRouteChange changed_source{.route = source.route};
        add_slice(changed_source, source.route, 0, source.position);
        add_slice(changed_source, source.route, source.position + 1,
                  current[source.route].size());
        DescribedRouteChange singleton{.route = current.size()};
        add_slice(singleton, source.route, source.position,
                  source.position + 1);
        consider_move(
            two_changes(std::move(changed_source), std::move(singleton)),
            ReusableRelocateSource{.route = source.route,
                                   .customer_slot = source.customer_slot});
      } else {
        DescribedRouteChange changed_source{.route = source.route};
        add_slice(changed_source, source.route, 0, source.position);
        add_slice(changed_source, source.route, source.position + 1,
                  current[source.route].size());
        DescribedRouteChange changed_target{.route = target_route};
        add_slice(changed_target, target_route, 0, original_gap);
        add_slice(changed_target, source.route, source.position,
                  source.position + 1);
        add_slice(changed_target, target_route, original_gap,
                  current[target_route].size());
        consider_move(
            two_changes(std::move(changed_source), std::move(changed_target)),
            ReusableRelocateSource{.route = source.route,
                                   .customer_slot = source.customer_slot});
      }
    }
  } else if (neighborhood == Neighborhood::TwoOptStar) {
    using TwoOptMove =
        std::tuple<std::size_t, std::size_t, std::size_t, std::size_t, bool>;
    std::vector<TwoOptMove>& moves = two_opt_moves_;
    moves.clear();
    moves.reserve(granular_.promising_pairs().size() * 4 +
                  current.size() * (current.size() - 1) / 2);
    for (std::size_t left_route = 0; left_route < current.size(); ++left_route) {
      for (std::size_t right_route = left_route + 1;
           right_route < current.size(); ++right_route) {
        if (support_is_closed(left_route, right_route, false)) continue;
        for (std::size_t left_position = 0;
             left_position < current[left_route].size(); ++left_position) {
          const int left_customer = current[left_route][left_position];
          for (const int right_customer :
               granular_.promising_neighbors_of(left_customer)) {
            const CustomerLocation right =
                customer_location(right_customer);
            if (right.route != right_route) continue;
            moves.emplace_back(left_route, right_route, left_position + 1,
                               right.position, false);
            moves.emplace_back(left_route, right_route, left_position,
                               right.position + 1, false);
            moves.emplace_back(left_route, right_route, left_position + 1,
                               right.position + 1, true);
            moves.emplace_back(left_route, right_route, left_position,
                               right.position, true);
          }
        }
        moves.emplace_back(left_route, right_route, 0,
                           current[right_route].size(), true);
      }
    }
    std::sort(moves.begin(), moves.end());
    moves.erase(std::unique(moves.begin(), moves.end()), moves.end());
    for (const auto& [left_route, right_route, left_cut, right_cut,
                      reverse_reconnection] : moves) {
      DescribedRouteChange changed_left{.route = left_route};
      DescribedRouteChange changed_right{.route = right_route};
      if (!reverse_reconnection) {
        add_slice(changed_left, left_route, 0, left_cut);
        add_slice(changed_left, right_route, right_cut,
                  current[right_route].size());
        add_slice(changed_right, right_route, 0, right_cut);
        add_slice(changed_right, left_route, left_cut,
                  current[left_route].size());
      } else {
        add_slice(changed_left, left_route, 0, left_cut);
        add_slice(changed_left, right_route, 0, right_cut, true);
        add_slice(changed_right, left_route, left_cut,
                  current[left_route].size(), true);
        add_slice(changed_right, right_route, right_cut,
                  current[right_route].size());
      }
      consider(two_changes(std::move(changed_left),
                           std::move(changed_right)));
    }
  }

  const bool neighborhood_completed = !expired(deadline);
  std::optional<MoveMemoryEntry> best_move;
  double best_delta = std::numeric_limits<double>::infinity();
  for (std::size_t support = 0; support < support_count; ++support) {
    SupportWorkingMemory& memory = support_working_memory_[support];
    if (memory.epoch != support_epoch_ || !memory.active) continue;
    if (memory.state != SupportState::Closed && neighborhood_completed) {
      move_memory_.put(memory.key, memory.entry);
    }
    double delta = memory.entry.candidate_cost;
    for (std::size_t index = 0; index < memory.key.route_count; ++index) {
      delta -= current_route_costs[memory.key.route_indices[index]];
    }
    if (!std::isfinite(delta) || delta >= -kCostTolerance) {
      continue;
    }
    if (!best_move.has_value() || delta < best_delta - kCostTolerance ||
        (std::abs(delta - best_delta) <= kCostTolerance &&
         move_less(memory.entry.move, best_move->move))) {
      best_move = memory.entry;
      best_delta = delta;
    }
  }
  if (!best_move.has_value()) return std::nullopt;

  const RouteChanges best_changes = materialize_changes(best_move->move);

  struct PreparedChange {
    std::size_t route{};
    std::optional<Plan> plan;
  };
  std::vector<PreparedChange> prepared;
  prepared.reserve(best_changes.size());
  for (const auto& [route, sequence] : best_changes) {
    std::optional<Plan> plan;
    if (!sequence.empty()) {
      plan = prepared_plan(sequence);
      if (completes_vnd_winner(profile_)) {
        ++statistics_.linked_completion_attempts;
        const auto direct_gaps = plan->gaps;
        plan = path_sampler_.repair_new_gaps(
            std::move(*plan), reference_gaps, penalty_lambda, deadline);
        if (plan->gaps != direct_gaps) {
          ++statistics_.linked_completion_changed_attempts;
        }
      }
    }
    prepared.push_back(PreparedChange{
        .route = route,
        .plan = std::move(plan)});
  }
  Solution best;
  best.plans.reserve(current.size() + 1);
  std::size_t change_index = 0;
  for (std::size_t route = 0; route < current.size(); ++route) {
    if (change_index < prepared.size() &&
        prepared[change_index].route == route) {
      if (prepared[change_index].plan.has_value()) {
        best.plans.push_back(*prepared[change_index].plan);
      }
      ++change_index;
    } else {
      best.plans.push_back(solution.plans[route]);
    }
  }
  if (change_index < prepared.size() &&
      prepared[change_index].route == current.size()) {
    if (prepared[change_index].plan.has_value()) {
      best.plans.push_back(*prepared[change_index].plan);
    }
    ++change_index;
  }
  if (change_index != prepared.size()) {
    throw std::logic_error("cached candidate route change is out of range");
  }
  return best;
}

std::optional<Solution> VariableNeighborhoodDescent::replace_path_neighbor(
    const Solution& solution, double penalty_lambda) {
  ++statistics_.replace_path_calls;
  if (solution.plans.empty()) return std::nullopt;
  const double current_cost =
      generalized_cost(solution, instance_, penalty_lambda);
  auto candidate = path_sampler_.random_replacement(solution);
  if (candidate.has_value() &&
      generalized_cost(*candidate, instance_, penalty_lambda) <
          current_cost - kCostTolerance) {
    ++statistics_.replace_path_accepted;
    return candidate;
  }
  // Section 5.1.5 samples one current path and one alternative path.
  return std::nullopt;
}

Solution VariableNeighborhoodDescent::inserted_replace_path(
    Solution solution, double penalty_lambda,
    std::span<const AnchorPair> dirty_anchor_pairs,
    std::optional<std::chrono::steady_clock::time_point> deadline) {
  if (dirty_anchor_pairs.empty()) return solution;
  LimitedRandomPathResult result = path_sampler_.limited_random_replacement(
      std::move(solution), penalty_lambda, dirty_anchor_pairs, deadline);
  statistics_.replace_path_calls += result.attempts;
  statistics_.replace_path_accepted += result.accepted;
  statistics_.accepted_moves += result.accepted;
  statistics_.adaptive_inserted_replace_path_failures +=
      result.failed_attempts;
  return std::move(result.solution);
}

Solution VariableNeighborhoodDescent::final_replace_path(
    Solution solution, double penalty_lambda,
    std::optional<std::chrono::steady_clock::time_point> deadline) {
  ++statistics_.adaptive_final_replace_path_calls;
  std::unordered_set<AnchorPair, AnchorPairHash> eligible;
  for (const Plan& plan : solution.plans) {
    std::vector<int> anchors;
    anchors.reserve(plan.customer_ids.size() + 2);
    anchors.push_back(instance_.depot().id);
    anchors.insert(anchors.end(), plan.customer_ids.begin(),
                   plan.customer_ids.end());
    anchors.push_back(instance_.depot().id);
    for (std::size_t gap = 0; gap < plan.gaps.size(); ++gap) {
      // Destroy/repair changes customer topology only.  The final path phase
      // therefore owns every direct and CS gap in the VND endpoint.
      eligible.insert(AnchorPair{anchors[gap], anchors[gap + 1]});
    }
  }
  const std::vector<AnchorPair> eligible_pairs(eligible.begin(),
                                               eligible.end());
  QuasiExhaustivePathResult result =
      path_sampler_.quasi_exhaustive_replacement(
          std::move(solution), penalty_lambda, eligible_pairs, deadline);
  statistics_.adaptive_final_replace_path_gaps_scanned +=
      result.gaps_scanned;
  statistics_.adaptive_final_replace_path_accepted += result.accepted;
  statistics_.replace_path_calls += result.gaps_scanned;
  statistics_.replace_path_accepted += result.accepted;
  statistics_.accepted_moves += result.accepted;
  return std::move(result.solution);
}

Solution build_solution_from_sequences(
    const std::vector<std::vector<int>>& sequences, PlanFactory& factory,
    const Solution* reference) {
  Solution result;
  result.plans.reserve(sequences.size());
  std::unordered_map<std::vector<int>, const Plan*, IntVectorHash>
      unchanged_plans;
  if (reference != nullptr) {
    unchanged_plans.reserve(reference->plans.size());
    for (const Plan& plan : reference->plans) {
      unchanged_plans.emplace(plan.customer_ids, &plan);
    }
  }
  const auto reference_gaps =
      reference == nullptr ? PlanFactory::GapMap{} : factory.gap_map(*reference);
  for (const auto& sequence : sequences) {
    if (!sequence.empty()) {
      const auto unchanged = unchanged_plans.find(sequence);
      if (unchanged != unchanged_plans.end()) {
        result.plans.push_back(*unchanged->second);
      } else {
        result.plans.push_back(reference == nullptr
                                   ? factory.plan_sequence(sequence)
                                   : factory.plan_sequence(sequence,
                                                           reference_gaps));
      }
    }
  }
  return result;
}

}  // namespace ils_sp
