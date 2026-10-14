#include "ils_sp/search.hpp"
#include "route_opening_rp_reserve.hpp"
#include "bounded_rp_proxy.hpp"

#ifndef ILS_SP_PACKED_DFS_MASK
#define ILS_SP_PACKED_DFS_MASK 0
#endif

#ifndef ILS_SP_TYPED_DELETE_ENTRY
#define ILS_SP_TYPED_DELETE_ENTRY 0
#endif

static_assert(ILS_SP_TYPED_DELETE_ENTRY == 0,
              "Uniform RP candidate requires ILS_SP_TYPED_DELETE_ENTRY=0");

#ifndef ILS_SP_FORWARD_ONLY_TWO_OPT
#define ILS_SP_FORWARD_ONLY_TWO_OPT 0
#endif

#ifndef ILS_SP_PHYSICAL_REVERSE_EMPTY_CHILD
#define ILS_SP_PHYSICAL_REVERSE_EMPTY_CHILD 0
#endif
#include "ils_sp/detail/xxhash_api.hpp"

#include <algorithm>
#include <bit>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fcntl.h>
#include <linux/fs.h>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <sys/vfs.h>
#include <sys/stat.h>
#include <sys/syscall.h>
#include <unistd.h>
#if defined(__x86_64__)
#include <xmmintrin.h>
#endif

namespace ils_sp {
namespace {

#if ILS_SP_PACKED_DFS_MASK
// The caller checks nonzero. Clearing the lowest set bit preserves the
// original ascending station-index DFS order, including shared prefixes.
[[nodiscard]] std::size_t take_lowest_station_index(
    std::uint32_t& remaining) noexcept {
  const auto index = static_cast<std::size_t>(std::countr_zero(remaining));
  remaining &= remaining - 1;
  return index;
}
#endif

// Decision-local physical labels. Customer cut k is immediately after customer
// k-1 (or the departure depot for k=0), leaving the following CS chain in the
// suffix. No charging amounts are inherited; reverse edges are built afresh.
class PhysicalReverseTwoOptProfile {
 public:
  PhysicalReverseTwoOptProfile(const Instance& instance, const Plan& plan,
                              const NonlinearSequenceEvaluator& evaluator)
      : visits_(plan.route.visits), cuts_{0},
        prefix_(visits_.size() + 1),
        reverse_prefix_(ILS_SP_FORWARD_ONLY_TWO_OPT ? 0 : visits_.size() + 1),
        suffix_(visits_.size() + 1),
        reverse_suffix_(ILS_SP_FORWARD_ONLY_TWO_OPT ? 0 : visits_.size() + 1) {
    for (std::size_t p = 0; p < visits_.size(); ++p) {
      if (instance.node(visits_[p]).kind == NodeKind::Customer) cuts_.push_back(p + 1);
    }
    if (cuts_.size() != plan.customer_ids.size() + 1) {
      throw std::logic_error("physical 2OPT customer metadata differs");
    }
    // Ordinary customer cut n remains before its trailing station chain.
    // Only the structural depot-edge descriptor uses n+1: the physical end.
    cuts_.push_back(visits_.size());
    for (std::size_t end = 2; end <= visits_.size(); ++end) {
      prefix_[end] = evaluator.concatenate(prefix_[end - 1],
          evaluator.edge_label(visits_[end - 2], {}, visits_[end - 1]));
#if !ILS_SP_FORWARD_ONLY_TWO_OPT
      reverse_prefix_[end] = evaluator.concatenate(
          evaluator.edge_label(visits_[end - 1], {}, visits_[end - 2]),
          reverse_prefix_[end - 1]);
#endif
    }
    for (std::size_t begin = visits_.size(); begin-- > 0;) {
      if (begin + 1 == visits_.size()) continue;
      suffix_[begin] = evaluator.concatenate(
          evaluator.edge_label(visits_[begin], {}, visits_[begin + 1]),
          suffix_[begin + 1]);
#if !ILS_SP_FORWARD_ONLY_TWO_OPT
      reverse_suffix_[begin] = evaluator.concatenate(reverse_suffix_[begin + 1],
          evaluator.edge_label(visits_[begin + 1], {}, visits_[begin]));
#endif
    }
  }

  void append_label(std::size_t customer_cut, bool prefix, bool reversed,
                    ChargingSequenceLabel& label, int& previous,
                    const NonlinearSequenceEvaluator& evaluator) const {
    const std::size_t cut = cuts_.at(customer_cut);
    const std::size_t begin = prefix ? 0 : cut;
    const std::size_t end = prefix ? cut : visits_.size();
    if (begin == end) return;
    const int first = visits_[reversed ? end - 1 : begin];
    evaluator.append(label, evaluator.edge_label(previous, {}, first));
    evaluator.append(label, prefix
        ? (reversed ? reverse_prefix_[end] : prefix_[end])
        : (reversed ? reverse_suffix_[begin] : suffix_[begin]));
    previous = visits_[reversed ? begin : end - 1];
  }

  void append_visits(std::size_t customer_cut, bool prefix, bool reversed,
                     std::vector<int>& result) const {
    const std::size_t cut = cuts_.at(customer_cut);
    const std::size_t begin = prefix ? 0 : cut;
    const std::size_t end = prefix ? cut : visits_.size();
    if (reversed) {
      for (std::size_t p = end; p > begin; --p) result.push_back(visits_[p - 1]);
    } else {
      result.insert(result.end(), visits_.begin() + static_cast<std::ptrdiff_t>(begin),
                    visits_.begin() + static_cast<std::ptrdiff_t>(end));
    }
  }

 private:
  std::vector<int> visits_;
  std::vector<std::size_t> cuts_;
  std::vector<ChargingSequenceLabel> prefix_, reverse_prefix_, suffix_, reverse_suffix_;
};

constexpr double kLinkedBoundSafety = 100.0 * kCostTolerance;

bool supported_bound_fp_environment() noexcept {
#if defined(__x86_64__) && !defined(__FAST_MATH__) && \
    (!defined(__FINITE_MATH_ONLY__) || __FINITE_MATH_ONLY__ == 0)
  unsigned short control;
  asm volatile("fnstcw %0" : "=m"(control));
  const unsigned csr = _mm_getcsr();
  return (control & 0x003f) == 0x003f && (control & 0x0c00) == 0 &&
         (csr & 0x1f80) == 0x1f80 && (csr & 0xe040) == 0;
#else
  return false;
#endif
}

// This is a certificate for the actual double implementation, not merely for
// the mathematical PWL curve. Unsupported inputs keep the original evaluator.
bool certifies_nonlinear_travel_bound(const Instance& instance) {
  if (!supported_bound_fp_environment()) return false;
  const auto& vehicle = instance.vehicle();
  const double capacity = vehicle.battery_capacity_wh;
  if (capacity < 1e-6 || capacity > 1e6 ||
      vehicle.speed_kmph < 1e-6 || vehicle.consumption_wh_per_km < 1e-6 ||
      vehicle.max_route_duration_h < 1e-6 ||
      vehicle.max_route_duration_h > 1e6) return false;
  for (const auto& node : instance.nodes()) {
    if (node.service_time_h > 1e6) return false;
    for (const auto& target : instance.nodes()) {
      const double energy = instance.energy_wh(node.id, target.id);
      if (!std::isfinite(energy) || energy < 0.0 || energy > 1e6) return false;
      const double travel = instance.travel_time_h(node.id, target.id);
      if (!std::isfinite(travel) || travel < 0.0 || travel > 1e9) return false;
    }
  }
  for (const int id : instance.station_ids()) {
    const auto& curve = instance.curve_for_station(id);
    const auto& points = curve.points();
    if (curve.capacity_wh() != capacity || points.front().time_h != 0.0 ||
        points.front().energy_wh != 0.0) return false;
    for (std::size_t index = 0; index < points.size(); ++index) {
      const auto& point = points[index];
      if (point.time_h > 1e6 ||
          !std::isfinite(curve.time_at_energy(point.energy_wh))) return false;
      if (index != 0) {
        const double energy_width = point.energy_wh - points[index - 1].energy_wh;
        const double time_width = point.time_h - points[index - 1].time_h;
        if ((energy_width > 0.0 && energy_width < 1e-6) ||
            (time_width > 0.0 && time_width < 1e-12)) return false;
      }
      // Zero and the clamped capacity have no downward join. Within a segment
      // all positive-width interpolation operations are monotone. Test every
      // other join using the original function and the next representable E.
      if (point.energy_wh > 0.0 && point.energy_wh < capacity &&
          curve.time_at_energy(point.energy_wh) > curve.time_at_energy(
              std::nextafter(point.energy_wh,
                             std::numeric_limits<double>::infinity()))) {
        return false;
      }
    }
  }
  return true;
}

std::optional<double> nonlinear_travel_bound(
    double travel, double service, std::size_t station_count,
    const Vehicle& vehicle, double lambda) {
  if (station_count > 62 ||
      !std::isfinite(travel) || travel < 0.0 || travel > 1e12 ||
      !std::isfinite(service) || service < 0.0 || service > 1e12 ||
      !std::isfinite(lambda) ||
      (lambda != 0.0 && (lambda < 1e-6 || lambda > 1e4))) return std::nullopt;
  const double no_charge_duration = travel + service;
  const double excess =
      std::max(0.0, no_charge_duration - vehicle.max_route_duration_h);
  if (lambda == 0.0 || excess == 0.0) return travel;
  const double downward = -std::numeric_limits<double>::infinity();
  const double penalty = std::max(0.0, std::nextafter(lambda * excess, downward));
  if (penalty == 0.0) return travel;
  return std::max(0.0, std::nextafter(travel + penalty, downward));
}

std::optional<double> linked_triangle_slack(const Instance& instance,
                                            bool data_certified) {
  if (!data_certified || instance.station_ids().size() > 31) return std::nullopt;
  const double down = -std::numeric_limits<double>::infinity();
  const double up = std::numeric_limits<double>::infinity();
  double slack = 0.0;
  for (const auto& a : instance.nodes()) {
    for (const auto& b : instance.nodes()) {
      const double direct = instance.travel_time_h(a.id, b.id);
      for (const int station : instance.station_ids()) {
        const double via = std::nextafter(
            instance.travel_time_h(a.id, station) +
                instance.travel_time_h(station, b.id), down);
        slack = std::max(slack, std::nextafter(direct - via, up));
      }
    }
  }
  return slack;
}

// Literal customer-only opening: the existing route retains every station
// occurrence in parent order; the new empty route receives only the customer.
// Charging and infeasibility are evaluated afresh by the unchanged evaluator.
std::vector<int> project_opening_visits(const Instance& instance,
                                      const Plan& parent, int moved_customer,
                                      bool singleton_child) {
  if (singleton_child) return {moved_customer};
  std::vector<int> visits;
  visits.reserve(parent.route.visits.size());
  for (int id : parent.route.visits) {
    if (instance.node(id).kind == NodeKind::Station || id != moved_customer) {
      visits.push_back(id);
    }
  }
  return visits;
}

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

// No energy/rate conservation claim is valid for tolerance-clamped Algorithm 3.
// The per-instance certificate below bounds only travel and unavoidable horizon.
[[nodiscard]] double linked_route_lower_bound(
    const RouteCostSummary& summary, const Instance&,
    double penalty_lambda, const PathSampler* sampler,
    std::size_t source_visits = 0, std::size_t source_customers = 0) {
  return sampler ? sampler->linked_insertion_lower_bound(
                       summary, penalty_lambda, source_visits, source_customers)
                 : -std::numeric_limits<double>::infinity();
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
    bool calculate_linked_lower_bound, const PathSampler* linked_bound_sampler,
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
                                     linked_bound_sampler, old_plan.route.visits.size(),
                                     old_plan.customer_ids.size()) -
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
        linked_bound_sampler);
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
    double penalty_lambda, const PathSampler* linked_bound_sampler) {
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
                               linked_bound_sampler, old_plan.route.visits.size(),
                               old_plan.customer_ids.size()) -
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
    const PathSampler* linked_bound_sampler, std::mt19937_64& random,
    InsertionCandidateStatistics& statistics) {
  statistics.granular_escape_candidates_seen += rejected.size();
  if (rejected.empty()) return;
  const RejectedInsertionDescriptor& selected =
      rejected[random_index(rejected, random)];
  candidates.push_back(rejected_insertion_move(
      current, selected, instance, scorer, penalty_lambda,
      linked_bound_sampler));
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
    bool repair_affected_paths = false,
    bool single_step_completion = false) {
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
          std::move(completed), reference_gaps, penalty_lambda, std::nullopt,
          single_step_completion);
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
    InsertionCompletionMemo& completion_memo,
    bool single_step_completion = false) {
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
                            statistics, completion_memo, false,
                            single_step_completion);
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
                              statistics, completion_memo, false,
                              single_step_completion);
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
    bool require_complete_ranking, bool single_step_completion) {
  statistics.candidates_generated += moves.size();
  const std::size_t direct_best_index = *direct.indices.best;
  InsertionMove direct_best = direct.best;

  BestTwoInsertionIndices exact;
  const auto evaluate_and_rank = [&](std::size_t index) {
    evaluate_linked_insertion(moves[index], current, instance, factory,
                              path_sampler, reference_gaps, penalty_lambda,
                              statistics, completion_memo, false,
                              single_step_completion);
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
         ((candidate_regret == incumbent_regret ||
           std::abs(candidate_regret - incumbent_regret) <= kCostTolerance) &&
          insertion_better(candidate, incumbent));
}

[[nodiscard]] InsertionMove choose_linked_two_regret(
    const Solution& current, const std::vector<int>& unassigned,
    const Instance& instance, PlanFactory& factory, PathSampler& path_sampler,
    const RouteSequenceScorer& scorer,
    const PlanFactory::GapMap& reference_gaps, double penalty_lambda,
    const PathSampler* linked_bound_sampler,
    InsertionCandidateStatistics& statistics,
    InsertionCompletionMemo& completion_memo,
    bool single_step_completion = false) {
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
        linked_bound_sampler);
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
        index == direct_winner_index, single_step_completion);
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

// A conservative lower bound, not an exact label travel value. The catalog's
// distance includes prefix-first subtraction and reconstruction. For <=63
// legs, 512 epsilon covers both that rounding and per-leg division/summing;
// 512 denorm_min covers gradual underflow. See the audit's suffix proof.
double PathSampler::catalog_travel_lower_bound(double distance, double speed) {
  const double downward = -std::numeric_limits<double>::infinity();
  const double factor = 1.0 - 512.0 * std::numeric_limits<double>::epsilon();
  const double scaled = std::nextafter((distance / speed) * factor, downward);
  return std::max(0.0, std::nextafter(
      scaled - 512.0 * std::numeric_limits<double>::denorm_min(), downward));
}

double PathSampler::linked_insertion_lower_bound(
    const RouteCostSummary& summary, double lambda,
    std::size_t source_visits, std::size_t source_customers) const {
  const double disabled = -std::numeric_limits<double>::infinity();
  const std::size_t stations = instance_.station_ids().size();
  // New insertions create exactly two direct gaps. Each catalogue path is
  // simple, hence can add at most m stations per gap. Inherited gaps are fixed.
  // Bound both the LabelTree and completed sequential sums by physical legs.
  if (!complete_paths_ || implicit_replacement_ || !linked_triangle_slack_h_ ||
      !supported_bound_fp_environment() ||
      source_customers > source_visits || source_visits > 317 || stations > 31)
    return disabled;
  const std::size_t inserted_stations = 2 * stations;
  const std::size_t legs = source_visits + 2 + inserted_stations;
  const std::size_t total_stations =
      source_visits - source_customers + inserted_stations;
  if (legs > 319 || total_stations > 62 ||
      !nonlinear_travel_bound(summary.travel_time_h, summary.service_time_h,
                              total_stations, instance_.vehicle(), lambda))
    return disabled;
  const double down = -std::numeric_limits<double>::infinity();
  const double up = std::numeric_limits<double>::infinity();
  // For nonnegative sums with <=L additions, gamma_L covers any addition tree.
  // (1-gamma_L)/(1+gamma_L) >= 1-4*L*epsilon in this L<=319 envelope.
  // The explicit denormal allowance also covers gradual-underflow additions.
  const double factor = std::nextafter(
      1.0 - (4.0 * static_cast<double>(legs)) *
                  std::numeric_limits<double>::epsilon(), down);
  const double underflow = std::nextafter(
      (4.0 * static_cast<double>(legs)) *
          std::numeric_limits<double>::denorm_min(), up);
  const double detour_slack = std::nextafter(
      static_cast<double>(inserted_stations) * *linked_triangle_slack_h_, up);
  const auto lower_sum = [&](double value, double removed_slack) {
    return std::max(0.0, std::nextafter(
        std::nextafter(std::nextafter(value * factor, down) - removed_slack,
                       down) - underflow, down));
  };
  const double travel = lower_sum(summary.travel_time_h, detour_slack);
  const double service = lower_sum(summary.service_time_h, 0.0);
  // Certified actual charging is nonnegative. The existing directed-rounding
  // kernel bounds travel + lambda*max(0,travel+service-Tmax); energy is omitted.
  return nonlinear_travel_bound(travel, service, total_stations,
                                 instance_.vehicle(), lambda).value_or(disabled);
}

std::optional<double> PathSampler::zero_shift_energy_cost_bound(
    double travel, double service, const Vehicle& vehicle, double lambda) {
  // Caller certifies original finite-rate zero-shift propagation, canonical
  // <=319 physical legs, <=62 CS, and the existing data/FP envelope. Restrict
  // both conversion factors as well to bound gradual-underflow amplification.
  if (vehicle.speed_kmph > 1e6 || vehicle.consumption_wh_per_km > 1e6 ||
      !nonlinear_travel_bound(travel, service, 62, vehicle, lambda))
    return std::nullopt;
  const double no_charge_duration = travel + service;
  if (no_charge_duration < vehicle.max_route_duration_h) return std::nullopt;
  const double energy = (travel * vehicle.speed_kmph) * vehicle.consumption_wh_per_km;
  if (!std::isfinite(energy) || energy > 4e8) return std::nullopt;
  // In this scope there is no actual charging. Original tolerance and all
  // travel/energy summation, conversion, SOC subtraction and deficit summation
  // errors total <2e-4 Wh (audit). Never use exact conservation or charge/rate.
  constexpr double kEnergyRoundoffAllowanceWh = 1e-3;
  const double downward = -std::numeric_limits<double>::infinity();
  const double deficit = std::nextafter(
      std::nextafter(energy - vehicle.battery_capacity_wh, downward) -
          kEnergyRoundoffAllowanceWh, downward);
  const double energy_penalty =
      std::max(0.0, deficit) / vehicle.consumption_wh_per_km;
  return travel + lambda * (energy_penalty +
      std::max(0.0, no_charge_duration - vehicle.max_route_duration_h));
}

GranularNeighborhood::GranularNeighborhood(const Instance& instance,
                                           std::size_t maximum_neighbors) {
  if (maximum_neighbors == 0) {
    throw std::invalid_argument("promising arc count must be positive");
  }
  const std::size_t count = instance.customer_ids().size();
  customer_rank_by_id_.reserve(count);
  consecutive_customer_ids_ = true;
  for (std::size_t rank = 0; rank < count; ++rank) {
    const int customer = instance.customer_ids()[rank];
    if (customer <= 0 || static_cast<std::size_t>(customer) != rank + 1) {
      consecutive_customer_ids_ = false;
    }
    customer_rank_by_id_.emplace(customer, rank);
  }
  promising_matrix_.assign(count * count, 0);
  if (count <= 1) {
    gamma_size_ = 0;
    for (const int customer : instance.customer_ids()) {
      directed_neighbors_.try_emplace(customer);
      promising_neighbors_.try_emplace(customer);
    }
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
  if (consecutive_customer_ids_ && customer_id > 0 &&
      static_cast<std::size_t>(customer_id) <= customer_rank_by_id_.size()) {
    return static_cast<std::size_t>(customer_id - 1);
  }
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
  if (consecutive_customer_ids_) {
    const std::size_t count = customer_rank_by_id_.size();
    if (left_id <= 0 || right_id <= 0 ||
        static_cast<std::size_t>(left_id) > count ||
        static_cast<std::size_t>(right_id) > count) return false;
    return promising_matrix_[(static_cast<std::size_t>(left_id) - 1) * count +
                             static_cast<std::size_t>(right_id) - 1] != 0;
  }
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
  // Generic new-arc query, not an operator-specific VND predicate. Each VND
  // customer operator supplies its own entry/reconnection condition before
  // scoring; empty-segment/depot mappings are documented in audit section 5.15.
  // Depot-only changes have no customer arc to satisfy this guard.
  return false;
}

std::uint64_t PathSampler::PathCatalog::encode_path(
    std::span<const int> stations, std::uint32_t flat_offset) const {
  if (!packed) {
    return static_cast<std::uint64_t>(flat_offset) |
        (static_cast<std::uint64_t>(stations.size()) << 32);
  }
  if (stations.size() > packed_station_ids.size()) {
    throw std::length_error("packed station path is too long");
  }
  const auto known =
      std::span<const int>{packed_station_ids}.first(packed_station_count);
  std::uint64_t encoded = 0;
  for (std::size_t index = 0; index < stations.size(); ++index) {
    const auto found =
        std::lower_bound(known.begin(), known.end(), stations[index]);
    if (found == known.end() || *found != stations[index]) {
      throw std::invalid_argument("packed path contains an unknown station");
    }
    const auto code = static_cast<std::uint64_t>(found - known.begin() + 1);
    encoded |= code << (60 - 4 * index);
  }
  return encoded;
}

std::size_t PathSampler::PathCatalog::path_size(
    const PathDescriptor& path) const noexcept {
  if (!packed) return static_cast<std::uint32_t>(path.station_storage >> 32);
  // All used nibbles are nonzero and the lowest nibble is always a sentinel.
  // countr_zero(0) is 64, so the empty path needs no special representation.
  return 16 - static_cast<std::size_t>(std::countr_zero(path.station_storage)) / 4;
}

std::span<const int> PathSampler::PathCatalog::stations(
    const PathDescriptor& path, DecodeBuffer& decoded) const {
  const std::size_t size = path_size(path);
  if (!packed) {
    return std::span<const int>{station_ids}.subspan(
        static_cast<std::uint32_t>(path.station_storage), size);
  }
  for (std::size_t index = 0; index < size; ++index) {
    const auto code = static_cast<std::size_t>(
        (path.station_storage >> (60 - 4 * index)) & 15);
    decoded[index] = packed_station_ids[code - 1];
  }
  return std::span<const int>{decoded}.first(size);
}

bool PathSampler::PathCatalog::lexically_less(
    const PathDescriptor& left, const PathDescriptor& right) const {
  // Sorted station ranks and zero trailing nibbles preserve the original
  // integer-ID lexicographic order, including a shorter shared prefix.
  if (packed) return left.station_storage < right.station_storage;
  DecodeBuffer left_buffer, right_buffer;
  const auto first = stations(left, left_buffer);
  const auto second = stations(right, right_buffer);
  return std::lexicographical_compare(first.begin(), first.end(),
                                      second.begin(), second.end());
}

std::size_t PathSampler::PathCatalog::retained_bytes() const noexcept {
  if (disk_backing) {
    const auto& header = disk_backing->header;
    return sizeof(PathCatalogPayload) + header[11] * sizeof(PathDescriptor) +
        header[12] * sizeof(int) + header[13] * sizeof(std::uint32_t) + 256;
  }
  return sizeof(PathCatalogPayload) + paths.capacity() * sizeof(PathDescriptor) +
      station_ids.capacity() * sizeof(int) +
      lexical_order.capacity() * sizeof(std::uint32_t) + 256;
}

std::size_t PathSampler::PathCatalog::size() const noexcept {
  return disk_backing ? static_cast<std::size_t>(disk_backing->header[8]) : paths.size();
}

#include "catalog_disk_cache.inc"
#include "bounded_rp_proxy.inc"

PathSampler::PathSampler(const Instance& instance, PlanFactory& factory,
                         std::mt19937_64& random, bool complete_paths,
                         bool implicit_replacement)
    : instance_(instance), factory_(factory), random_(random),
      complete_paths_(complete_paths), implicit_replacement_(implicit_replacement),
      maximum_charging_rate_wh_per_h_(maximum_charging_rate(instance)) {
  build_curve_order();
  build_reachable_station_index();
  if (!complete_paths_) build_shortest_station_paths();
  replacement_summary_statistics_.bound_certified =
      certifies_nonlinear_travel_bound(instance_);
  linked_triangle_slack_h_ = linked_triangle_slack(
      instance_, complete_paths_ && !implicit_replacement_ &&
                     replacement_summary_statistics_.bound_certified);
  initialize_catalog_disk_cache();
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
  // Pairwise epsilon equivalence is not transitive and cannot be used by
  // std::sort. Compare exact finite signatures on one common energy grid.
  // Tolerance remains only in the separate dominance validation below.
  std::vector<double> energy_breakpoints;
  for (const auto& type : curve_types) {
    for (const ChargingPoint& point : instance_.charging_curves().at(type).points()) {
      energy_breakpoints.push_back(point.energy_wh);
    }
  }
  std::sort(energy_breakpoints.begin(), energy_breakpoints.end());
  energy_breakpoints.erase(
      std::unique(energy_breakpoints.begin(), energy_breakpoints.end()),
      energy_breakpoints.end());
  std::map<std::string, std::vector<double>> curve_signatures;
  for (const auto& type : curve_types) {
    const ChargingCurve& curve = instance_.charging_curves().at(type);
    auto& signature = curve_signatures[type];
    for (std::size_t index = 1; index < energy_breakpoints.size(); ++index) {
      const double duration = curve.charge_duration(
          energy_breakpoints[index - 1], energy_breakpoints[index]);
      if (!std::isfinite(duration)) {
        throw std::invalid_argument("charging curve ordering requires finite durations");
      }
      signature.push_back(duration);
    }
  }
  std::sort(curve_types.begin(), curve_types.end(),
            [&](const std::string& left, const std::string& right) {
              const auto& first = curve_signatures.at(left);
              const auto& second = curve_signatures.at(right);
              if (first != second) return first < second;
              return left < right;
            });
  // A rank is consumed as a total dominance chain. Validate every pair:
  // adjacent epsilon comparisons alone would not establish transitivity.
  for (std::size_t left = 0; left < curve_types.size(); ++left) {
    for (std::size_t right = left + 1; right < curve_types.size(); ++right) {
      if (!charging_curve_no_slower(
              instance_.charging_curves().at(curve_types[left]),
              instance_.charging_curves().at(curve_types[right]))) {
        throw std::invalid_argument(
            "paper path dominance requires a consistent charging curve order");
      }
    }
  }
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

std::optional<std::vector<int>> PathSampler::sample_explicit_path(
    int origin_id, int target_id, const std::vector<int>& current_path) {
  const auto catalog = path_catalog(origin_id, target_id);
  CatalogReader reader(*catalog);
  PathCatalog::DecodeBuffer decoded;
  const auto [first, last] = reader.equal_range(current_path);
  // Same sorted multiset, all copies excluded, same one uniform-size_t draw.
  const std::size_t excluded = last - first;
  const std::size_t available = catalog->size() - excluded;
  if (available == 0) return std::nullopt;
  // A lazy read can fail after the ordinal is selected. Commit the engine
  // only after the selected physical path is fully read and materialized.
  // Resident catalogues retain the original direct-RNG fast path.
  std::optional<std::mt19937_64> pending_random;
  if (catalog->disk_backing) pending_random.emplace(random_);
  auto& engine = pending_random ? *pending_random : random_;
  std::size_t ordinal =
      std::uniform_int_distribution<std::size_t>(0, available - 1)(engine);
  if (ordinal >= first) ordinal += excluded;
  const auto stations = catalog->stations(
      reader.descriptor(reader.lexical_index(ordinal)), decoded);
  std::vector<int> result(stations.begin(), stations.end());
  if (pending_random) random_ = *pending_random;
  return result;
}

std::shared_ptr<const PathSampler::PathCatalog> PathSampler::path_catalog(
    int origin_id, int target_id, RepairTimingStatistics* timing,
    PathCatalogStatistics* statistics) {
  if (statistics != nullptr) *statistics = {};
  const AnchorPair key{origin_id, target_id};
  if (auto cached = catalog_cache_.get(key); cached.has_value()) {
    if (timing != nullptr) ++timing->path_option_cache_hits;
    if (statistics != nullptr) statistics->cache_hit = true;
    return *cached;
  }
  if (timing != nullptr) ++timing->path_option_cache_misses;
  const auto disk_started = catalog_disk_cache_.enabled &&
      (timing != nullptr || statistics != nullptr)
      ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
  if (auto restored = disk_cache_load(key)) {
    // Same object lifetime and RAM admission as an ordinary rebuild. In
    // particular, old weak_ptr physical-summary rows are not resurrected.
    const auto evictions_before = catalog_cache_.evictions();
    const std::size_t bytes = restored->retained_bytes();
    if (bytes <= catalog_cache_budget_) {
      while (catalog_cache_.size() != 0 &&
             catalog_cache_bytes_ + bytes > catalog_cache_budget_) {
        const auto evicted = catalog_cache_.discard_oldest();
        catalog_cache_bytes_ -= (*evicted)->retained_bytes();
      }
      catalog_cache_.put(key, restored);
      catalog_cache_bytes_ += bytes;
    }
    if (timing != nullptr) {
      timing->path_option_cache_evictions +=
          catalog_cache_.evictions() - evictions_before;
      timing->path_option_cache_miss_build_s += elapsed_seconds(disk_started);
    }
    // Optional catalogue setup timing includes a backing-store restoration.
    if (statistics != nullptr) statistics->setup_s += elapsed_seconds(disk_started);
    return restored;
  }
  if (catalog_disk_cache_.enabled && (timing != nullptr || statistics != nullptr)) {
    const double seconds = elapsed_seconds(disk_started);
    if (timing != nullptr) timing->path_option_cache_miss_build_s += seconds;
    if (statistics != nullptr) statistics->setup_s += seconds;
  }
  // Allocation failures used to report only bad_alloc.  Trace the failing
  // catalog phase without temporary path/string buffers or success logs.
  struct PathBuildTrace {
    int origin;
    int target;
    int exceptions = std::uncaught_exceptions();
    const char* phase = "setup";
    std::uint64_t generated = 0;
    std::size_t retained = 0;
    PathCatalogStatistics* statistics;
    double PathCatalogStatistics::* phase_seconds =
        &PathCatalogStatistics::setup_s;
    std::chrono::steady_clock::time_point phase_started;
    PathBuildTrace(int origin_id, int target_id,
                   PathCatalogStatistics* profile)
        : origin(origin_id), target(target_id), statistics(profile),
          phase_started(profile != nullptr
              ? std::chrono::steady_clock::now()
              : std::chrono::steady_clock::time_point{}) {}
    void finish_phase() noexcept {
      if (statistics != nullptr) {
        const auto now = std::chrono::steady_clock::now();
        (statistics->*phase_seconds) +=
            std::chrono::duration<double>(now - phase_started).count();
        phase_started = now;
      }
    }
    void enter(const char* next,
               double PathCatalogStatistics::* seconds) noexcept {
      finish_phase();
      phase = next;
      phase_seconds = seconds;
    }
    ~PathBuildTrace() noexcept {
      finish_phase();
      if (std::uncaught_exceptions() > exceptions) {
        std::fprintf(stderr,
            "path_catalog_build_failure origin=%d target=%d phase=%s "
            "generated=%llu retained=%zu\n",
            origin, target, phase,
            static_cast<unsigned long long>(generated), retained);
      }
    }
  } trace{origin_id, target_id, statistics};
  if (option_keys_built_.insert(key).second) {
    ++option_first_builds_;
  } else {
    ++option_rebuilds_;
  }
  const auto build_started = std::chrono::steady_clock::now();
  // Only the construction sequence representation changes. The live DFS path
  // and every prefix-distance operation remain the original ones. A single
  // build uses one pointer alternative, so CandidatePath remains 32 bytes.
  struct PackedStationPath {
    double distance{};
    std::uint64_t stations{};
    std::size_t fastest_curve_rank{};
    std::size_t slowest_curve_rank{};
  };
  const bool packed_construction = complete_paths_ &&
      packed_construction_enabled_ && packed_catalogs_enabled_ &&
      instance_.station_ids().size() <= PathCatalog::DecodeBuffer{}.size();
  struct CandidatePath {
    union {
      const StationPath* station_path{};
      const PackedStationPath* packed_station_path;
    };
    double first_leg{};
    double last_leg{};
    double total_distance{};
  };
  static_assert(sizeof(CandidatePath) == 32);
  const auto fastest_rank = [&](const CandidatePath& candidate) {
    return packed_construction ? candidate.packed_station_path->fastest_curve_rank
                               : candidate.station_path->fastest_curve_rank;
  };
  const auto slowest_rank = [&](const CandidatePath& candidate) {
    return packed_construction ? candidate.packed_station_path->slowest_curve_rank
                               : candidate.station_path->slowest_curve_rank;
  };
  const auto sequence_less = [&](const CandidatePath& left,
                                 const CandidatePath& right) {
    return packed_construction
        ? left.packed_station_path->stations < right.packed_station_path->stations
        : left.station_path->stations < right.station_path->stations;
  };
  const auto dominates = [&](const CandidatePath& left,
                            const CandidatePath& right) {
    return slowest_rank(left) <= fastest_rank(right) &&
           left.first_leg <= right.first_leg + kCostTolerance &&
           left.last_leg <= right.last_leg + kCostTolerance &&
           left.total_distance <= right.total_distance + kCostTolerance;
  };
  std::vector<CandidatePath> candidates;
  // Owners may move during enumeration: no CandidatePath pointing into a
  // bucket escapes a predicate until the whole DFS has finished. The later
  // candidate views are created only after all owner mutations have stopped.
  struct StoredStationPath {
    using Path = StationPath;
    StationPath station_path;
    double first_leg{};
    double last_leg{};
  };
  struct PackedStoredStationPath {
    using Path = PackedStationPath;
    PackedStationPath station_path;
    double first_leg{};
    double last_leg{};
  };
  static_assert(sizeof(PackedStoredStationPath) == 48);
  static_assert(std::is_nothrow_move_constructible_v<StoredStationPath> &&
                std::is_nothrow_move_assignable_v<StoredStationPath> &&
                std::is_nothrow_move_constructible_v<PackedStoredStationPath> &&
                std::is_nothrow_move_assignable_v<PackedStoredStationPath>);
  const auto described = []<typename Stored>(const Stored& stored) {
    // Preserve first + (prefix - first) + last exactly; no shared middle sum.
    CandidatePath result;
    if constexpr (std::is_same_v<Stored, PackedStoredStationPath>) {
      result.packed_station_path = &stored.station_path;
    } else {
      result.station_path = &stored.station_path;
    }
    result.first_leg = stored.first_leg;
    result.last_leg = stored.last_leg;
    result.total_distance = stored.first_leg + stored.station_path.distance +
                            stored.last_leg;
    return result;
  };
  const std::size_t curve_types = curve_rank_by_type_.size();
  const auto& stations = instance_.station_ids();
  // Within a first-station group the original stored first_leg is identical.
  // Check that mandatory dominance predicate once before traversing its paths;
  // retain bucket order, distance arithmetic and the original dominance test.
  const std::size_t first_groups = stations.size();
  // Contiguous owners avoid per-path list nodes and pointer chasing in the
  // hot dominance scans. erase_if keeps survivor order; StationPath's vector
  // moves transfer its buffer rather than copying stations during compaction.
  const std::size_t bucket_count = curve_types * curve_types * first_groups;
  std::vector<std::vector<StoredStationPath>> complete_storage(
      !packed_construction ? bucket_count : 0);
  std::vector<std::vector<PackedStoredStationPath>> packed_complete_storage(
      packed_construction ? bucket_count : 0);
  std::uint64_t& generated_count = trace.generated;
  const auto& reachable_first =
      stations_reachable_from_anchor_.at(origin_id);
  const auto& reachable_last = stations_reaching_anchor_.at(target_id);
  candidates.reserve(reachable_first.size() * reachable_last.size());
  if (complete_paths_) {
    trace.enter("elementary_enumeration", &PathCatalogStatistics::enumeration_s);
    std::vector<bool> used(stations.size(), false);
    std::vector<int> path;
    std::vector<double> first_distances;
    first_distances.reserve(first_groups);
    for (const int station : stations) {
      first_distances.push_back(instance_.distance_km(origin_id, station));
    }
    const double maximum_rate = maximum_charging_rate_wh_per_h_;
    const double capacity = instance_.vehicle().battery_capacity_wh;
    const double consumption = instance_.vehicle().consumption_wh_per_km;
    const double horizon = instance_.vehicle().max_route_duration_h;
    const double minimum_rest_distance =
        instance_.distance_km(instance_.depot().id, origin_id) +
        instance_.distance_km(target_id, instance_.depot().id);
    const double minimum_service = instance_.node(origin_id).service_time_h +
        (origin_id == target_id ? 0.0 : instance_.node(target_id).service_time_h);
    const auto time_lower_bound = [&](double distance) {
      // Any complete depot route containing this path also has to reach its
      // origin and return from its target.  Triangle inequality lower-bounds
      // those portions even when other customers/stations occur there.
      distance += minimum_rest_distance;
      return minimum_service + distance / instance_.vehicle().speed_kmph +
             std::max(0.0, distance * consumption - capacity) / maximum_rate;
    };
    const auto build_complete = [&]<typename Stored>(
        std::vector<std::vector<Stored>>& storage) {
      constexpr bool packed = std::is_same_v<Stored, PackedStoredStationPath>;
      const auto complete_bucket = [&](std::size_t fast, std::size_t slow,
                                       std::size_t first) -> auto& {
        return storage[(fast * curve_types + slow) * first_groups + first];
      };
      const auto enumerate = [&](auto&& self, int previous,
                                  double prefix_distance, double first_leg,
                                  std::size_t fastest, std::size_t slowest,
                                  std::size_t first_index,
                                  std::uint64_t encoded_prefix
#if ILS_SP_PACKED_DFS_MASK
                                  , std::uint32_t remaining_stations
#endif
                                  ) -> void {
#if ILS_SP_PACKED_DFS_MASK
        std::uint32_t choices = remaining_stations;
#endif
        for (std::size_t cursor = 0; cursor < stations.size();) {
          std::size_t index;
#if ILS_SP_PACKED_DFS_MASK
          if constexpr (packed) {
            if (choices == 0) break;
            index = take_lowest_station_index(choices);
          } else
#endif
          {
            index = cursor++;
            if (used[index]) continue;
          }
          const int station = stations[index];
          const double leg = instance_.distance_km(previous, station);
          if (leg * consumption > capacity + kEnergyToleranceWh) continue;
          const double distance = prefix_distance + leg;
          const double last_leg = instance_.distance_km(station, target_id);
          // Triangle inequality and globally fastest charging give a lower
          // bound for every continuation.  No feasible route of length Tmax
          // can contain a path failing this bound; no hop/Top-K cap is used.
          if (time_lower_bound(distance + last_leg) > horizon + kTimeToleranceH)
            continue;
          const double first = path.empty() ? leg : first_leg;
          const std::size_t next_first_index = path.empty() ? index : first_index;
          const std::size_t rank = station_curve_ranks_[index];
          const std::size_t next_fastest = path.empty() ? rank : std::min(fastest, rank);
          const std::size_t next_slowest = path.empty() ? rank : std::max(slowest, rank);
          used[index] = true;
          path.push_back(station);
          std::uint64_t encoded_path = 0;
          if constexpr (packed) {
            // sorted station IDs map to nonzero nibbles, most significant first;
            // depth <= alphabet <= 15 leaves the low nibble as a zero sentinel.
            encoded_path = encoded_prefix |
                (static_cast<std::uint64_t>(index + 1) << (64 - 4 * path.size()));
          }
          if (last_leg * consumption <= capacity + kEnergyToleranceWh) {
            typename Stored::Path proposed;
            proposed.distance = distance - first;
            if constexpr (packed) proposed.stations = encoded_path;
            else proposed.stations = path;
            proposed.fastest_curve_rank = next_fastest;
            proposed.slowest_curve_rank = next_slowest;
            CandidatePath proposal;
            if constexpr (packed) proposal.packed_station_path = &proposed;
            else proposal.station_path = &proposed;
            proposal.first_leg = first;
            proposal.last_leg = last_leg;
            proposal.total_distance = distance + last_leg;
            ++generated_count;
            bool rejected = false;
            for (std::size_t slow = 0; slow <= next_fastest && !rejected; ++slow) {
              for (std::size_t fast = 0; fast <= slow && !rejected; ++fast) {
                for (std::size_t group_first = 0;
                     group_first < first_groups && !rejected; ++group_first) {
                  if (!(first_distances[group_first] <=
                        proposal.first_leg + kCostTolerance)) continue;
                  for (const Stored& stored :
                       complete_bucket(fast, slow, group_first)) {
                    const CandidatePath known = described(stored);
                    if (dominates(known, proposal) &&
                        (!dominates(proposal, known) ||
                         sequence_less(known, proposal))) {
                      rejected = true;
                      break;
                    }
                  }
                }
              }
            }
            if (!rejected) {
              for (std::size_t fast = next_slowest; fast < curve_types; ++fast) {
                for (std::size_t slow = fast; slow < curve_types; ++slow) {
                  for (std::size_t group_first = 0; group_first < first_groups;
                       ++group_first) {
                    if (!(proposal.first_leg <=
                          first_distances[group_first] + kCostTolerance)) continue;
                    auto& bucket = complete_bucket(fast, slow, group_first);
                    trace.retained -= std::erase_if(bucket, [&](const Stored& stored) {
                      return dominates(proposal, described(stored));
                    });
                  }
                }
              }
              const double stored_first = instance_.distance_km(
                  origin_id, path.front());
              const double stored_last = instance_.distance_km(
                  path.back(), target_id);
              complete_bucket(next_fastest, next_slowest, next_first_index).push_back(
                  Stored{std::move(proposed), stored_first, stored_last});
              ++trace.retained;
            }
          }
          self(self, station, distance, first, next_fastest, next_slowest,
               next_first_index, encoded_path
#if ILS_SP_PACKED_DFS_MASK
               , [&] {
                   if constexpr (packed) {
                     return remaining_stations & ~(std::uint32_t{1} << index);
                   } else {
                     return std::uint32_t{0};
                   }
                 }()
#endif
               );
          path.pop_back();
          used[index] = false;
        }
      };
      enumerate(enumerate, origin_id, 0.0, 0.0, 0, 0, 0, 0
#if ILS_SP_PACKED_DFS_MASK
                , [&] {
                    if constexpr (packed) {
                      // Packed construction is selected only for <=15 stations.
                      return (std::uint32_t{1} << stations.size()) - 1;
                    } else {
                      return std::uint32_t{0};
                    }
                  }()
#endif
                );
      trace.enter("candidate_views", &PathCatalogStatistics::candidate_views_s);
      for (const auto& bucket : storage) {
        for (const Stored& retained : bucket) {
          candidates.push_back(described(retained));
        }
      }
    };
    if (packed_construction) build_complete(packed_complete_storage);
    else build_complete(complete_storage);
  } else {
    trace.enter("compressed_projection", &PathCatalogStatistics::enumeration_s);
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
  }
  if (!complete_paths_) generated_count = candidates.size();
  option_candidates_generated_ += generated_count;
  if (timing != nullptr) {
    timing->path_option_candidates_generated += generated_count;
  }

  trace.enter("candidate_sort", &PathCatalogStatistics::candidate_sort_s);
  std::sort(candidates.begin(), candidates.end(), [&](const CandidatePath& left,
                                                       const CandidatePath& right) {
    const auto left_key = std::tuple{left.first_leg, left.last_leg,
        left.total_distance, fastest_rank(left), slowest_rank(left)};
    const auto right_key = std::tuple{right.first_leg, right.last_leg,
        right.total_distance, fastest_rank(right), slowest_rank(right)};
    if (left_key != right_key) return left_key < right_key;
    return sequence_less(left, right);
  });
  struct RankedCandidate {
    CandidatePath candidate;
    bool active{true};
  };
  const std::size_t curve_count = curve_rank_by_type_.size();
  trace.enter("dominance_filter", &PathCatalogStatistics::dominance_filter_s);
  std::vector<RankedCandidate> ranked;
  ranked.reserve(candidates.size());
  std::vector<std::vector<std::size_t>> by_fastest(curve_count);
  std::vector<std::vector<std::size_t>> by_slowest(curve_count);
  for (CandidatePath& candidate : candidates) {
    bool is_dominated = false;
    for (std::size_t slowest = 0;
         slowest <= fastest_rank(candidate);
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
    for (std::size_t fastest = slowest_rank(candidate);
         fastest < curve_count; ++fastest) {
      const auto& indexes = by_fastest[fastest];
      // candidates were sorted by first_leg before these append-only indexes
      // were built. Inactive entries retain the same distance and position.
      // Every entry before this suffix fails the original first-leg predicate;
      // use that predicate verbatim, not a rounded first_leg - tolerance key.
      const auto first = catalog_first_leg_index_enabled_
          ? std::partition_point(indexes.begin(), indexes.end(),
                [&](std::size_t index) {
                  return !(candidate.first_leg <=
                      ranked[index].candidate.first_leg + kCostTolerance);
                })
          : indexes.begin();
      for (auto position = first; position != indexes.end(); ++position) {
        const std::size_t index = *position;
        if (ranked[index].active &&
            dominates(candidate, ranked[index].candidate)) {
          ranked[index].active = false;
        }
      }
    }
    const std::size_t index = ranked.size();
    const std::size_t fastest = fastest_rank(candidate);
    const std::size_t slowest = slowest_rank(candidate);
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

  // Surviving views refer to complete_storage or shortest_station_paths_,
  // not these consumed sorting/ranking buffers.  Release their allocations
  // before constructing the additional station vectors and charging labels.
  trace.enter("workspace_release", &PathCatalogStatistics::workspace_release_s);
  std::vector<CandidatePath>().swap(candidates);
  std::vector<RankedCandidate>().swap(ranked);
  std::vector<std::vector<std::size_t>>().swap(by_fastest);
  std::vector<std::vector<std::size_t>>().swap(by_slowest);

  trace.enter("descriptor_materialization",
              &PathCatalogStatistics::descriptor_materialization_s);
  auto result = std::make_shared<PathCatalog>();
  result->packed = packed_catalogs_enabled_ &&
      instance_.station_ids().size() <= result->packed_station_ids.size();
  if (result->packed) {
    result->packed_station_count =
        static_cast<std::uint8_t>(instance_.station_ids().size());
    std::copy(instance_.station_ids().begin(), instance_.station_ids().end(),
              result->packed_station_ids.begin());
  }
  const auto maximum_offset = std::numeric_limits<std::uint32_t>::max();
  if (nondominated.size() >= maximum_offset) {
    throw std::length_error("station catalog descriptor index overflow");
  }
  std::size_t station_count = 0;
  for (const CandidatePath& candidate : nondominated) {
    const auto size = packed_construction
        ? result->path_size(PathDescriptor{candidate.packed_station_path->stations, 0.0})
        : candidate.station_path->stations.size();
    if (size > maximum_offset - station_count) {
      throw std::length_error("station catalog flat offset overflow");
    }
    station_count += size;
  }
  result->paths.reserve(nondominated.size() + 1);
  if (!result->packed) result->station_ids.reserve(station_count);
  result->paths.push_back(PathDescriptor{
      .station_storage = 0,
      .total_distance_km = instance_.distance_km(origin_id, target_id)});
  for (const CandidatePath& candidate : nondominated) {
    if (packed_construction) {
      result->paths.push_back(PathDescriptor{
          .station_storage = candidate.packed_station_path->stations,
          .total_distance_km = candidate.total_distance});
    } else {
      const auto& stations = candidate.station_path->stations;
      result->paths.push_back(PathDescriptor{
          .station_storage = result->encode_path(
              stations, static_cast<std::uint32_t>(result->station_ids.size())),
          .total_distance_km = candidate.total_distance,
      });
      if (!result->packed) {
        result->station_ids.insert(result->station_ids.end(), stations.begin(),
                                   stations.end());
      }
    }
  }
  trace.enter("option_sort", &PathCatalogStatistics::option_sort_s);
  std::sort(result->paths.begin(), result->paths.end(),
            [&](const PathDescriptor& left, const PathDescriptor& right) {
              const auto left_key =
                  std::pair{left.total_distance_km, result->path_size(left)};
              const auto right_key =
                  std::pair{right.total_distance_km, result->path_size(right)};
              if (left_key != right_key) return left_key < right_key;
              return result->lexically_less(left, right);
            });
  result->lexical_order.resize(result->paths.size());
  for (std::size_t index = 0; index < result->paths.size(); ++index) {
    result->lexical_order[index] = static_cast<std::uint32_t>(index);
  }
  std::sort(result->lexical_order.begin(), result->lexical_order.end(),
            [&](std::uint32_t left, std::uint32_t right) {
              return result->lexically_less(result->paths[left],
                                             result->paths[right]);
            });
  const std::size_t evictions_before = catalog_cache_.evictions();
  trace.enter("cache_admission", &PathCatalogStatistics::cache_admission_s);
  disk_cache_store(key, *result);
  const std::size_t bytes = result->retained_bytes();
  if (bytes <= catalog_cache_budget_) {
    while (catalog_cache_.size() != 0 &&
           catalog_cache_bytes_ + bytes > catalog_cache_budget_) {
      const auto evicted = catalog_cache_.discard_oldest();
      catalog_cache_bytes_ -= (*evicted)->retained_bytes();
    }
    catalog_cache_.put(key, result);
    catalog_cache_bytes_ += bytes;
  }
  if (timing != nullptr) {
    timing->path_option_cache_evictions +=
        catalog_cache_.evictions() - evictions_before;
    timing->path_option_cache_miss_build_s += elapsed_seconds(build_started);
  }
  trace.enter("cleanup", &PathCatalogStatistics::cleanup_s);
  return result;
}

std::shared_ptr<const PathSampler::PathOptions> PathSampler::path_options(
    int origin_id, int target_id, RepairTimingStatistics* timing,
    PathCatalogStatistics* statistics) {
  const auto catalog = path_catalog(origin_id, target_id, timing, statistics);
  const AnchorPair key{origin_id, target_id};
  if (auto cached = option_cache_.get(key); cached.has_value()) return *cached;
  const auto started = std::chrono::steady_clock::now();
  auto options = std::make_shared<PathOptions>();
  options->reserve(catalog->size());
  NonlinearSequenceEvaluator evaluator(instance_);
  PathCatalog::DecodeBuffer decoded;
  CatalogReader reader(*catalog);
  for (std::size_t index = 0; index < catalog->size(); ++index) {
    const PathDescriptor descriptor = reader.descriptor(index);
    const auto path = catalog->stations(descriptor, decoded);
    options->push_back(PathOption{
        .stations = std::vector<int>(path.begin(), path.end()),
        .total_distance_km = descriptor.total_distance_km,
        .label = evaluator.edge_label(origin_id, path, target_id)});
  }
  const auto retained_bytes = [](const PathOptions& values) {
    std::size_t bytes = sizeof(PathOptions) + values.capacity() * sizeof(PathOption) + 256;
    for (const PathOption& option : values) {
      bytes += option.stations.capacity() * sizeof(int) +
          option.label.station_legs.capacity() * sizeof(ChargingSequenceLeg) + 64;
    }
    return bytes;
  };
  constexpr std::size_t label_cache_budget = 128 * 1024 * 1024;
  const std::size_t bytes = retained_bytes(*options);
  if (bytes <= label_cache_budget) {
    while (option_cache_.size() != 0 && option_cache_bytes_ + bytes > label_cache_budget) {
      const auto evicted = option_cache_.discard_oldest();
      option_cache_bytes_ -= retained_bytes(**evicted);
    }
    option_cache_.put(key, options);
    option_cache_bytes_ += bytes;
  }
  const double elapsed = elapsed_seconds(started);
  if (statistics != nullptr) statistics->label_materialization_s += elapsed;
  if (timing != nullptr) timing->path_option_cache_miss_build_s += elapsed;
  return options;
}

std::vector<std::vector<int>> PathSampler::paths_between(int origin_id,
    int target_id, PathCatalogStatistics* statistics) {
  const auto catalog = path_catalog(origin_id, target_id, nullptr, statistics);
  const auto export_started = statistics != nullptr
      ? std::chrono::steady_clock::now()
      : std::chrono::steady_clock::time_point{};
  std::vector<std::vector<int>> paths;
  paths.reserve(catalog->size());
  PathCatalog::DecodeBuffer decoded;
  CatalogReader reader(*catalog);
  for (std::size_t ordinal = 0; ordinal < catalog->size(); ++ordinal) {
    const auto path = catalog->stations(
        reader.descriptor(reader.lexical_index(ordinal)), decoded);
    paths.emplace_back(path.begin(), path.end());
  }
  if (statistics != nullptr) {
    statistics->export_s = elapsed_seconds(export_started);
  }
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

bool PathSampler::path_dominated(int origin, int target,
                                 const std::vector<int>& path) const {
  if (path.empty()) return false;  // Direct paths are never RP-dominated.
  const auto& stations = instance_.station_ids();
  const double first = instance_.distance_km(origin, path.front());
  const double last = instance_.distance_km(path.back(), target);
  double total = last;
  int previous = origin;
  std::size_t fastest = curve_rank_by_type_.size(), slowest = 0;
  for (const int node : path) {
    total += instance_.distance_km(previous, node);
    const auto rank = curve_rank_by_type_.at(instance_.node(node).station_type);
    fastest = std::min(fastest, rank);
    slowest = std::max(slowest, rank);
    previous = node;
  }
  const auto reference_key = std::tie(first, last, total, fastest, slowest, path);
  std::vector<bool> used(stations.size(), false);
  std::vector<int> witness;
  const auto search = [&](auto&& self, int from, double distance,
                           double witness_first, std::size_t best,
                           std::size_t worst) -> bool {
    for (std::size_t index = 0; index < stations.size(); ++index) {
      if (used[index] || station_curve_ranks_[index] > fastest) continue;
      const int next = stations[index];
      const double leg = instance_.distance_km(from, next);
      if (instance_.energy_wh(from, next) >
          instance_.vehicle().battery_capacity_wh + kEnergyToleranceWh) continue;
      const double leading = witness.empty() ? leg : witness_first;
      if (leading > first + kCostTolerance) continue;
      const double prefix = distance + leg;
      const double trailing = instance_.distance_km(next, target);
      if (prefix + trailing > total + kCostTolerance) continue;
      const std::size_t next_best = witness.empty() ? station_curve_ranks_[index]
                                                    : std::min(best, station_curve_ranks_[index]);
      const std::size_t next_worst = witness.empty() ? station_curve_ranks_[index]
                                                     : std::max(worst, station_curve_ranks_[index]);
      used[index] = true;
      witness.push_back(next);
      const double witness_total = prefix + trailing;
      if (trailing <= last + kCostTolerance &&
          instance_.energy_wh(next, target) <=
              instance_.vehicle().battery_capacity_wh + kEnergyToleranceWh &&
          witness != path) {
        const bool mutual = slowest <= next_best &&
            first <= leading + kCostTolerance && last <= trailing + kCostTolerance &&
            total <= witness_total + kCostTolerance;
        if (!mutual || std::tie(leading, trailing, witness_total,
                               next_best, next_worst, witness) < reference_key) return true;
      }
      if (self(self, next, prefix, leading, next_best, next_worst)) return true;
      witness.pop_back();
      used[index] = false;
    }
    return false;
  };
  return search(search, origin, 0.0, 0.0, 0, 0);
}

std::optional<std::vector<int>> PathSampler::sample_implicit_path(
    int origin, int target, const std::vector<int>& current) {
  const auto& stations = instance_.station_ids();
  std::vector<bool> used(stations.size(), false);
  std::vector<int> path;
  int from = origin;
  double distance = 0.0;
  // Class A: geometric stopping, uniform feasible next station. Every finite
  // elementary path has positive proposal probability; no hop/Top-K cap.
  if (!current.empty() && std::bernoulli_distribution(0.5)(random_)) return path;
  while (path.size() < stations.size()) {
    std::vector<std::size_t> reachable;
    for (std::size_t index = 0; index < stations.size(); ++index) {
      if (!used[index] && instance_.energy_wh(from, stations[index]) <=
          instance_.vehicle().battery_capacity_wh + kEnergyToleranceWh)
        reachable.push_back(index);
    }
    if (reachable.empty()) return std::nullopt;
    const std::size_t next = reachable[random_index(reachable, random_)];
    distance += instance_.distance_km(from, stations[next]);
    from = stations[next];
    used[next] = true;
    path.push_back(from);
    const double tail = instance_.distance_km(from, target);
    const double route_distance = instance_.distance_km(instance_.depot().id, origin) +
        distance + tail + instance_.distance_km(target, instance_.depot().id);
    const double lower_time = route_distance / instance_.vehicle().speed_kmph +
        std::max(0.0, route_distance * instance_.vehicle().consumption_wh_per_km -
                          instance_.vehicle().battery_capacity_wh) /
            maximum_charging_rate_wh_per_h_ + instance_.node(origin).service_time_h +
        (origin == target ? 0.0 : instance_.node(target).service_time_h);
    if (lower_time > instance_.vehicle().max_route_duration_h + kTimeToleranceH)
      return std::nullopt;
    if (instance_.energy_wh(from, target) <=
            instance_.vehicle().battery_capacity_wh + kEnergyToleranceWh &&
        (path.size() == stations.size() || std::bernoulli_distribution(0.5)(random_))) {
      if (path == current || path_dominated(origin, target, path)) return std::nullopt;
      return path;
    }
  }
  return std::nullopt;
}

std::optional<PathSampler::RandomReplacementProposal>
PathSampler::random_replacement_proposal(
    const Solution& solution) {
  std::size_t gap_count = 0;
  for (const Plan& plan : solution.plans) gap_count += plan.gaps.size();
  if (gap_count == 0) return std::nullopt;
  std::optional<std::mt19937_64> before_catalogue_read;
  if (catalog_disk_cache_.enabled && !implicit_replacement_)
    before_catalogue_read.emplace(random_);
  std::size_t selected_ordinal =
      std::uniform_int_distribution<std::size_t>(0, gap_count - 1)(random_);
  std::size_t selected_route = 0;
  while (selected_ordinal >= solution.plans[selected_route].gaps.size()) {
    selected_ordinal -= solution.plans[selected_route].gaps.size();
    ++selected_route;
  }
  const std::size_t selected_gap = selected_ordinal;
  const Plan& plan = solution.plans[selected_route];
  std::vector<int> anchors;
  anchors.reserve(plan.customer_ids.size() + 2);
  anchors.push_back(instance_.depot().id);
  anchors.insert(anchors.end(), plan.customer_ids.begin(),
                 plan.customer_ids.end());
  anchors.push_back(instance_.depot().id);
  auto gaps = plan.gaps;
  if (implicit_replacement_) {
    auto sampled = sample_implicit_path(anchors[selected_gap],
        anchors[selected_gap + 1], plan.gaps[selected_gap]);
    if (!sampled.has_value()) return std::nullopt;
    gaps[selected_gap] = std::move(*sampled);
  } else {
    std::optional<std::vector<int>> sampled;
    try {
      sampled = sample_explicit_path(anchors[selected_gap], anchors[selected_gap + 1],
                                     plan.gaps[selected_gap]);
    } catch (...) {
      // A failed catalogue read must also roll back this operation's earlier
      // gap-selection draw. Successful and null proposals retain the old RNG.
      if (before_catalogue_read) random_ = *before_catalogue_read;
      throw;
    }
    if (!sampled.has_value()) return std::nullopt;
    gaps[selected_gap] = std::move(*sampled);
  }
  return RandomReplacementProposal{
      .route_index = selected_route,
      .plan = factory_.make_plan(plan.customer_ids, std::move(gaps))};
}

// Compatibility callers still receive a complete candidate. The proposal
// owns only its changed Plan; no unchanged route is copied for a null draw.
std::optional<Solution> PathSampler::random_replacement(
    const Solution& solution) {
  auto proposal = random_replacement_proposal(solution);
  if (!proposal.has_value()) return std::nullopt;
  Solution candidate = solution;
  candidate.plans[proposal->route_index] = std::move(proposal->plan);
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
  const auto setup_started = statistics != nullptr
      ? std::chrono::steady_clock::now()
      : std::chrono::steady_clock::time_point{};
  std::vector<int> anchors;
  anchors.reserve(plan.customer_ids.size() + 2);
  anchors.push_back(instance_.depot().id);
  anchors.insert(anchors.end(), plan.customer_ids.begin(),
                 plan.customer_ids.end());
  anchors.push_back(instance_.depot().id);
  NonlinearSequenceEvaluator sequence_evaluator(instance_);
  std::vector<ChargingSequenceLabel> edges;
  std::vector<ChargingSequenceLabel> prefixes;
  std::vector<ChargingSequenceLabel> suffixes;
  bool labels_prepared = false;
  const auto prepare_labels = [&] {
    if (labels_prepared) return true;
    edges.reserve(plan.gaps.size());
    for (std::size_t gap_index = 0; gap_index < plan.gaps.size(); ++gap_index) {
      if (deadline_reached(deadline)) return false;
      edges.push_back(sequence_evaluator.edge_label(
          anchors[gap_index], plan.gaps[gap_index], anchors[gap_index + 1]));
    }
    for (const std::size_t gap : eligible_gaps) {
      if (gap >= plan.gaps.size()) {
        throw std::out_of_range("eligible path-repair gap is invalid");
      }
    }
    prefixes.resize(plan.gaps.size() + 1);
    for (std::size_t gap_index = 0; gap_index < plan.gaps.size(); ++gap_index) {
      if (deadline_reached(deadline)) return false;
      prefixes[gap_index + 1] = sequence_evaluator.concatenate(
          prefixes[gap_index], edges[gap_index]);
    }
    suffixes.resize(plan.gaps.size() + 1);
    for (std::size_t gap_index = plan.gaps.size(); gap_index-- > 0;) {
      if (deadline_reached(deadline)) return false;
      suffixes[gap_index] = sequence_evaluator.concatenate(
          edges[gap_index], suffixes[gap_index + 1]);
    }
    labels_prepared = true;
    return true;
  };
  // The full physical interior identifies customer/gap metadata only for a
  // canonical plan. Partial customer partitions are valid here. Timed calls
  // retain the original eager setup and every existing deadline check.
  const bool bound_fp_supported = supported_bound_fp_environment();
  if (!replacement_summary_statistics_.fp_fallback &&
      replacement_summary_statistics_.bound_certified && !bound_fp_supported) {
    // Preserve the permanent fallback after the certified FP environment
    // changes; restoring the FP mode must not silently re-enable bounds.
    replacement_summary_statistics_.fp_fallback = true;
  }
  bool cacheable = !replacement_summary_statistics_.fp_fallback &&
      !deadline.has_value() && statistics == nullptr &&
      replacement_summary_budget_ != 0 &&
      plan.route.visits.size() <= 256 &&
      plan.gaps.size() == plan.customer_ids.size() + 1;
  std::size_t visit = 0;
  if (cacheable) {
    const auto matches = [&](int id) {
      return visit < plan.route.visits.size() && plan.route.visits[visit++] == id;
    };
    for (std::size_t gap = 0; gap < plan.gaps.size(); ++gap) {
      for (int station : plan.gaps[gap]) {
        if (instance_.node(station).kind != NodeKind::Station || !matches(station)) {
          cacheable = false;
        }
      }
      if (gap < plan.customer_ids.size()) {
        const int customer = plan.customer_ids[gap];
        if (instance_.node(customer).kind != NodeKind::Customer || !matches(customer) ||
            std::find(plan.customer_ids.begin(), plan.customer_ids.begin() + gap,
                      customer) != plan.customer_ids.begin() + gap) {
          cacheable = false;
        }
      }
    }
    cacheable = cacheable && visit == plan.route.visits.size();
  }
  if (!prepare_labels()) return std::nullopt;
  for (const std::size_t gap_index : eligible_gaps) {
    if (gap_index >= plan.gaps.size()) {
      throw std::out_of_range("eligible path-repair gap is invalid");
    }
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
  const bool bound_enabled = cacheable &&
      replacement_summary_statistics_.bound_certified &&
      bound_fp_supported;
  const bool shortfall_bound_enabled = bound_enabled &&
      shortfall_replacement_bound_enabled_ &&
      NonlinearSequenceEvaluator::bounded_evaluation_compatible() &&
      nonlinear_travel_bound(0.0, 0.0, 0, instance_.vehicle(), penalty_lambda).has_value();
  // Already computed once for this immutable Instance. Under the certified
  // data/FP scope its finite conservative upper rate excludes original Inf
  // marginal rates; do not rescan curves for each physical candidate.
  const bool all_charging_rates_finite =
      std::isfinite(maximum_charging_rate_wh_per_h_);
  if (statistics != nullptr) {
    statistics->setup_s += elapsed_seconds(setup_started);
  }
  for (const std::size_t gap_index : eligible_gaps) {
    if (deadline_reached(deadline)) break;
    const auto catalog_started = statistics != nullptr
        ? std::chrono::steady_clock::now()
        : std::chrono::steady_clock::time_point{};
    const std::size_t cache_hits_before = option_cache_hits();
    const std::size_t cache_misses_before = option_cache_misses();
    // Untimed completion builds option labels on demand from descriptors.
    // Timed/profile requests retain eager option labels and deadline checks.
    // Both scans use the same catalog descriptor order.
    const bool scan_descriptors = !deadline.has_value() && statistics == nullptr;
    const auto catalog = scan_descriptors
        ? path_catalog(anchors[gap_index], anchors[gap_index + 1])
        : std::shared_ptr<const PathCatalog>{};
    const auto options = scan_descriptors
        ? std::shared_ptr<const PathOptions>{}
        : path_options(anchors[gap_index], anchors[gap_index + 1]);
    std::optional<CatalogReader> reader;
    if (catalog) reader.emplace(*catalog);
    const std::size_t option_count = scan_descriptors
        ? catalog->size() : options->size();
    if (statistics != nullptr) {
      statistics->option_catalog_s += elapsed_seconds(catalog_started);
      statistics->option_cache_hits +=
          option_cache_hits() - cache_hits_before;
      statistics->option_cache_misses +=
          option_cache_misses() - cache_misses_before;
    }
    if (deadline_reached(deadline)) break;
    const auto scoring_started = statistics != nullptr
        ? std::chrono::steady_clock::now()
        : std::chrono::steady_clock::time_point{};
    double option_service = 0.0;
    if (bound_enabled) {
      const auto& target = instance_.node(anchors[gap_index + 1]);
      if (target.kind == NodeKind::Customer) option_service = target.service_time_h;
    }
    bool bulk_bound = bound_enabled && bulk_replacement_bound_enabled_;
    std::size_t tail_station_upper = 0;
    double tail_service = 0.0;
    if (bulk_bound) {
      if (!prepare_labels()) return std::nullopt;
      const auto& prefix = prefixes[gap_index];
      const auto& suffix = suffixes[gap_index + 1];
      tail_station_upper = prefix.station_legs.size() +
          instance_.station_ids().size() + suffix.station_legs.size();
      bulk_bound = tail_station_upper <= 62;
      tail_service = prefix.initialized
          ? prefix.service_time_h + option_service : option_service;
      if (suffix.initialized) tail_service += suffix.service_time_h;
    }
    for (std::size_t option_index = 0; option_index < option_count; ++option_index) {
      const PathDescriptor descriptor = scan_descriptors
          ? reader->descriptor(option_index) : PathDescriptor{};
      if (bulk_bound) {
        const auto& prefix = prefixes[gap_index];
        const auto& suffix = suffixes[gap_index + 1];
        const double option_travel = catalog_travel_lower_bound(
            descriptor.total_distance_km,
            instance_.vehicle().speed_kmph);
        double travel = prefix.initialized
            ? prefix.travel_time_h + option_travel : option_travel;
        if (suffix.initialized) travel += suffix.travel_time_h;
        auto lower = nonlinear_travel_bound(
                travel, tail_service, tail_station_upper, instance_.vehicle(),
                penalty_lambda);
        if (lower && complete_paths_ && shortfall_bound_enabled && all_charging_rates_finite &&
            zero_shift_energy_bound_enabled_) {
          if (const auto energy_lower = zero_shift_energy_cost_bound(
                  travel, tail_service, instance_.vehicle(), penalty_lambda))
            *lower = std::max(*lower, *energy_lower);
        }
        if (lower) {
          ++replacement_summary_statistics_.tail_checks;
          if (*lower >= incumbent_cost - kCostTolerance ||
              *lower > minimum_cost + kCostTolerance) {
            // The cutoff is already certified. Count raw descriptors, including
            // the cutoff descriptor and any current/excluded path occurrences.
            // This intentionally differs from the former eligible-score metric:
            // diagnostic accounting must not issue lexical catalogue queries.
            replacement_summary_statistics_.raw_tail_descriptors_skipped +=
                option_count - option_index;
            break;
          }
        }
      }
      PathCatalog::DecodeBuffer decoded;
      const std::span<const int> stations = scan_descriptors
          ? catalog->stations(descriptor, decoded)
          : std::span<const int>{(*options)[option_index].stations};
      if (deadline_reached(deadline)) break;
      if (std::ranges::equal(stations, plan.gaps[gap_index])) continue;
      if (std::any_of(excluded.begin(), excluded.end(),
                      [&](const ExcludedPathOption& item) {
                        return item.gap == gap_index &&
                               std::ranges::equal(item.stations, stations);
      })) {
        continue;
      }
      if (statistics != nullptr) ++statistics->options_scored;
      RouteCostSummary summary;
      if (!prepare_labels()) return std::nullopt;
      bool candidate_shortfall_bound = false;
      if (bound_enabled) {
        const auto& prefix = prefixes[gap_index];
        const auto& suffix = suffixes[gap_index + 1];
        const std::size_t station_count = prefix.station_legs.size() +
            stations.size() + suffix.station_legs.size();
        if (station_count != 0 && station_count <= 62) {
          double travel;
          // Reproduce edge_label/append scalar arithmetic exactly. In
          // particular, do not substitute total_distance / speed, add an
          // empty prefix, or regroup the precomposed suffix.
          double option_travel = 0.0;
          int previous = anchors[gap_index];
          for (const int station : stations) {
            option_travel += instance_.travel_time_h(previous, station);
            previous = station;
          }
          option_travel += instance_.travel_time_h(
              previous, anchors[gap_index + 1]);
          travel = prefix.initialized
              ? prefix.travel_time_h + option_travel : option_travel;
          if (suffix.initialized) travel += suffix.travel_time_h;
          double service = prefix.initialized
              ? prefix.service_time_h + option_service : option_service;
          if (suffix.initialized) service += suffix.service_time_h;
          // Canonical <=256 visits plus <=62 candidate CS and certified
          // pair energy <=1e6 bound every energy segment well below 1e13.
          // Original prefix/suffix validation and catalog CS identity remain.
          if (const auto lower = nonlinear_travel_bound(
                  travel, service, station_count, instance_.vehicle(),
                  penalty_lambda)) {
            candidate_shortfall_bound = shortfall_bound_enabled;
            ++replacement_summary_statistics_.bound_checks;
            if (*lower >= incumbent_cost - kCostTolerance ||
                *lower > minimum_cost + kCostTolerance) {
              ++replacement_summary_statistics_.bound_pruned;
              continue;
            }
          }
        }
      }
      ChargingSequenceLabel candidate = prefixes[gap_index];
      if (scan_descriptors) {
        const ChargingSequenceLabel option_label = sequence_evaluator.edge_label(
            anchors[gap_index], stations, anchors[gap_index + 1]);
        ++replacement_summary_statistics_.labels_built;
        sequence_evaluator.append(candidate, option_label);
      } else {
        sequence_evaluator.append(candidate, (*options)[option_index].label);
      }
      sequence_evaluator.append(candidate, suffixes[gap_index + 1]);
      if (candidate_shortfall_bound) {
        auto evaluated = sequence_evaluator.evaluate_bounded(
            candidate, penalty_lambda, incumbent_cost - kCostTolerance,
            minimum_cost + kCostTolerance, all_charging_rates_finite);
        ++replacement_summary_statistics_.misses;
        if (std::get_if<NonlinearCostLowerBound>(&evaluated) != nullptr) {
          // A partial bound is not an exact score and must not enter ranking.
          ++replacement_summary_statistics_.propagation_pruned;
          continue;
        }
        summary = std::get<RouteCostSummary>(evaluated);
      } else {
        summary = sequence_evaluator.evaluate(candidate);
        ++replacement_summary_statistics_.misses;
      }
      const double candidate_cost = generalized_cost(summary, instance_, penalty_lambda);
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
            .stations = std::vector<int>(stations.begin(), stations.end())});
      }
    }
    if (statistics != nullptr) {
      statistics->fast_scoring_s += elapsed_seconds(scoring_started);
    }
  }
  if (deadline_reached(deadline)) return std::nullopt;
  if (tolerance_band.empty()) return std::nullopt;

  const auto materialization_started = statistics != nullptr
      ? std::chrono::steady_clock::now()
      : std::chrono::steady_clock::time_point{};
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
      const auto acceptance_started = statistics != nullptr
          ? std::chrono::steady_clock::now()
          : std::chrono::steady_clock::time_point{};
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
    std::optional<std::chrono::steady_clock::time_point> deadline,
    bool single_step) {
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
  if (single_step) {
    if (introduced.empty()) return plan;
    // Keep the original catalog, lazy labels, certified bounds and proposal
    // tie ordering. No fixed-point iteration or full-rejection exclusion retry.
    auto proposed = best_replacement_in_gaps(
        plan, penalty_lambda, introduced, {}, deadline, nullptr);
    if (proposed &&
        generalized_cost(*proposed->evaluation, instance_, penalty_lambda) <
            generalized_cost(*plan.evaluation, instance_, penalty_lambda) -
                kCostTolerance) {
      return std::move(*proposed);
    }
    return plan;
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
  const auto total_started = statistics != nullptr
      ? std::chrono::steady_clock::now()
      : std::chrono::steady_clock::time_point{};
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
                                               std::size_t granular_threshold,
                                               bool rebuild_station_paths,
                                               bool preserve_station_visits)
    : instance_(instance),
      factory_(factory),
      path_sampler_(path_sampler),
      random_(random),
      profile_(profile),
      rebuild_station_paths_(rebuild_station_paths),
      preserve_station_visits_(preserve_station_visits),
      granular_(profile == SearchProfile::GranularAdaptive &&
                        instance.customer_ids().size() > granular_threshold
                    ? std::make_unique<GranularNeighborhood>(
                          instance, promising_arc_count)
                    : nullptr) {
  if ((rebuild_station_paths_ || preserve_station_visits_) &&
      profile_ != SearchProfile::LinkedInsertionCandidates) {
    throw std::invalid_argument("station initialization requires candidate completion");
  }
  if (rebuild_station_paths_ && preserve_station_visits_) {
    throw std::invalid_argument("rebuilt and physical initialization are exclusive");
  }
}

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
  if (preserve_station_visits_) {
    // Random first customer is unchanged. The one shuffle is the deterministic
    // base continuation and tie order for initialization-only fortified rollout.
    std::shuffle(unassigned.begin(), unassigned.end(), random_);
    std::unordered_map<int, Plan> singletons;
    std::unordered_map<int, double> singleton_bounds;
    NonlinearSequenceEvaluator direct_scorer(instance_);
    for (const int customer : unassigned) {
      singleton_bounds.emplace(customer, linked_route_lower_bound(
          direct_scorer.evaluate(direct_scorer.route_label(Route{{customer}})),
          instance_, penalty_lambda, &path_sampler_));
    }
    // Fixed instance/lambda/direct-rebuild policy throughout this construction.
    // Detached trials share deterministic completion results, never RNG draws.
    std::unordered_map<std::vector<int>, Plan, IntVectorHash> rebuilt_sequences;
    const auto insert_best = [&](const Solution& source, int customer) {
      ++statistics_.physical_insertion_decisions;
      std::optional<InsertionMove> best;
      std::optional<Plan> best_plan;
      bool best_non_retained = true;
      // Exact-cost ties preserve existing charging structure across placements.
      // This lexicographic key changes no strict generalized-cost comparison.
      const auto consider = [&](InsertionMove move, const Plan& plan,
                                bool non_retained) {
        const auto rank = std::tuple{move.delta, non_retained, move.customer_id,
                                     move.route_index, move.position};
        if (!best || rank < std::tuple{best->delta, best_non_retained,
                                      best->customer_id, best->route_index,
                                      best->position}) {
          best = std::move(move);
          best_plan = plan;
          best_non_retained = non_retained;
        }
      };
      for (std::size_t route = 0; route < source.plans.size(); ++route) {
        const double old_cost = generalized_cost(
            *source.plans[route].evaluation, instance_, penalty_lambda);
        const std::size_t slots = source.plans[route].route.visits.size() + 1;
        for (std::size_t slot = 0; slot < slots; ++slot) {
          auto visits = source.plans[route].route.visits;
          visits.insert(visits.begin() + static_cast<std::ptrdiff_t>(slot), customer);
          Plan retained = factory_.plan_visits(std::move(visits));
          const double retained_cost = generalized_cost(
              *retained.evaluation, instance_, penalty_lambda);
          auto completed = rebuilt_sequences.find(retained.customer_ids);
          if (completed == rebuilt_sequences.end()) {
            ++statistics_.completion_memo_misses;
            ++statistics_.candidates_evaluated;
            Plan rebuilt = path_sampler_.repair_new_gaps(
                factory_.plan_sequence(retained.customer_ids), {}, penalty_lambda);
            completed = rebuilt_sequences.emplace(
                retained.customer_ids, std::move(rebuilt)).first;
          } else {
            ++statistics_.completion_memo_hits;
          }
          const double rebuilt_cost = generalized_cost(
              *completed->second.evaluation, instance_, penalty_lambda);
          const bool use_rebuilt = rebuilt_cost < retained_cost;
          const Plan& selected = use_rebuilt ? completed->second : retained;
          InsertionMove move;
          move.customer_id = customer;
          move.route_index = route;
          move.position = slot;
          move.delta = (use_rebuilt ? rebuilt_cost : retained_cost) - old_cost;
          ++statistics_.physical_insertion_candidates;
          consider(std::move(move), selected, use_rebuilt);
        }
      }
      ++statistics_.physical_insertion_candidates;
      if (singleton_bounds.at(customer) <= best->delta + kLinkedBoundSafety) {
        auto cached = singletons.find(customer);
        if (cached == singletons.end()) {
          cached = singletons.emplace(customer, path_sampler_.repair_new_gaps(
              factory_.plan_sequence({customer}), {}, penalty_lambda)).first;
        }
        InsertionMove move;
        move.customer_id = customer;
        move.route_index = source.plans.size();
        move.position = 0;
        move.delta = generalized_cost(*cached->second.evaluation, instance_, penalty_lambda);
        if (singleton_bounds.at(customer) > move.delta + kLinkedBoundSafety) {
          throw std::logic_error("singleton lower bound exceeded completed cost");
        }
        // A completed singleton has no incumbent physical structure to retain.
        consider(std::move(move), cached->second, true);
      }
      if (!best_plan) throw std::logic_error("missing retained/rebuilt rollout winner");
      Solution result = source;
      if (best->route_index == result.plans.size()) {
        result.plans.push_back(std::move(*best_plan));
      } else {
        result.plans[best->route_index] = std::move(*best_plan);
      }
      return result;
    };
    struct ConstructionAction {
      int customer{};
      Solution before;
      Solution after;
    };
    struct ConstructionContinuation {
      std::vector<ConstructionAction> actions;
      double cost{};
    };
    const auto complete = [&](Solution state, const std::vector<int>& order) {
      ConstructionContinuation result;
      result.actions.reserve(order.size());
      for (const int customer : order) {
        Solution after = insert_best(state, customer);
        result.actions.push_back({customer, state, after});
        state = std::move(after);
      }
      state.validate_partition(instance_);
      result.cost = generalized_cost(state, instance_, penalty_lambda);
      return result;
    };
    const auto same_physical_state = [](const Solution& left, const Solution& right) {
      if (left.plans.size() != right.plans.size()) return false;
      for (std::size_t i = 0; i < left.plans.size(); ++i) {
        if (left.plans[i].route.visits != right.plans[i].route.visits) return false;
      }
      return true;
    };
    // Default policy: finite full breadth, with no new time or restart budget.
    // Only a fully returned continuation can replace the saved base fallback.
    // Preserve the complete actual one-action rollout as the deeper policy's
    // incumbent, rather than only its less informed shuffled base completion.
    const auto fortify = [&](Solution state, std::vector<int> remaining,
                             ConstructionContinuation incumbent,
                             bool use_root_pair) {
      std::vector<ConstructionAction> committed;
      committed.reserve(remaining.size());
      bool first_decision = use_root_pair;
      while (!remaining.empty()) {
        const auto consider_order = [&](const std::vector<int>& order) {
          bool same_order = order.size() == incumbent.actions.size();
          for (std::size_t i = 0; same_order && i < order.size(); ++i) {
            same_order = order[i] == incumbent.actions[i].customer;
          }
          if (same_order) return;
          ConstructionContinuation challenger = complete(state, order);
          if (challenger.cost < incumbent.cost - kCostTolerance) {
            incumbent = std::move(challenger);
          }
        };
        for (const int customer : remaining) {
          if (first_decision && remaining.size() > 1) {
            // Only the first post-random decision looks ahead two customer choices.
            // Every trial still uses true cheapest placement and the same base suffix.
            for (const int second : remaining) {
              if (second == customer) continue;
              std::vector<int> order{customer, second};
              order.reserve(remaining.size());
              for (const int other : remaining) {
                if (other != customer && other != second) order.push_back(other);
              }
              consider_order(order);
            }
          } else {
            std::vector<int> order{customer};
            order.reserve(remaining.size());
            for (const int other : remaining) if (other != customer) order.push_back(other);
            consider_order(order);
          }
        }
        first_decision = false;
        if (incumbent.actions.empty() ||
            !same_physical_state(state, incumbent.actions.front().before)) {
          throw std::logic_error("rollout continuation lost its committed prefix");
        }
        const int customer = incumbent.actions.front().customer;
        ConstructionAction action = std::move(incumbent.actions.front());
        state = action.after;
        committed.push_back(std::move(action));
        incumbent.actions.erase(incumbent.actions.begin());
        std::erase(remaining, customer);
      }
      state.validate_partition(instance_);
      if (std::abs(generalized_cost(state, instance_, penalty_lambda) -
                   incumbent.cost) > kCostTolerance) {
        throw std::logic_error("fortified constructor saved cost mismatch");
      }
      return ConstructionContinuation{std::move(committed), incumbent.cost};
    };
    ConstructionContinuation one_action = fortify(
        current, unassigned, complete(current, unassigned), false);
    const double baseline_cost = one_action.cost;
    ConstructionContinuation selected = fortify(
        current, unassigned, std::move(one_action), true);
    if (selected.cost > baseline_cost + kCostTolerance) {
      throw std::logic_error("deeper construction worsened one-action baseline");
    }
    if (!selected.actions.empty()) {
      current = std::move(selected.actions.back().after);
    }
    current.validate_partition(instance_);
    return current;
  }
  InsertionCompletionMemo completion_memo{
      kInsertionCompletionMemoCapacity};
  VirtualInsertionCompletionMemo virtual_completion_memo{
      kVirtualInsertionCompletionMemoCapacity};
  std::unique_ptr<RouteSequenceScorer> sequence_scorer;
  while (!unassigned.empty()) {
    const PlanFactory::GapMap reference_gaps = rebuild_station_paths_
        ? PlanFactory::GapMap{} : factory_.gap_map(current);
    Solution direct_source;
    if (rebuild_station_paths_) {
      direct_source.plans.reserve(current.plans.size());
      for (const Plan& plan : current.plans) {
        direct_source.plans.push_back(factory_.plan_sequence(plan.customer_ids));
      }
    }
    // Rebuilt candidates may discard any old station detour. Their lower
    // bounds must use direct slices too; an empty gap map alone is insufficient.
    // insertion_moves still subtracts the actual incumbent route cost.
    const Solution& scoring_source = rebuild_station_paths_ ? direct_source : current;
    if (sequence_scorer == nullptr) {
      sequence_scorer = std::make_unique<RouteSequenceScorer>(
          instance_, scoring_source, reference_gaps);
    } else {
      sequence_scorer->reset(scoring_source, reference_gaps);
    }
    const RouteSequenceScorer& scorer = *sequence_scorer;
    std::optional<InsertionMove> best;
    const bool complete_all_candidates =
        completes_insertion_candidates(profile_);
    // Rebuilt initialization may change every gap, outside the two-gap proof.
    const PathSampler* linked_bound_sampler =
        rebuild_station_paths_ ? nullptr : &path_sampler_;
    if (uses_granular_topology(profile_)) {
      std::vector<InsertionMove> candidates;
      std::vector<RejectedInsertionDescriptor> rejected;
      for (const int customer_id : unassigned) {
        auto moves = insertion_moves(
            current, customer_id, instance_, scorer, penalty_lambda, true,
            linked_bound_sampler, false, granular_.get(),
            &statistics_.candidates_rejected_by_granular,
            granular_ != nullptr ? &rejected : nullptr);
        candidates.insert(candidates.end(),
                          std::make_move_iterator(moves.begin()),
                          std::make_move_iterator(moves.end()));
      }
      if (granular_ != nullptr) {
        append_granular_escape(
            candidates, rejected, current, instance_, scorer,
            penalty_lambda, linked_bound_sampler, random_,
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
            true, linked_bound_sampler, false, granular_.get(),
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
                 false, nullptr, false, granular_.get(),
                 &statistics_.candidates_rejected_by_granular)) {
          if (!best.has_value() || insertion_better(move, *best)) {
            best = std::move(move);
          }
        }
      }
    }
    const int inserted = best->customer_id;
    // Materialize the scored completion. Without completion, changed gaps are
    // direct; rebuilt construction starts every gap direct before completion.
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

// One current-state epoch; stores insertion decisions, never output Plans or
// owning evaluation handles. All actual winner materialization is replayed.
struct Perturbation::RepairDecisionMemo {
  struct Step { int customer{}; std::size_t route{}, slot{}; };
  struct Certificate {
    std::weak_ptr<const RouteEvaluation> value;
    bool exact{};
    bool matches(const Plan& p) const {
      const auto locked = value.lock();
      return locked && locked == p.evaluation && exact == p.exact_charging;
    }
  };
  struct AnchorPlan {
    std::vector<int> customers, visits;
    std::vector<std::vector<int>> gaps;
    Certificate certificate;
    bool matches(const Plan& p) const {
      return certificate.matches(p) && customers == p.customer_ids &&
          visits == p.route.visits && gaps == p.gaps;
    }
  };
  struct Entry {
    std::vector<int> removed;
    RepairOperator operation{};
    std::vector<Certificate> certificates;
    std::vector<Step> steps;
    std::uint64_t candidates{}, profiles{};
    std::size_t bytes() const {
      return sizeof(Entry) + removed.capacity()*sizeof(int) +
          certificates.capacity()*(sizeof(Certificate)+sizeof(RouteEvaluation)+64) +
          steps.capacity()*sizeof(Step);
    }
  };
  static constexpr std::size_t capacity = 16, byte_limit = 128*1024;
  std::vector<AnchorPlan> anchor;
  std::vector<Entry> entries;
  std::uint64_t lambda_bits{};
  std::size_t anchor_bytes{};
  bool active{};
  RepairDecisionMemo() { entries.reserve(capacity); }
  void reset() { anchor.clear();entries.clear();anchor_bytes=0;active=false; }
  void begin(const Solution& s, double lambda, InsertionCandidateStatistics& stats) {
    const auto bits=std::bit_cast<std::uint64_t>(lambda);
    if (!std::isfinite(lambda) || lambda<0 || !supported_bound_fp_environment()) {
      reset();return;
    }
    bool same=active && lambda_bits==bits && anchor.size()==s.plans.size();
    for(std::size_t i=0;same && i<anchor.size();++i) same=anchor[i].matches(s.plans[i]);
    if(same)return;
    reset();++stats.repair_memo_epoch_changes;lambda_bits=bits;
    // Upper-bound the source vector payload before copying an oversized state.
    std::size_t estimate=sizeof(*this)+capacity*sizeof(Entry);
    for(const auto& p:s.plans) {
      estimate+=sizeof(AnchorPlan)+sizeof(RouteEvaluation)+64+
          (p.customer_ids.size()+p.route.visits.size())*sizeof(int)+
          p.gaps.size()*sizeof(std::vector<int>);
      for(const auto& g:p.gaps)estimate+=g.size()*sizeof(int);
      if(estimate>byte_limit){++stats.repair_memo_byte_bypasses;return;}
    }
    for(const auto& p:s.plans) anchor.push_back({p.customer_ids,p.route.visits,p.gaps,{p.evaluation,p.exact_charging}});
    anchor_bytes=sizeof(*this)+capacity*sizeof(Entry)+anchor.capacity()*sizeof(AnchorPlan);
    for(const auto& p:anchor) {
      anchor_bytes+=sizeof(RouteEvaluation)+64+(p.customers.capacity()+p.visits.capacity())*sizeof(int)+p.gaps.capacity()*sizeof(std::vector<int>);
      for(const auto& g:p.gaps)anchor_bytes+=g.capacity()*sizeof(int);
    }
    if(anchor_bytes>byte_limit){reset();++stats.repair_memo_byte_bypasses;return;}
    active=true;
  }
  static bool valid_tape(const Entry& e,const DestroyedSolution& d) {
    if(e.steps.size()!=d.removed_customers.size())return false;
    std::vector<std::size_t> lengths,shells;
    for(const auto& p:d.remaining.plans)lengths.push_back(p.route.visits.size());
    for(const auto& s:d.station_shells)shells.push_back(s.size());
    auto pending=d.removed_customers;
    for(const auto& step:e.steps) {
      if(std::find(pending.begin(),pending.end(),step.customer)==pending.end())return false;
      const auto required=std::count_if(shells.begin(),shells.end(),[](auto n){return n!=0;});
      const bool must_activate=pending.size()==static_cast<std::size_t>(required);
      if(step.route<lengths.size()) {
        if(must_activate || step.slot>lengths[step.route])return false;
        ++lengths[step.route];
      } else {
        const auto shell=step.route-lengths.size();
        if(shell<shells.size()) {
          if(shells[shell]==0 || step.slot>shells[shell])return false;
          lengths.push_back(shells[shell]+1);shells[shell]=0;
        } else {
          if(shell!=shells.size() || must_activate || step.slot!=0)return false;
          lengths.push_back(1);
        }
      }
      std::erase(pending,step.customer);
    }
    return pending.empty() && std::none_of(shells.begin(),shells.end(),[](auto n){return n!=0;});
  }
  std::optional<Entry> find(const DestroyedSolution& d,RepairOperator op,InsertionCandidateStatistics& stats) {
    if(!active)return std::nullopt;
    for(std::size_t j=0;j<entries.size();++j) {
      const auto& e=entries[j];
      if(e.operation!=op || e.removed!=d.removed_customers)continue;
      bool same=e.certificates.size()==d.remaining.plans.size();
      for(std::size_t i=0;same && i<e.certificates.size();++i)same=e.certificates[i].matches(d.remaining.plans[i]);
      if(!same || !valid_tape(e,d)) {
        entries.erase(entries.begin()+static_cast<std::ptrdiff_t>(j));
        ++stats.repair_memo_invalid_entries;break;
      }
      Entry result=e;
      std::rotate(entries.begin(),entries.begin()+static_cast<std::ptrdiff_t>(j),entries.begin()+static_cast<std::ptrdiff_t>(j+1));
      ++stats.repair_memo_hits;
      stats.repair_memo_candidates_saved+=result.candidates;
      stats.repair_memo_profiles_saved+=result.profiles;
      return result;
    }
    ++stats.repair_memo_misses;return std::nullopt;
  }
  Entry input(const DestroyedSolution& d,RepairOperator op) const {
    Entry e;e.removed=d.removed_customers;e.operation=op;
    for(const auto& p:d.remaining.plans)e.certificates.push_back({p.evaluation,p.exact_charging});
    e.steps.reserve(e.removed.size());return e;
  }
  void put(Entry entry,InsertionCandidateStatistics& stats) {
    if(!active)return;
    if(anchor_bytes+entry.bytes()>byte_limit){++stats.repair_memo_byte_bypasses;return;}
    std::size_t bytes=anchor_bytes+entry.bytes();for(const auto& e:entries)bytes+=e.bytes();
    while(!entries.empty() && (entries.size()>=capacity || bytes>byte_limit)) {
      bytes-=entries.back().bytes();entries.pop_back();++stats.repair_memo_evictions;
    }
    entries.insert(entries.begin(),std::move(entry));
  }
};

Perturbation::~Perturbation() = default;

Perturbation::Perturbation(const Instance& instance, PlanFactory& factory,
                           PathSampler& path_sampler,
                           std::mt19937_64& random,
                           SearchProfile profile,
                           std::size_t promising_arc_count,
                           std::size_t granular_threshold,
                           PartialRouteSubset partial_route_subset,
                           bool preserve_station_visits)
    : instance_(instance),
      factory_(factory),
      path_sampler_(path_sampler),
      random_(random),
      profile_(profile),
      partial_route_subset_(partial_route_subset),
      preserve_station_visits_(preserve_station_visits),
      granular_(uses_granular_topology(profile) &&
                        instance.customer_ids().size() >
                            granular_threshold
                    ? std::make_unique<GranularNeighborhood>(
                          instance, promising_arc_count)
                    : nullptr) {
  if (preserve_station_visits_ && profile_ != SearchProfile::Paper) {
    throw std::invalid_argument("station-preserving LNS requires the paper VND profile");
  }
  if (preserve_station_visits_) repair_decision_memo_ = std::make_unique<RepairDecisionMemo>();
}

DestroyedSolution Perturbation::destroy(const Solution& solution,
                                        DestroyOperator operation,
                                        std::size_t removal_count,
                                        double penalty_lambda) {
  std::vector<int> customers = all_customers(solution);
  removal_count = std::min(removal_count, customers.size());
  if (removal_count == 0) return DestroyedSolution{solution, {}};
  std::vector<int> removed;
  removed.reserve(removal_count);

  // TRC 5.1.4 removes exactly kappa customers.  Route-oriented operators
  // keep their route selection. Only their last oversize route has an
  // explicit A subset choice; Closest keeps its ordered correlation prefix.
  const auto append_up_to_quota = [&](const std::vector<int>& route,
                                     bool route_oriented = false) {
    const std::size_t count =
        std::min(route.size(), removal_count - removed.size());
    if (count == 0) return;
    if (route_oriented && count < route.size()) {
      ++statistics_.partial_route_events;
      if (partial_route_subset_ == PartialRouteSubset::Random) {
        ++statistics_.randomized_partial_routes;
        std::vector<std::size_t> positions(route.size());
        std::iota(positions.begin(), positions.end(), 0);
        // Partial Fisher-Yates chooses a uniform subset without replacement.
        // Normalize its output order to avoid separately permuting repair input.
        for (std::size_t index = 0; index < count; ++index) {
          const std::size_t selected =
              std::uniform_int_distribution<std::size_t>(index, route.size() - 1)(random_);
          std::swap(positions[index], positions[selected]);
        }
        std::sort(positions.begin(), positions.begin() + static_cast<std::ptrdiff_t>(count));
        for (std::size_t index = 0; index < count; ++index) removed.push_back(route[positions[index]]);
        return;
      }
    }
    removed.insert(removed.end(), route.begin(),
                   route.begin() + static_cast<std::ptrdiff_t>(count));
  };

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
      append_up_to_quota(route_customers, true);
      if (removed.size() >= removal_count) break;
    }
  } else if (operation == DestroyOperator::ClosestRemoval) {
    const int seed = customers[random_index(customers, random_)];
    removed.push_back(seed);
    std::erase(customers, seed);
    std::sort(customers.begin(), customers.end(), [&](int left, int right) {
      return std::tuple{instance_.distance_km(seed, left), left} <
             std::tuple{instance_.distance_km(seed, right), right};
    });
    append_up_to_quota(customers);
  } else {
    struct TargetContribution {
      double value{};
      int customer_id{};
      std::size_t route_index{};
    };
    std::vector<TargetContribution> contribution;
    const auto reference_gaps = preserve_station_visits_
        ? PlanFactory::GapMap{} : factory_.gap_map(solution);
    if (!preserve_station_visits_ && sequence_scorer_ == nullptr) {
      sequence_scorer_ = std::make_unique<RouteSequenceScorer>(
          instance_, solution, reference_gaps);
    } else if (!preserve_station_visits_) {
      sequence_scorer_->reset(solution, reference_gaps);
    }
    for (std::size_t route_index = 0; route_index < solution.plans.size();
         ++route_index) {
      const Plan& plan = solution.plans[route_index];
      const double current_cost =
          generalized_cost(*plan.evaluation, instance_, penalty_lambda);
      std::optional<PhysicalRouteScorer> physical_scorer;
      if (preserve_station_visits_) physical_scorer.emplace(instance_, plan);
      for (const int customer_id : plan.customer_ids) {
        double reduced_cost = 0.0;
        if (plan.customer_ids.size() > 1 && physical_scorer.has_value()) {
          const auto position = std::find(plan.route.visits.begin(),
                                          plan.route.visits.end(), customer_id);
          reduced_cost = generalized_cost(physical_scorer->evaluate_removal(
              static_cast<std::size_t>(position - plan.route.visits.begin())),
              instance_, penalty_lambda);
        } else if (plan.customer_ids.size() > 1) {
          auto reduced = plan.customer_ids;
          std::erase(reduced, customer_id);
          reduced_cost = sequence_cost(reduced, *sequence_scorer_, penalty_lambda);
        }
        contribution.push_back(TargetContribution{
            .value = current_cost - reduced_cost,
            .customer_id = customer_id,
            .route_index = route_index});
      }
    }
    // Only the target is needed. Pairwise tolerance ties are not transitive
    // and cannot be used as a std::sort comparator. Anchor the tie band to
    // the actual maximum, then choose its smallest customer/route key.
    const auto maximum = std::max_element(
        contribution.begin(), contribution.end(),
        [](const auto& left, const auto& right) {
          return left.value < right.value;
        });
    auto selected = maximum;
    for (auto candidate = contribution.begin(); candidate != contribution.end();
         ++candidate) {
      if (maximum->value - candidate->value <= kCostTolerance &&
          std::tie(candidate->customer_id, candidate->route_index) <
              std::tie(selected->customer_id, selected->route_index)) {
        selected = candidate;
      }
    }
    const int target = selected->customer_id;
    const std::size_t target_route = selected->route_index;
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
      if (removed.size() == removal_count) break;
      const auto& route_customers =
          solution.plans[candidate.route_index].customer_ids;
      append_up_to_quota(route_customers, true);
    }
    // The target route supplies any remaining quota after other routes.
    if (removed.size() < removal_count) {
      std::vector<int> eligible;
      for (const int customer_id :
           solution.plans[target_route].customer_ids) {
        if (customer_id != target) eligible.push_back(customer_id);
      }
      append_up_to_quota(eligible, true);
    }
  }

  if (preserve_station_visits_) {
    return destroy_joint_occurrences(solution, std::move(removed));
  }
  auto sequences = remove_customers(solution, removed);
  Solution remaining =
      build_solution_from_sequences(sequences, factory_, &solution);
  return DestroyedSolution{.remaining = std::move(remaining),
                           .removed_customers = std::move(removed)};
}

// User-authorized experimental extension: destroy only samples customers;
// stations in the damaged customer blocks are subsequently movable tokens.
DestroyedSolution Perturbation::destroy_joint_occurrences(
    const Solution& solution, std::vector<int> removed_customers) {
  const std::unordered_set<int> removed(removed_customers.begin(),
                                        removed_customers.end());
  DestroyedSolution result;
  result.removed_customers = std::move(removed_customers);
  for (std::size_t route = 0; route < solution.plans.size(); ++route) {
    const auto& plan = solution.plans[route];
    if (std::none_of(plan.customer_ids.begin(), plan.customer_ids.end(),
                     [&](int id) { return removed.contains(id); })) {
      result.remaining.plans.push_back(plan);
      continue;
    }
    std::vector<int> visits;
    visits.reserve(plan.route.visits.size());
    std::size_t customer_gap = 0;
    for (std::size_t visit = 0; visit < plan.route.visits.size(); ++visit) {
      const int id = plan.route.visits[visit];
      if (instance_.node(id).kind == NodeKind::Customer) {
        if (!removed.contains(id)) visits.push_back(id);
        ++customer_gap;
      } else {
        // A gap is inside a damaged consecutive-customer block exactly when
        // either adjacent projected customer was removed. Depot boundaries
        // are retained. Union the blocks by visiting each raw occurrence once.
        const bool left = customer_gap != 0 &&
            removed.contains(plan.customer_ids[customer_gap - 1]);
        const bool right = customer_gap < plan.customer_ids.size() &&
            removed.contains(plan.customer_ids[customer_gap]);
        if (left || right) {
          result.removed_stations.push_back({id, route, visit});
        } else {
          visits.push_back(id);
        }
      }
    }
    // All CS of a customer-empty route were extracted before normalization.
    if (visits.empty()) continue;
    if (factory_.has_adjacent_station_visits(visits)) {
      throw std::logic_error("joint destruction created adjacent identical CS");
    }
    result.remaining.plans.push_back(factory_.plan_visits(std::move(visits)));
  }
  return result;
}

Solution Perturbation::repair_joint_occurrences(
    const Solution& original, DestroyedSolution destroyed,
    RepairOperator operation, double penalty_lambda) {
  // Do not use the old customer-only decision tape: it has neither occurrence
  // identity nor the additional station insertion steps.
  last_joint_failure_.clear();
  ++statistics_.joint_repair_calls;
  statistics_.joint_removed_station_tokens += destroyed.removed_stations.size();
  statistics_.joint_max_removed_station_tokens = std::max<std::uint64_t>(
      statistics_.joint_max_removed_station_tokens, destroyed.removed_stations.size());
  statistics_.joint_max_pending_tokens = std::max<std::uint64_t>(
      statistics_.joint_max_pending_tokens,
      destroyed.removed_customers.size() + destroyed.removed_stations.size());
  struct Token {
    int id;
    bool station;
    std::size_t source_route, source_visit;
    std::size_t ordinal{};
  };
  std::vector<Token> pending;
  for (int id : destroyed.removed_customers)
    pending.push_back({id, false, 0, 0});
  for (const auto& token : destroyed.removed_stations)
    pending.push_back({token.station_id, true, token.source_route,
                       token.source_visit});
  for (std::size_t i = 0; i < pending.size(); ++i) pending[i].ordinal = i;
  Solution current = std::move(destroyed.remaining);
  // Active partial indices are stable under insertion. They need not equal
  // original route indices after a customer-empty source route disappeared.
  std::vector<std::vector<int>> retained;
  for (const auto& p : current.plans) retained.push_back(p.route.visits);
  const auto fallback = [&](const char* reason) {
    last_joint_failure_ = reason;
    if (last_joint_failure_ == "all_pending_tokens_blocked")
      ++statistics_.joint_deadlock_rollbacks;
    else ++statistics_.joint_invariant_rollbacks;
    ++statistics_.station_conservation_checks;
    ++statistics_.station_conservation_rollbacks;
    ++statistics_.station_conservation_noops;
    return original;
  };
  if (!destroyed.station_shells.empty()) return fallback("unexpected_station_shell");
  std::vector<std::unique_ptr<PhysicalRouteScorer>> scorers;
  scorers.reserve(current.plans.size() + destroyed.removed_customers.size());
  for (const auto& p : current.plans) {
    scorers.push_back(std::make_unique<PhysicalRouteScorer>(instance_, p));
    ++statistics_.physical_route_profiles_built;
  }
  NonlinearSequenceEvaluator sequence(instance_);
  // Per-repair ranks only: an unchanged route has the same full certificate,
  // raw slots and lambda. No Plan/owner/key is retained by this workspace.
  struct RouteRanks {
    std::size_t version{};
    std::size_t finite_count{};
    std::size_t size{};
    std::array<InsertionMove, 2> best;
  };
  const std::size_t route_capacity = current.plans.size() +
                                     destroyed.removed_customers.size();
  std::vector<RouteRanks> route_ranks;
  if (route_capacity && pending.size() > route_ranks.max_size() / route_capacity)
    throw std::length_error("joint insertion rank workspace is too large");
  route_ranks.resize(pending.size() * route_capacity);
  std::vector<std::size_t> route_versions(route_capacity, 0);
  for (std::size_t r = 0; r < current.plans.size(); ++r) route_versions[r] = 1;
  while (!pending.empty()) {
    ++statistics_.physical_insertion_decisions;
    std::optional<InsertionMove> chosen;
    std::size_t chosen_token = 0;
    double best_regret = -1.0;
    for (std::size_t token_index = 0; token_index < pending.size(); ++token_index) {
      const auto& token = pending[token_index];
      std::vector<InsertionMove> moves;
      moves.reserve(current.plans.size() * 2 + !token.station);
      std::size_t logical_candidates = 0;
      const auto make_move = [&](std::size_t route, std::size_t slot,
                                 double incumbent, double cost) {
        const double delta = cost - incumbent;
        return InsertionMove{
            .customer_id = token.id, .route_index = route, .position = slot,
            .incumbent_route_cost = incumbent, .direct_delta = delta,
            .linked_lower_bound_delta = delta, .delta = delta};
      };
      for (std::size_t route = 0; route < current.plans.size(); ++route) {
        auto& ranks = route_ranks[token.ordinal * route_capacity + route];
        if (ranks.version != route_versions[route]) {
          ++statistics_.joint_rank_cache_misses;
          ranks.finite_count = ranks.size = 0;
          const auto append = [&](std::size_t slot, double incumbent, double cost) {
            auto move = make_move(route, slot, incumbent, cost);
            if (!std::isfinite(move.delta)) return;
            ++ranks.finite_count;
            ++statistics_.joint_rank_fresh_descriptors;
            // Two distinct original positions survive even at exactly equal
            // cost. The strict original tuple also keeps their slot order.
            if (!ranks.size || insertion_better(move, ranks.best[0])) {
              if (ranks.size) ranks.best[1] = ranks.best[0];
              ranks.best[0] = std::move(move);
              ranks.size = std::min<std::size_t>(2, ranks.size + 1);
            } else if (ranks.size == 1 || insertion_better(move, ranks.best[1])) {
              ranks.best[1] = std::move(move);
              ranks.size = 2;
            }
          };
          const auto& plan = current.plans[route];
          const double incumbent = generalized_cost(*plan.evaluation, instance_,
                                                    penalty_lambda);
          std::optional<RouteCostSummary> summary;
          std::size_t representative = plan.route.visits.size() + 1;
          for (std::size_t slot = 0; slot <= plan.route.visits.size(); ++slot) {
            if (token.station) {
              if ((slot && plan.route.visits[slot - 1] == token.id) ||
                  (slot < plan.route.visits.size() &&
                   plan.route.visits[slot] == token.id)) continue;
              auto visits = plan.route.visits;
              visits.insert(visits.begin() + static_cast<std::ptrdiff_t>(slot),
                            token.id);
              ++statistics_.joint_rank_score_calls;
              append(slot, incumbent, sequence.generalized_cost(
                  sequence.route_label(Route{std::move(visits)}), penalty_lambda));
            } else {
              const auto rep = scorers[route]->insertion_score_representative(slot);
              if (!summary || rep != representative) {
                ++statistics_.joint_rank_score_calls;
                summary = scorers[route]->evaluate_insertion(slot, token.id);
                representative = rep;
              }
              append(slot, incumbent,
                     generalized_cost(*summary, instance_, penalty_lambda));
            }
          }
          ranks.version = route_versions[route];
        } else {
          ++statistics_.joint_rank_cache_hits;
          statistics_.joint_rank_reused_descriptors += ranks.finite_count;
        }
        logical_candidates += ranks.finite_count;
        for (std::size_t i = 0; i < ranks.size; ++i) moves.push_back(ranks.best[i]);
      }
      if (!token.station) {
        // Keep singleton scoring and its changing append-route index exactly
        // where they were; only existing unchanged routes share ranks.
        ++statistics_.joint_rank_score_calls;
        auto singleton = make_move(current.plans.size(), 0, 0.0,
            sequence.generalized_cost(sequence.route_label(Route{{token.id}}),
                                       penalty_lambda));
        if (std::isfinite(singleton.delta)) {
          ++logical_candidates;
          ++statistics_.joint_rank_fresh_descriptors;
          moves.push_back(std::move(singleton));
        }
      }
      statistics_.physical_insertion_candidates += logical_candidates;
      if (moves.empty()) {
        // A later customer insertion can create a separator for this station.
        ++joint_blocked_token_scans_;
        ++statistics_.joint_blocked_token_scans;
        continue;
      }
      if (operation == RepairOperator::GreedyInsertion) {
        const auto best = std::min_element(moves.begin(), moves.end(), insertion_better);
        if (!chosen || insertion_better(*best, *chosen)) {
          chosen = *best;
          chosen_token = token_index;
        }
      } else {
        const auto n = std::min<std::size_t>(2, moves.size());
        std::partial_sort(moves.begin(), moves.begin() +
            static_cast<std::ptrdiff_t>(n), moves.end(), insertion_better);
        const double regret = n == 1 ? std::numeric_limits<double>::infinity()
                                     : moves[1].delta - moves[0].delta;
        if (!chosen || regret_better(regret, moves[0], best_regret, *chosen)) {
          chosen = moves[0]; chosen_token = token_index; best_regret = regret;
        }
      }
    }
    if (!chosen) return fallback("all_pending_tokens_blocked");
    const auto started = std::chrono::steady_clock::now();
    const auto token = pending[chosen_token];
    const auto route = chosen->route_index;
    std::vector<int> visits;
    if (route < current.plans.size()) visits = current.plans[route].route.visits;
    else if (token.station) return fallback("station_singleton_forbidden");
    visits.insert(visits.begin() + static_cast<std::ptrdiff_t>(chosen->position), token.id);
    if (factory_.has_adjacent_station_visits(visits)) return fallback("adjacent_station_candidate");
    apply_insertion(current, route, factory_.plan_visits(std::move(visits)));
    ++route_versions.at(route);
    if (token.station) {
      ++joint_station_insertions_;
      ++statistics_.joint_station_insertion_attempts;
    } else {
      ++joint_customer_insertions_;
      ++statistics_.joint_customer_insertion_attempts;
    }
    // Erase exactly one occurrence, even if several tokens share a station ID.
    pending.erase(pending.begin() + static_cast<std::ptrdiff_t>(chosen_token));
    if (!pending.empty()) {
      auto updated = std::make_unique<PhysicalRouteScorer>(instance_, current.plans[route]);
      if (route == scorers.size()) scorers.push_back(std::move(updated));
      else scorers[route] = std::move(updated);
      ++statistics_.physical_route_profiles_built;
    }
    statistics_.repair_timing.winner_materialization_s += elapsed_seconds(started);
  }
  const auto counts = [&](const Solution& s) {
    std::map<int, std::size_t> result;
    for (const auto& p : s.plans)
      for (int id : p.route.visits) ++result[id];
    return result;
  };
  if (counts(original) != counts(current)) return fallback("occurrence_conservation");
  for (std::size_t route = 0; route < current.plans.size(); ++route) {
    const auto& p = current.plans[route];
    if (p.customer_ids.empty() || factory_.has_adjacent_station_visits(p.route.visits))
      return fallback("invalid_completed_route");
    if (route < retained.size()) {
      std::size_t i = 0;
      for (int id : p.route.visits)
        if (i < retained[route].size() && id == retained[route][i]) ++i;
      if (i != retained[route].size()) return fallback("retained_subsequence_changed");
    }
  }
  current.validate_partition(instance_);
  ++statistics_.station_conservation_checks;
  bool unchanged = current.plans.size() == original.plans.size();
  for (std::size_t r = 0; unchanged && r < current.plans.size(); ++r)
    unchanged = current.plans[r].route.visits == original.plans[r].route.visits;
  if (unchanged) ++statistics_.station_conservation_noops;
  else ++statistics_.station_conservation_effective;
  return current;
}

Solution Perturbation::repair(DestroyedSolution destroyed,
                              RepairOperator operation,
                              double penalty_lambda,
                              std::vector<InsertionPathHint>* hints,
                              RepairDecisionMemo* memo) {
  if (preserve_station_visits_) {
    return repair_preserving_stations(std::move(destroyed), operation,
                                      penalty_lambda, memo);
  }
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
            &path_sampler_, false, granular_.get(),
            &statistics_.candidates_rejected_by_granular,
            granular_ != nullptr ? &rejected : nullptr);
        candidates.insert(candidates.end(),
                          std::make_move_iterator(moves.begin()),
                          std::make_move_iterator(moves.end()));
      }
      if (granular_ != nullptr) {
        append_granular_escape(
            candidates, rejected, current, instance_, scorer,
            penalty_lambda, &path_sampler_, random_,
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
            true, &path_sampler_, false, granular_.get(),
            &statistics_.candidates_rejected_by_granular);
        candidates.insert(candidates.end(),
                          std::make_move_iterator(moves.begin()),
                          std::make_move_iterator(moves.end()));
      }
      chosen = choose_linked_greedy(
          std::move(candidates), current, instance_, factory_, path_sampler_,
          reference_gaps, penalty_lambda, statistics_, completion_memo,
          single_step_repair_completion_);
    } else if (completes_insertion_candidates(profile_)) {
      chosen = choose_linked_two_regret(
          current, unassigned, instance_, factory_, path_sampler_, scorer,
          reference_gaps, penalty_lambda, &path_sampler_,
          statistics_, completion_memo, single_step_repair_completion_);
    } else if (operation == RepairOperator::GreedyInsertion) {
      for (const int customer_id : unassigned) {
        auto moves = insertion_moves(
            current, customer_id, instance_, scorer, penalty_lambda,
            false, nullptr, false,
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
                                     penalty_lambda, false, nullptr, false,
                                     granular_.get(),
                                     &statistics_.candidates_rejected_by_granular);
        std::sort(moves.begin(), moves.end(), insertion_better);
        const double regret = moves.size() == 1
                                  ? std::numeric_limits<double>::infinity()
                                  : moves[1].delta - moves[0].delta;
        if (!chosen.has_value() ||
            regret_better(regret, moves[0], best_regret, *chosen)) {
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

Solution Perturbation::repair_preserving_stations(
    DestroyedSolution destroyed, RepairOperator operation,
    double penalty_lambda, RepairDecisionMemo* memo) {
  const auto replay = memo ? memo->find(destroyed, operation, statistics_) : std::nullopt;
  std::optional<RepairDecisionMemo::Entry> recording;
  if (memo && memo->active && !replay) recording = memo->input(destroyed, operation);
  const auto candidates_before = statistics_.physical_insertion_candidates;
  const auto profiles_before = statistics_.physical_route_profiles_built;
  std::size_t replay_index = 0;
  Solution current = std::move(destroyed.remaining);
  std::vector<int> unassigned = std::move(destroyed.removed_customers);
  auto station_shells = std::move(destroyed.station_shells);
  std::vector<std::unique_ptr<PhysicalRouteScorer>> scorers;
  scorers.reserve(current.plans.size() + unassigned.size());
  if (!replay) for (const Plan& plan : current.plans) {
    scorers.push_back(std::make_unique<PhysicalRouteScorer>(instance_, plan));
    ++statistics_.physical_route_profiles_built;
  }
  NonlinearSequenceEvaluator singleton_scorer(instance_);
  while (!unassigned.empty()) {
    // Each nonempty dormant shell must receive a customer before the final
    // station-conservation guard. Once customers and shells are equally
    // numerous, any nonactivation would make that completion impossible.
    // Shells remain inactive at zero cost; activation still pays the full route.
    const std::size_t pending_shells = std::count_if(
        station_shells.begin(), station_shells.end(),
        [](const auto& shell) { return !shell.empty(); });
    const bool must_activate_shell = unassigned.size() == pending_shells;
    ++statistics_.physical_insertion_decisions;
    std::optional<InsertionMove> chosen;
    double best_regret = -1.0;
    if (replay) {
      const auto& step = replay->steps[replay_index++];
      chosen = InsertionMove{.customer_id = step.customer,
                             .route_index = step.route, .position = step.slot};
    } else for (const int customer : unassigned) {
      std::vector<InsertionMove> moves;
      for (std::size_t route = 0;
           !must_activate_shell && route < current.plans.size(); ++route) {
        const Plan& plan = current.plans[route];
        const double incumbent = generalized_cost(
            *plan.evaluation, instance_, penalty_lambda);
        // A score is retained only within this fixed route/customer scan.
        // Keep every descriptor, so Greedy ties and TwoRegret multiplicity
        // are exactly those of the original full physical-slot enumeration.
        std::optional<RouteCostSummary> insertion_summary;
        std::size_t summary_representative = plan.route.visits.size()+1;
        for (std::size_t slot = 0; slot <= plan.route.visits.size(); ++slot) {
          const std::size_t representative =
              scorers[route]->insertion_score_representative(slot);
          if (!insertion_summary || representative != summary_representative) {
            insertion_summary = scorers[route]->evaluate_insertion(slot, customer);
            summary_representative = representative;
          }
          const double delta = generalized_cost(
              *insertion_summary, instance_, penalty_lambda) - incumbent;
          moves.push_back(InsertionMove{.customer_id = customer,
                                        .route_index = route,
                                        .position = slot,
                                        .incumbent_route_cost = incumbent,
                                        .direct_delta = delta,
                                        .linked_lower_bound_delta = delta,
                                        .delta = delta});
        }
      }
      // A dormant shell is charged only after its first customer activates it.
      // Score the full activated route, never a delta from a phantom CS route.
      for (std::size_t shell = 0; shell < station_shells.size(); ++shell) {
        if (station_shells[shell].empty()) continue;
        std::unordered_set<std::vector<int>, IntVectorHash> distinct_activations;
        for (std::size_t slot = 0; slot <= station_shells[shell].size(); ++slot) {
          auto visits = station_shells[shell];
          visits.insert(visits.begin() + static_cast<std::ptrdiff_t>(slot), customer);
          // Preserve every inherited CS occurrence, matching destruction and
          // the exact conservation guard. Deduplicate only identical vectors.
          if (!distinct_activations.insert(visits).second) continue;
          const double cost = singleton_scorer.generalized_cost(
              singleton_scorer.route_label(Route{std::move(visits)}), penalty_lambda);
          moves.push_back(InsertionMove{.customer_id = customer,
              .route_index = current.plans.size() + shell, .position = slot,
              .incumbent_route_cost = 0.0, .direct_delta = cost,
              .linked_lower_bound_delta = cost, .delta = cost});
        }
      }
      if (!must_activate_shell) {
        const double singleton_cost = singleton_scorer.generalized_cost(
            singleton_scorer.route_label(Route{{customer}}), penalty_lambda);
        moves.push_back(InsertionMove{.customer_id = customer,
                                      .route_index = current.plans.size() + station_shells.size(),
                                      .position = 0,
                                      .incumbent_route_cost = 0.0,
                                      .direct_delta = singleton_cost,
                                      .linked_lower_bound_delta = singleton_cost,
                                      .delta = singleton_cost});
      }
      statistics_.physical_insertion_candidates += moves.size();
      if (operation == RepairOperator::GreedyInsertion) {
        const auto best = std::min_element(moves.begin(), moves.end(),
                                           insertion_better);
        if (!chosen.has_value() || insertion_better(*best, *chosen)) chosen = *best;
      } else {
        const auto count = std::min<std::size_t>(2, moves.size());
        std::partial_sort(moves.begin(), moves.begin() +
            static_cast<std::ptrdiff_t>(count), moves.end(), insertion_better);
        const double regret = count == 1 ? std::numeric_limits<double>::infinity()
                                        : moves[1].delta - moves[0].delta;
        if (!chosen.has_value() ||
            regret_better(regret, moves[0], best_regret, *chosen)) {
          chosen = moves[0];
          best_regret = regret;
        }
      }
    }
    if (recording) recording->steps.push_back(
        {chosen->customer_id, chosen->route_index, chosen->position});
    const auto started = std::chrono::steady_clock::now();
    std::size_t route = chosen->route_index;
    const int customer = chosen->customer_id;
    std::vector<int> visits;
    if (route < current.plans.size()) {
      visits = current.plans[route].route.visits;
    } else {
      const std::size_t shell = route - current.plans.size();
      if (shell < station_shells.size()) {
        visits = std::move(station_shells[shell]);
        station_shells[shell].clear();
        ++statistics_.station_shells_activated;
      }
      route = current.plans.size();
    }
    visits.insert(visits.begin() + static_cast<std::ptrdiff_t>(chosen->position),
                  customer);
    apply_insertion(current, route, factory_.plan_visits(std::move(visits)));
    std::erase(unassigned, customer);
    if (!replay && !unassigned.empty()) {
      auto updated = std::make_unique<PhysicalRouteScorer>(instance_, current.plans[route]);
      if (route == scorers.size()) scorers.push_back(std::move(updated));
      else scorers[route] = std::move(updated);
      ++statistics_.physical_route_profiles_built;
    }
    statistics_.repair_timing.winner_materialization_s += elapsed_seconds(started);
  }
  statistics_.station_shells_discarded += std::count_if(
      station_shells.begin(), station_shells.end(),
      [](const auto& shell) { return !shell.empty(); });
  current.validate_partition(instance_);
  if (recording) {
    recording->candidates = statistics_.physical_insertion_candidates - candidates_before;
    recording->profiles = statistics_.physical_route_profiles_built - profiles_before;
    memo->put(std::move(*recording), statistics_);
  }
  return current;
}

Solution Perturbation::enforce_station_conservation(
    const Solution& original, Solution repaired,
    const std::vector<std::vector<int>>& expected_station_chains) {
  ++statistics_.station_conservation_checks;
  // Route indices can change when dormant shells reactivate. Compare a
  // multiset of full ordered CS sequences, with exact vector equality after
  // hashing; station-free routes do not carry a charging skeleton.
  const auto skeletons = [&](const Solution& input) {
    std::unordered_map<std::vector<int>, std::size_t, IntVectorHash> counts;
    counts.reserve(input.plans.size());
    for (const Plan& plan : input.plans) {
      std::vector<int> stations;
      for (const int id : plan.route.visits) {
        if (instance_.node(id).kind == NodeKind::Station) stations.push_back(id);
      }
      if (!stations.empty()) ++counts[std::move(stations)];
    }
    return counts;
  };
  std::unordered_map<std::vector<int>, std::size_t, IntVectorHash> expected;
  for(const auto& chain:expected_station_chains)if(!chain.empty())++expected[chain];
  if (expected != skeletons(repaired)) {
    // No resampling, CS replacement, phantom route or alternative pricing.
    // The caller's original fully evaluated solution is the bounded fallback.
    ++statistics_.station_conservation_rollbacks;
    ++statistics_.station_conservation_noops;
    return original;
  }
  bool unchanged = original.plans.size() == repaired.plans.size();
  for (std::size_t route = 0; unchanged && route < original.plans.size(); ++route) {
    unchanged = original.plans[route].route.visits == repaired.plans[route].route.visits;
  }
  if (unchanged) ++statistics_.station_conservation_noops;
  else ++statistics_.station_conservation_effective;
  return repaired;
}

Solution Perturbation::apply(const Solution& solution, double penalty_lambda) {
  if (!preserve_station_visits_ && repair_decision_memo_)
    repair_decision_memo_->begin(solution, penalty_lambda, statistics_);
  last_removal_count_ = 0;
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
  last_removal_count_ = destroyed.removed_customers.size();
  statistics_.repair_timing.destroy_s += elapsed_seconds(destroy_started);
  if (preserve_station_visits_) {
    return repair_joint_occurrences(solution, std::move(destroyed),
                                    repair_operation, penalty_lambda);
  }
  Solution repaired = repair(
      std::move(destroyed), repair_operation, penalty_lambda,
      uses_granular_topology(profile_) ? &hints : nullptr, repair_decision_memo_.get());
  if (!uses_granular_topology(profile_) || hints.empty()) return repaired;
  const auto post_hint_started = std::chrono::steady_clock::now();
  HintedPathResult hinted = path_sampler_.apply_insertion_hints(
      std::move(repaired), hints, penalty_lambda);
  statistics_.repair_timing.post_hint_s += elapsed_seconds(post_hint_started);
  statistics_.post_repair_hint_attempts += hinted.attempts;
  statistics_.post_repair_hint_accepted += hinted.accepted;
  return std::move(hinted.solution);
}

double VariableNeighborhoodDescent::PhysicalScoreSummary::generalized_cost(
    const Instance& instance, double penalty_lambda) const {
  if (!std::isfinite(penalty_lambda) || penalty_lambda < 0.0) {
    throw std::invalid_argument("penalty lambda must be finite and nonnegative");
  }
  // Identical operations/ordering to generalized_cost(RouteCostSummary).
  return raw_cost_h +
         penalty_lambda *
             (energy_shortfall_wh /
                  instance.vehicle().consumption_wh_per_km +
              duration_excess_h);
}

VariableNeighborhoodDescent::VariableNeighborhoodDescent(
    const Instance& instance, PlanFactory& factory, PathSampler& path_sampler,
    std::mt19937_64& random, std::size_t promising_arc_count,
    SearchProfile profile, std::size_t replace_path_attempts,
    bool project_opening_station_visits, bool preserve_swap_station_visits)
    : instance_(instance),
      factory_(factory),
      granular_(instance, promising_arc_count),
      path_sampler_(path_sampler),
      random_(random),
      profile_(profile),
      replace_path_attempts_(replace_path_attempts),
      project_opening_station_visits_(project_opening_station_visits),
      preserve_swap_station_visits_(preserve_swap_station_visits) {
  if (preserve_swap_station_visits_ && profile_ != SearchProfile::Paper) {
    throw std::invalid_argument("physical-slot Swap requires the paper profile");
  }
  if (replace_path_attempts_ != 1) {
    throw std::invalid_argument("only the single-proposal ReplacePath stage is active; the four-trial experiment was rejected");
  }
  if (project_opening_station_visits_ && profile_ != SearchProfile::Paper) {
    throw std::invalid_argument("projected route opening requires the paper profile");
  }
}

bool VariableNeighborhoodDescent::expired(
    std::optional<std::chrono::steady_clock::time_point> deadline) const {
  return deadline.has_value() && std::chrono::steady_clock::now() >= *deadline;
}

Solution VariableNeighborhoodDescent::improve(
    Solution solution, double penalty_lambda,
    std::optional<std::chrono::steady_clock::time_point> deadline,
    std::optional<std::size_t> cooperative_removal_count) {
  if (cooperative_removal_count.has_value() &&
      (profile_ != SearchProfile::Paper || path_sampler_.implicit_replacement_)) {
    throw std::invalid_argument(
        "cooperative RP requires the Paper profile and explicit path sampling");
  }
  // Reuse only complete customer-support winners at the identical lambda.
  // Absolute opening slots also require the same global route count.
  const bool reuse_scope = cooperative_removal_count.has_value() &&
      profile_ == SearchProfile::Paper && supported_bound_fp_environment();
  if (!reuse_scope || !move_memory_reusable_ ||
      !move_memory_lambda_.has_value() ||
      *move_memory_lambda_ != penalty_lambda ||
      move_memory_route_count_ != solution.plans.size()) {
    move_memory_.clear();
    route_versions_.clear();
    next_route_version_ = 1;
  }
  move_memory_usable_ = true;
  move_memory_reusable_ = reuse_scope;
  move_memory_lambda_ = penalty_lambda;
  move_memory_route_count_ = solution.plans.size();

  // All schedule state is local to this VND. Accepted moves invalidate only
  // the current layout's selector; they never replenish the proposal budget.
  constexpr std::size_t no_deficit_index =
      std::numeric_limits<std::size_t>::max();
  struct CooperativePathState {
    std::size_t remaining{};
    detail::RouteOpeningRpReserve opening_reserve;
    std::vector<std::pair<std::size_t, std::size_t>> gaps;
    std::size_t next_gap{};
    bool selector_valid{};
    // Positions in gaps, plus their inverse indices, for O(1) removal.
    std::vector<std::size_t> deficit_positions;
    std::vector<std::size_t> deficit_index;
  };
  std::optional<CooperativePathState> cooperative;
  std::size_t entry_allowance = 0;
  if (cooperative_removal_count.has_value()) {
    cooperative.emplace();
    std::size_t entry_gap_count = 0;
    for (const Plan& plan : solution.plans) entry_gap_count += plan.gaps.size();
    entry_allowance = std::min(*cooperative_removal_count, entry_gap_count);
    // min(G_entry, 2*kappa), without overflowing the doubling operation.
    cooperative->remaining = entry_allowance +
        std::min(*cooperative_removal_count, entry_gap_count - entry_allowance);
    cooperative->opening_reserve.initialize(cooperative->remaining);
    entry_allowance = std::min(entry_allowance,
        cooperative->opening_reserve.ordinary_allowance(cooperative->remaining));
    ++statistics_.cooperative_vnd_calls;
    statistics_.cooperative_proposal_budget += cooperative->remaining;
  }

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
      move_memory_route_count_ = solution.plans.size();
    }
    if (cooperative.has_value()) {
      cooperative->gaps.clear();
      cooperative->deficit_positions.clear();
      cooperative->deficit_index.clear();
      cooperative->next_gap = 0;
      cooperative->selector_valid = false;
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

  const auto cooperative_neighbor = [&](std::size_t allowance, bool entry_phase,
                                        bool opening_guidance = false)
      -> std::optional<Solution> {
    ++statistics_.replace_path_calls;
    auto& state = *cooperative;
    if (state.remaining == 0 || allowance == 0 || expired(deadline)) {
      return std::nullopt;
    }
#if ILS_SP_TYPED_DELETE_ENTRY
    std::vector<std::size_t> deletion_positions;
#endif
    if (!state.selector_valid) {
      for (std::size_t route = 0; route < solution.plans.size(); ++route) {
        const Plan& plan = solution.plans[route];
        const bool energy_deficient =
            plan.evaluation->energy_shortfall_wh > kEnergyToleranceWh;
        for (std::size_t gap = 0; gap < plan.gaps.size(); ++gap) {
          const std::size_t position = state.gaps.size();
          state.gaps.emplace_back(route, gap);  // Includes direct gaps.
#if ILS_SP_TYPED_DELETE_ENTRY
          // Every explicit catalogue includes its direct path unconditionally.
          // Build once at VND entry, with no catalogue lookup or extra scoring.
          if (entry_phase && !plan.gaps[gap].empty()) {
            deletion_positions.push_back(position);
          }
#endif
          state.deficit_index.push_back(energy_deficient
              ? state.deficit_positions.size() : no_deficit_index);
          if (energy_deficient) state.deficit_positions.push_back(position);
        }
      }
      state.selector_valid = true;
    }
    const double current_cost =
        generalized_cost(solution, instance_, penalty_lambda);
    const std::size_t attempts = std::min(allowance, state.remaining);
    for (std::size_t attempt = 0;
         attempt < attempts && state.next_gap < state.gaps.size(); ++attempt) {
      if (expired(deadline)) break;
      std::optional<std::mt19937_64> before_catalogue_read;
      if (path_sampler_.catalog_disk_cache_.enabled) {
        before_catalogue_read.emplace(random_);
      }
      // Only the first entry proposal and opening-owned token may use the proxy.
      // The uniform half reaches the original draw without any catalogue scan.
      std::optional<PathSampler::BoundedProxyChoice> guided;
      // Energy-repair guidance is inapplicable to fully feasible or time-only
      // states. Bypass before the coin so those proposals keep exact f0 RNG.
      if ((opening_guidance || entry_phase) && attempt == 0 &&
          std::any_of(solution.plans.begin(), solution.plans.end(),
              [](const Plan& p) {
                return p.evaluation->energy_shortfall_wh > kEnergyToleranceWh;
              })) {
        const auto before_hint = random_;
        try {
          if (std::bernoulli_distribution(0.5)(random_)) {
            auto scan = path_sampler_.bounded_proxy_choice(
                solution, state.gaps, state.next_gap, penalty_lambda);
            guided = std::move(scan.choice);
            // No partial/capped scan may select a proposal. Restore the coin
            // on unsupported scope so fallback is the literal old RNG draw.
            if (!guided) random_ = before_hint;
          }
        } catch (...) {
          random_ = before_hint;
          throw;
        }
      }
      // Uniform fallback: one uniformly random unvisited gap,
      // followed by one uniformly random noncurrent catalogue path. The typed
      // direct-deletion entry is disabled by this candidate's build contract.
      // BEGIN_DEFICIT_GAP_DRAW
#if ILS_SP_TYPED_DELETE_ENTRY
      const bool typed_deletion = entry_phase && attempt == 0 &&
          !deletion_positions.empty();
      std::size_t selected;
      if (typed_deletion) {
        selected = deletion_positions[std::uniform_int_distribution<std::size_t>(
            0, deletion_positions.size() - 1)(random_)];
      } else {
#endif
#if ILS_SP_TYPED_DELETE_ENTRY
      selected = guided ? guided->position : std::uniform_int_distribution<std::size_t>(
#else
      const std::size_t selected = guided ? guided->position : std::uniform_int_distribution<std::size_t>(
#endif
          state.next_gap, state.gaps.size() - 1)(random_);
#if ILS_SP_TYPED_DELETE_ENTRY
      }
#endif
#if ILS_SP_TYPED_DELETE_ENTRY
      const auto selected_gap = state.gaps[selected];
      // Typed proposals have their own one-shot allowance, not an ordinary
      // gap visit. A rejected deletion leaves every ordinary path reachable.
      if (!typed_deletion) {
#endif
      const std::size_t first = state.next_gap;
      const std::size_t first_deficit = state.deficit_index[first];
      const std::size_t selected_deficit = state.deficit_index[selected];
      std::swap(state.gaps[first], state.gaps[selected]);
      std::swap(state.deficit_index[first], state.deficit_index[selected]);
      if (first_deficit != no_deficit_index) {
        state.deficit_positions[first_deficit] = selected;
      }
      if (selected_deficit != no_deficit_index) {
        state.deficit_positions[selected_deficit] = first;
        const std::size_t last_position = state.deficit_positions.back();
        state.deficit_positions[selected_deficit] = last_position;
        state.deficit_index[last_position] = selected_deficit;
        state.deficit_positions.pop_back();
      }
      state.deficit_index[first] = no_deficit_index;
#if ILS_SP_TYPED_DELETE_ENTRY
      ++state.next_gap;
      }
      const auto [route, gap] = selected_gap;
#else
      const auto [route, gap] = state.gaps[state.next_gap++];
#endif
      // END_DEFICIT_GAP_DRAW
      --state.remaining;
      ++statistics_.replace_path_proposals;
      if (entry_phase) ++statistics_.cooperative_entry_proposals;
      const Plan& plan = solution.plans[route];
      const int origin = gap == 0 ? instance_.depot().id
                                  : plan.customer_ids[gap - 1];
      const int target = gap == plan.customer_ids.size() ? instance_.depot().id
                                                        : plan.customer_ids[gap];
      std::optional<std::vector<int>> sampled;
      try {
        // Guided winner is an existing noncurrent catalogue path. Otherwise
        // use the exact indexed uniform sampler. Only this path is full-scored.
#if ILS_SP_TYPED_DELETE_ENTRY
        if (typed_deletion) {
          sampled.emplace();  // The catalogue direct path, not a charging repair.
        } else
#endif
        {
        if (guided) sampled = std::move(guided->stations);
        else sampled = path_sampler_.sample_explicit_path(origin, target, plan.gaps[gap]);
        }
      } catch (...) {
        // Match legacy RP's rollback of the gap draw on a failed disk read.
        if (before_catalogue_read) random_ = *before_catalogue_read;
        throw;
      }
      if (!sampled.has_value()) {
        ++statistics_.replace_path_null_proposals;
        continue;
      }
      if (expired(deadline)) break;
      auto gaps = plan.gaps;
      gaps[gap] = std::move(*sampled);
      Plan proposed = factory_.make_plan(plan.customer_ids, std::move(gaps));
      if (expired(deadline)) break;
      // Preserve generalized_cost(Solution)'s route addition order; no delta
      // or subtract/add update, and no change to this VND's penalty lambda.
      double candidate_cost = 0.0;
      for (std::size_t index = 0; index < solution.plans.size(); ++index) {
        const Plan& candidate_plan = index == route ? proposed : solution.plans[index];
        candidate_cost += generalized_cost(
            *candidate_plan.evaluation, instance_, penalty_lambda);
      }
      if (candidate_cost < current_cost - kCostTolerance) {
        if (expired(deadline)) break;
        ++statistics_.replace_path_accepted;
        Solution candidate = solution;
        candidate.plans[route] = std::move(proposed);
        return candidate;  // First improvement, never best-of-budget.
      }
      ++statistics_.replace_path_nonimproving_proposals;
    }
    return std::nullopt;
  };

  // Entry is attempted exactly once. An accepted RP returns to Swap just as
  // a normal tail RP does; an unused entry allowance stays in the total budget.
  if (cooperative.has_value() && !expired(deadline)) {
    auto entry = cooperative_neighbor(entry_allowance, true);
    if (!expired(deadline)) accept(std::move(entry));
  }
  // Experimental remaining-list RVND for cooperative calls. Legacy calls
  // still choose the first remaining bit, preserving their fixed-order RNG.
  const auto active_neighborhoods = [&]() {
    unsigned int active = (1U << order_.size()) - 1U;
    if (cooperative.has_value() &&
        cooperative->opening_reserve.ordinary_allowance(cooperative->remaining) == 0) {
      active &= ~(1U << 3);  // Protect the one event-or-stagnation token.
    }
    return active;
  };
  unsigned int remaining = active_neighborhoods();
  while ((remaining != 0 ||
          (cooperative.has_value() && cooperative->opening_reserve.held)) &&
         !expired(deadline)) {
    if (remaining == 0) {
      // No route-opening event used the reserved token. Release it once after
      // normal stagnation; failure ends, strict acceptance restarts the same
      // invocation without repeating entry or replenishing any token.
      cooperative->opening_reserve.release();
      auto tail = cooperative_neighbor(1, false);
      if (expired(deadline) || !accept(std::move(tail))) break;
      remaining = active_neighborhoods();
      ++statistics_.rvnd_improvement_resets;
      continue;
    }
    const std::size_t neighborhood_index = cooperative.has_value()
        ? select_active_neighborhood(remaining)
        : static_cast<std::size_t>(std::countr_zero(remaining));
    if (cooperative.has_value()) ++statistics_.rvnd_neighborhood_selections;
    const std::size_t previous_route_count = solution.plans.size();
    auto best = cooperative.has_value() &&
                        order_[neighborhood_index] == Neighborhood::ReplacePath
        ? cooperative_neighbor(cooperative->opening_reserve.ordinary_allowance(
              cooperative->remaining), false)
        : best_neighbor(solution, order_[neighborhood_index],
                        penalty_lambda, deadline, cooperative.has_value());
    if (expired(deadline)) break;
    if (accept(std::move(best))) {
      remaining = active_neighborhoods();
      if (cooperative.has_value()) ++statistics_.rvnd_improvement_resets;
      if (cooperative.has_value() &&
          order_[neighborhood_index] != Neighborhood::ReplacePath &&
          cooperative->opening_reserve.release_on_opening(
              previous_route_count, solution.plans.size())) {
        // One proxy-guided or uniform RP proposal after the first strictly
        // accepted customer route opening. This consumes an existing token;
        // RP itself cannot trigger another opening event or refill the holder.
        auto opening = cooperative_neighbor(1, false, true);
        if (expired(deadline)) break;
        if (accept(std::move(opening))) {
          ++statistics_.rvnd_improvement_resets;
        }
        // A failed reserved probe retires its gap, not the ordinary RP bit.
        // The accepted customer move already reopened every neighborhood.
        remaining = active_neighborhoods();
      }
    } else {
      remaining &= ~(1U << neighborhood_index);
    }
  }
  solution.validate_partition(instance_);
  return solution;
}

std::size_t VariableNeighborhoodDescent::select_active_neighborhood(
    unsigned int active) {
  const unsigned int allowed = (1U << order_.size()) - 1U;
  if (active == 0 || (active & ~allowed) != 0) {
    throw std::invalid_argument("RVND active set is empty or invalid");
  }
  const unsigned int count = static_cast<unsigned int>(std::popcount(active));
  unsigned int ordinal = count == 1 ? 0
      : std::uniform_int_distribution<unsigned int>(0, count - 1)(random_);
  for (std::size_t index = 0; index < order_.size(); ++index) {
    if ((active & (1U << index)) != 0 && ordinal-- == 0) return index;
  }
  throw std::logic_error("RVND active ordinal is out of range");
}

std::optional<Solution> VariableNeighborhoodDescent::best_neighbor(
    const Solution& solution, Neighborhood neighborhood, double penalty_lambda,
    std::optional<std::chrono::steady_clock::time_point> deadline,
    bool cooperative_physical_customers) {
  if (neighborhood == Neighborhood::ReplacePath) {
    return replace_path_neighbor(solution, penalty_lambda, deadline);
  }
  return best_customer_neighbor(solution, neighborhood, penalty_lambda,
                                deadline, cooperative_physical_customers);
}

void VariableNeighborhoodDescent::order_relocate_moves(
    const std::vector<std::vector<int>>& current) {
  if (relocate_moves_.size() < 2) return;
  const auto& customers = instance_.customer_ids();
  // Montoya IDs are 1..N in order. Retain the original comparison for other
  // layouts, rather than assigning a new rank that would change tuple order.
  for (std::size_t i = 0; i < customers.size(); ++i) {
    if (customers[i] != static_cast<int>(i + 1)) {
      std::sort(relocate_moves_.begin(), relocate_moves_.end());
      relocate_moves_.erase(std::unique(relocate_moves_.begin(),
                                        relocate_moves_.end()),
                             relocate_moves_.end());
      return;
    }
  }
  order_gap_offsets_.resize(current.size() + 1);
  std::size_t columns = 0;
  for (std::size_t route = 0; route < current.size(); ++route) {
    order_gap_offsets_[route] = columns;
    columns += current[route].size() + 1;
  }
  order_gap_offsets_.back() = columns++;
  if (customers.size() > relocate_order_flags_.max_size() / columns) {
    throw std::length_error("Relocate ordering workspace is too large");
  }
  relocate_order_flags_.assign(customers.size() * columns, 0);
  for (const auto& [customer, route, gap] : relocate_moves_) {
    relocate_order_flags_[(static_cast<std::size_t>(customer) - 1) * columns +
                          order_gap_offsets_[route] + gap] = 1;
  }
  relocate_moves_.clear();
  for (std::size_t rank = 0; rank < customers.size(); ++rank) {
    for (std::size_t route = 0; route <= current.size(); ++route) {
      const std::size_t count = route == current.size()
                                    ? 1 : current[route].size() + 1;
      for (std::size_t gap = 0; gap < count; ++gap) {
        if (relocate_order_flags_[rank * columns + order_gap_offsets_[route] + gap]) {
          relocate_moves_.emplace_back(customers[rank], route, gap);
        }
      }
    }
  }
}

void VariableNeighborhoodDescent::order_two_opt_moves(
    const std::vector<std::vector<int>>& current) {
  if (two_opt_moves_.size() < 2) return;
  order_gap_offsets_.resize(current.size());
  std::size_t columns = 0;
  for (std::size_t route = 0; route < current.size(); ++route) {
    order_gap_offsets_[route] = columns;
    columns += current[route].size() + 1;
  }
  if (columns == 0 || columns > two_opt_order_flags_.max_size() / 2 / columns) {
    throw std::length_error("2-OPT ordering workspace is too large");
  }
  two_opt_order_flags_.assign(2 * columns * columns, 0);
  for (const auto& [left, right, left_cut, right_cut, reversed] : two_opt_moves_) {
    two_opt_order_flags_[2 * ((order_gap_offsets_[left] + left_cut) * columns +
                               order_gap_offsets_[right] + right_cut) + reversed] = 1;
  }
  two_opt_moves_.clear();
  // The nesting matches the original five-field tuple's lexicographic order.
  for (std::size_t left = 0; left < current.size(); ++left) {
    for (std::size_t right = left + 1; right < current.size(); ++right) {
      for (std::size_t left_cut = 0; left_cut <= current[left].size(); ++left_cut) {
        for (std::size_t right_cut = 0; right_cut <= current[right].size(); ++right_cut) {
          const std::size_t index = 2 * ((order_gap_offsets_[left] + left_cut) * columns +
                                        order_gap_offsets_[right] + right_cut);
          for (const bool reversed : {false, true}) {
            if (two_opt_order_flags_[index + reversed]) {
              two_opt_moves_.emplace_back(left, right, left_cut, right_cut, reversed);
            }
          }
        }
      }
    }
  }
}

bool VariableNeighborhoodDescent::described_move_less(
    const MoveDescriptor& left, const MoveDescriptor& right,
    const std::vector<std::vector<int>>& current) {
  // The original tie-break is lexicographic over (route slot, customer vector)
  // pairs. Traverse the same vectors through their slices without allocating.
  for (std::size_t change = 0;
       change < std::min(left.change_count, right.change_count); ++change) {
    const auto& a = left.changes[change];
    const auto& b = right.changes[change];
    if (a.route != b.route) return a.route < b.route;
    std::size_t a_slice = 0, b_slice = 0, a_offset = 0, b_offset = 0;
    while (true) {
      while (a_slice < a.slice_count &&
             a_offset == a.slices[a_slice].end - a.slices[a_slice].begin) {
        ++a_slice;
        a_offset = 0;
      }
      while (b_slice < b.slice_count &&
             b_offset == b.slices[b_slice].end - b.slices[b_slice].begin) {
        ++b_slice;
        b_offset = 0;
      }
      if (a_slice == a.slice_count || b_slice == b.slice_count) {
        if (a_slice == a.slice_count && b_slice == b.slice_count) break;
        return a_slice == a.slice_count;
      }
      const auto& as = a.slices[a_slice];
      const auto& bs = b.slices[b_slice];
      const int ac = current[as.route_index][
          as.reversed ? as.end - 1 - a_offset : as.begin + a_offset];
      const int bc = current[bs.route_index][
          bs.reversed ? bs.end - 1 - b_offset : bs.begin + b_offset];
      if (ac != bc) return ac < bc;
      ++a_offset;
      ++b_offset;
    }
  }
  return left.change_count < right.change_count;
}

std::optional<Solution> VariableNeighborhoodDescent::best_customer_neighbor(
    const Solution& solution, Neighborhood neighborhood, double penalty_lambda,
    std::optional<std::chrono::steady_clock::time_point> deadline,
    bool cooperative_physical_customers) {
  const bool physical_reverse_mode = cooperative_physical_customers &&
      profile_ == SearchProfile::Paper &&
      neighborhood == Neighborhood::TwoOptStar;
  const bool physical_relocate_mode = cooperative_physical_customers &&
      profile_ == SearchProfile::Paper && neighborhood == Neighborhood::Relocate;
  ++statistics_.customer_neighborhood_calls;
  const auto current = solution.customer_sequences();
  // Physical operators and support-memory hits do not need reconstructed
  // customer-gap labels. Materialize that legacy input only at its first use.
  std::optional<PlanFactory::GapMap> reference_gaps;
  const auto ensure_reference_gaps = [&]() -> const PlanFactory::GapMap& {
    if (!reference_gaps.has_value()) reference_gaps = factory_.gap_map(solution);
    return *reference_gaps;
  };
  bool sequence_scorer_ready = false;
  const auto ensure_sequence_scorer = [&]() -> const RouteSequenceScorer& {
    if (!sequence_scorer_ready) {
      const auto& gaps = ensure_reference_gaps();
      if (sequence_scorer_ == nullptr) {
        sequence_scorer_ = std::make_unique<RouteSequenceScorer>(
            instance_, solution, gaps);
      } else {
        // Reset on the first actual use in EVERY pass, even when the physical
        // calls in between left this persistent object pointing at old routes.
        sequence_scorer_->reset(solution, gaps);
      }
      sequence_scorer_ready = true;
    }
    return *sequence_scorer_;
  };
  auto prepared_plan = [&](const std::vector<int>& sequence) {
    const auto& gaps_by_anchor = ensure_reference_gaps();
    auto gaps = direct_gaps(sequence.size());
    std::vector<int> anchors;
    anchors.reserve(sequence.size() + 2);
    anchors.push_back(instance_.depot().id);
    anchors.insert(anchors.end(), sequence.begin(), sequence.end());
    anchors.push_back(instance_.depot().id);
    for (std::size_t gap = 0; gap + 1 < anchors.size(); ++gap) {
      const auto known =
          gaps_by_anchor.find(AnchorPair{anchors[gap], anchors[gap + 1]});
      if (known != gaps_by_anchor.end()) gaps[gap] = known->second;
    }
    return factory_.make_plan(sequence, std::move(gaps));
  };
  std::vector<double> current_route_costs;
  current_route_costs.reserve(solution.plans.size());
  for (const Plan& plan : solution.plans) {
    current_route_costs.push_back(
        generalized_cost(*plan.evaluation, instance_, penalty_lambda));
  }

  struct CustomerLocation {
    std::size_t route{};
    std::size_t position{};
    std::size_t customer_slot{};
  };
  const std::size_t customer_count = instance_.customer_ids().size();
  std::vector<CustomerLocation> location(customer_count);
  std::vector<bool> customer_seen(customer_count, false);
  bool canonical_physical_input = true;
  std::vector<std::vector<std::size_t>> current_ranks;
  current_ranks.reserve(current.size());
  for (std::size_t route = 0; route < current.size(); ++route) {
    std::vector<std::size_t> route_ranks;
    route_ranks.reserve(current[route].size());
    for (std::size_t position = 0; position < current[route].size();
         ++position) {
      const std::size_t rank =
          granular_.customer_rank(current[route][position]);
      if (customer_seen[rank]) canonical_physical_input = false;
      customer_seen[rank] = true;
      location[rank] = CustomerLocation{route, position, rank};
      route_ranks.push_back(rank);
    }
    current_ranks.push_back(std::move(route_ranks));
  }
  canonical_physical_input = canonical_physical_input &&
      std::all_of(customer_seen.begin(), customer_seen.end(),
                  [](bool seen) { return seen; });
  // A physical visit sequence identifies the scoring input only when its
  // customer/gap metadata agree and the partition has unique customers.
  // Unsupported public inputs retain the original uncached scoring path.
  std::vector<std::uint64_t> physical_ids(current.size(), 0);
  if (canonical_physical_input) {
    for (const Plan& plan : solution.plans) {
      if (plan.gaps.size() != plan.customer_ids.size() + 1) {
        canonical_physical_input = false;
        break;
      }
      // expand_route stores the interior visits; depot legs are implicit.
      std::size_t visit = 0;
      const auto matches = [&](int id) {
        return visit < plan.route.visits.size() &&
               plan.route.visits[visit++] == id;
      };
      for (std::size_t gap = 0; gap < plan.gaps.size(); ++gap) {
        for (const int station : plan.gaps[gap]) {
          if (!matches(station)) canonical_physical_input = false;
        }
        if (gap < plan.customer_ids.size() && !matches(plan.customer_ids[gap])) {
          canonical_physical_input = false;
        }
      }
      if (visit != plan.route.visits.size()) canonical_physical_input = false;
      if (!canonical_physical_input) break;
    }
    if (canonical_physical_input && physical_summary_budget_ != 0) {
      for (std::size_t route = 0; route < current.size(); ++route) {
        const auto& visits = solution.plans[route].route.visits;
        // Bound registry payload too; longer routes are evaluated normally.
        if (visits.size() > 256) continue;
        if (const auto id = physical_route_ids_.get(visits); id.has_value()) {
          physical_ids[route] = *id;
        } else if (next_physical_route_id_ != 0) {
          physical_ids[route] = next_physical_route_id_++;
          physical_route_ids_.put(visits, physical_ids[route]);
        }
      }
    }
  }
  // Unsupported metadata still invalidates the whole pass. Physical-length
  // limits are isolated below to the supports that actually touch a long route.
  // Weak evaluation owners retain no charge/SOC vectors.
  if (!canonical_physical_input) {
    move_memory_.clear();
    route_versions_.clear();
    next_route_version_ = 1;
    move_memory_usable_ = false;
    move_memory_reusable_ = false;
  }
  std::vector<std::uint64_t> current_route_versions(current.size(), 0);
  if (move_memory_usable_) {
    // Reserve enough IDs for this whole pass before assigning any of them.
    if (next_route_version_ > std::numeric_limits<std::uint64_t>::max() -
                                  current.size()) {
      move_memory_.clear();
      route_versions_.clear();
      next_route_version_ = 1;
    }
    for (std::size_t route = 0; route < current.size(); ++route) {
      const Plan& plan = solution.plans[route];
      // Zero is an uncacheable support version, never a shared long-route ID.
      // Keep the existing bounded raw key domain; do not retain long keys.
      if (plan.route.visits.size() > 256) continue;
      const auto known = route_versions_.get(plan.route.visits);
      if (known.has_value() && known->evaluation.lock() == plan.evaluation) {
        current_route_versions[route] = known->id;
      } else {
        current_route_versions[route] = next_route_version_++;
        route_versions_.put(plan.route.visits, RouteStateVersion{
            .id = current_route_versions[route], .evaluation = plan.evaluation});
      }
    }
  }
  const auto customer_location = [&](int customer_id) -> CustomerLocation {
    return location[granular_.customer_rank(customer_id)];
  };
  const auto add_slice = [&](DescribedRouteChange& change,
                             std::size_t source_route, std::size_t begin,
                             std::size_t end, bool reversed = false)
#if defined(__GNUC__) || defined(__clang__)
      __attribute__((always_inline))
#endif
  {
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

  // Only the active changes/slices are read. Reuse initialized storage within
  // this pass instead of transferring full arrays for every descriptor. A
  // remembered winner still owns its complete copy before the next reset.
  MoveDescriptor candidate_move;
  const auto one_change = [&](std::size_t route) -> DescribedRouteChange& {
    candidate_move.change_count = 1;
    candidate_move.physical_reverse_two_opt = false;
    candidate_move.physical_forward_two_opt = false;
    candidate_move.physical_relocate = false;
    auto& change = candidate_move.changes[0];
    change.route = route;
    change.slice_count = 0;
    return change;
  };
  const auto two_changes = [&](std::size_t first_route,
                               std::size_t second_route) {
    if (first_route == second_route) {
      throw std::logic_error("move descriptor changes one route twice");
    }
    candidate_move.change_count = 2;
    candidate_move.physical_reverse_two_opt = false;
    candidate_move.physical_forward_two_opt = false;
    candidate_move.physical_relocate = false;
    auto& first = candidate_move.changes[first_route < second_route ? 0 : 1];
    auto& second = candidate_move.changes[first_route < second_route ? 1 : 0];
    first.route = first_route;
    second.route = second_route;
    first.slice_count = second.slice_count = 0;
    return std::tie(first, second);
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
    return described_move_less(left, right, current);
  };

  const std::size_t route_count = current.size();
  const std::size_t support_count =
      route_count * route_count + 2 * route_count;
  // Pin each queried support table once per pass, rather than locking a weak
  // pointer for every descriptor. LRU payload accounting still describes only
  // retained cache entries; evicted tables can live until this pass returns.
  std::vector<std::shared_ptr<PhysicalSummaryTable>> physical_tables_in_pass;
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
      memory.physical_summaries = nullptr;
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
                      .opens_route = opens_route,
                      .physical_customer_mode = physical_reverse_mode || physical_relocate_mode};
    key.route_indices[0] = first;
    key.route_versions[0] = current_route_versions[first];
    if (second.has_value()) {
      key.route_indices[key.route_count] = *second;
      key.route_versions[key.route_count] = current_route_versions[*second];
      ++key.route_count;
    }
    return key;
  };

  const auto can_remember_support = [&](const MoveMemoryKey& key) {
    if (!move_memory_usable_) return false;
    for (std::size_t index = 0; index < key.route_count; ++index) {
      if (key.route_versions[index] == 0) return false;
    }
    return true;
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
      const auto cached = can_remember_support(key)
          ? move_memory_.get(key) : std::nullopt;
      if (cached.has_value()) {
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
  std::vector<std::optional<PhysicalScoreSummary>> relocate_source_summaries(
      neighborhood == Neighborhood::Relocate ? customer_count : 0);
  const auto physical_table = [&](SupportWorkingMemory& memory,
                                  const MoveSupport& support)
      -> PhysicalSummaryTable* {
    if (physical_ids[support.first] == 0 ||
        (support.second.has_value() && physical_ids[*support.second] == 0)) {
      return nullptr;
    }
    if (memory.physical_summaries != nullptr) return memory.physical_summaries;
    const PhysicalSupportKey key{
        .neighborhood = neighborhood,
        .routes = {physical_ids[support.first],
                   support.second.has_value() ? physical_ids[*support.second] : 0},
        .opens_route = support.opens_route,
        .physical_customer_mode = physical_reverse_mode || physical_relocate_mode};
    if (auto cached = physical_summaries_.get(key); cached.has_value()) {
      physical_tables_in_pass.push_back(std::move(*cached));
      memory.physical_summaries = physical_tables_in_pass.back().get();
      return memory.physical_summaries;
    }
    const std::size_t first_size = current[support.first].size();
    const std::size_t second_size = support.second.has_value()
        ? current[*support.second].size() : first_size;
    std::size_t row_count{};
    if (neighborhood == Neighborhood::Swap) {
      row_count = first_size * second_size;
    } else if (neighborhood == Neighborhood::Relocate) {
      row_count = support.opens_route ? first_size
          : support.second.has_value()
              ? first_size * (second_size + 1) + second_size * (first_size + 1)
              : first_size * (first_size + 1);
    } else {
      row_count = 2 * (first_size + 1) * (second_size + 1) +
          static_cast<std::size_t>(physical_reverse_mode);
    }
    const std::size_t bytes = row_count * sizeof(PhysicalSummaryRow);
    if (bytes == 0 || bytes > physical_summary_budget_) return nullptr;
    while (physical_summary_bytes_ + bytes > physical_summary_budget_ ||
           physical_summaries_.size() == physical_summaries_.capacity()) {
      auto evicted = physical_summaries_.discard_oldest();
      if (!evicted.has_value()) break;
      physical_summary_bytes_ -= (*evicted)->rows.size() * sizeof(PhysicalSummaryRow);
    }
    auto table = std::make_shared<PhysicalSummaryTable>();
    table->rows.resize(row_count);
    physical_summary_bytes_ += bytes;
    physical_summaries_.put(key, table);
    statistics_.physical_summary_cache_groups = physical_summaries_.size();
    statistics_.physical_summary_cache_bytes = physical_summary_bytes_;
    statistics_.physical_summary_cache_peak_bytes = std::max(
        statistics_.physical_summary_cache_peak_bytes,
        static_cast<std::uint64_t>(physical_summary_bytes_));
    statistics_.physical_summary_cache_evictions = physical_summaries_.evictions();
    physical_tables_in_pass.push_back(std::move(table));
    memory.physical_summaries = physical_tables_in_pass.back().get();
    return memory.physical_summaries;
  };
  NonlinearSequenceEvaluator opening_evaluator(instance_);
  // Profiles live only in this neighborhood pass. A later accepted move gets
  // fresh directional labels; no lambda-dependent score is stored here.
  std::vector<std::unique_ptr<PhysicalReverseTwoOptProfile>> two_opt_profiles(
      physical_reverse_mode ? route_count : 0);
  const auto physical_two_opt_profile = [&](std::size_t route)
      -> const PhysicalReverseTwoOptProfile& {
    auto& profile = two_opt_profiles.at(route);
    if (!profile) {
      profile = std::make_unique<PhysicalReverseTwoOptProfile>(
          instance_, solution.plans.at(route), opening_evaluator);
#if ILS_SP_FORWARD_ONLY_TWO_OPT
      statistics_.two_opt_reverse_label_slots_saved +=
          2 * (solution.plans.at(route).route.visits.size() + 1);
#endif
    }
    return *profile;
  };
  // Reuse only accumulator storage within this sequential neighborhood pass.
  // The first initialized append copy-assigns every field before any use.
  ChargingSequenceLabel physical_two_opt_scratch;
  const auto physical_two_opt_score = [&](const MoveDescriptor& move,
                                          std::size_t child) {
    ChargingSequenceLabel& label = physical_two_opt_scratch;
    label.initialized = false;
    int previous = instance_.depot().id;
    if (move.physical_forward_two_opt) {
      // Child 0 is A-prefix + B-suffix; child 1 is B-prefix + A-suffix.
      // Ordinary cut n stays before the parent's trailing station chain.
      physical_two_opt_profile(move.changes[child].route).append_label(
          move.two_opt_cuts[child], true, false,
          label, previous, opening_evaluator);
      physical_two_opt_profile(move.changes[1 - child].route).append_label(
          move.two_opt_cuts[1 - child], false, false,
          label, previous, opening_evaluator);
    } else {
      physical_two_opt_profile(move.changes[0].route).append_label(
          move.two_opt_cuts[0], child == 0, child == 1,
          label, previous, opening_evaluator);
      physical_two_opt_profile(move.changes[1].route).append_label(
          move.two_opt_cuts[1], child == 0, child == 0,
          label, previous, opening_evaluator);
    }
    opening_evaluator.append(label, opening_evaluator.edge_label(
        previous, {}, instance_.depot().id));
    return opening_evaluator.evaluate(label);
  };
  const auto physical_two_opt_visits = [&](const MoveDescriptor& move,
                                           std::size_t child) {
    std::vector<int> visits;
    if (move.physical_forward_two_opt) {
      physical_two_opt_profile(move.changes[child].route).append_visits(
          move.two_opt_cuts[child], true, false, visits);
      physical_two_opt_profile(move.changes[1 - child].route).append_visits(
          move.two_opt_cuts[1 - child], false, false, visits);
    } else {
      physical_two_opt_profile(move.changes[0].route).append_visits(
          move.two_opt_cuts[0], child == 0, child == 1, visits);
      physical_two_opt_profile(move.changes[1].route).append_visits(
          move.two_opt_cuts[1], child == 0, child == 0, visits);
    }
    return visits;
  };
  const bool physical_swap = preserve_swap_station_visits_ &&
                             neighborhood == Neighborhood::Swap;
  const bool paper_two_opt = profile_ == SearchProfile::Paper &&
                             neighborhood == Neighborhood::TwoOptStar;
  // Reuse the immutable Instance certificate; only the FP environment and
  // lambda domain need checking in this pass. Timed calls keep full scoring.
  // Physical 2OPT children partition the two parents' physical fragments.
  // The existing support certificate bounds every child by their combined
  // visits/station occurrences, so its generalized cost is also nonnegative.
  const bool exact_child_rejection_supported =
      (physical_swap || paper_two_opt) &&
      !deadline.has_value() && &path_sampler_.instance_ == &instance_ &&
      path_sampler_.replacement_summary_statistics().bound_certified &&
      supported_bound_fp_environment() &&
      nonlinear_travel_bound(0.0, 0.0, 0, instance_.vehicle(),
                             penalty_lambda).has_value();
  // Construct a route profile only on a physical-summary miss. Reuse it while
  // this route slot's entire physical sequence is unchanged, including CS.
  // The solution and route_count are immutable throughout this best-neighbor
  // pass. Validate each persistent profile once, then reuse its stable slot.
  std::vector<unsigned char> swap_profile_ready(route_count, 0);
  const auto swap_profile = [&](std::size_t route) -> SwapPhysicalProfile& {
    if (swap_profile_ready.at(route)) {
      ++statistics_.physical_profile_repeat_checks_saved;
      return swap_profiles_.at(route);
    }
    if (swap_profiles_.size() != route_count) swap_profiles_.resize(route_count);
    auto& profile = swap_profiles_.at(route);
    const Plan& plan = solution.plans.at(route);
    if (!profile.scorer || profile.visits != plan.route.visits) {
      SwapPhysicalProfile fresh;
      fresh.visits = plan.route.visits;
      for (std::size_t slot = 0; slot < fresh.visits.size(); ++slot) {
        if (instance_.node(fresh.visits[slot]).kind == NodeKind::Customer) {
          fresh.customer_slots.push_back(slot);
        }
      }
      fresh.scorer = std::make_unique<PhysicalRouteScorer>(instance_, plan, true);
      profile = std::move(fresh);
    }
    swap_profile_ready[route] = 1;
    return profile;
  };
  const auto mark_physical_relocate = [&](std::size_t source_route,
                                           std::size_t source_position,
                                           std::size_t target_route,
                                           std::size_t customer_gap) {
    if (!physical_relocate_mode || target_route == route_count) return;
    candidate_move.physical_relocate = true;
    candidate_move.relocate_source_route = source_route;
    candidate_move.relocate_target_route = target_route;
    candidate_move.relocate_source_position =
        swap_profile(source_route).customer_slots.at(source_position);
    // Keep the original Gamma/customer-gap list. Only its physical realization
    // changes: insert immediately after the predecessor, before the CS chain.
    candidate_move.relocate_target_slot = customer_gap == 0 ? 0
        : swap_profile(target_route).customer_slots.at(customer_gap - 1) + 1;
  };
  auto consider_move = [&](const MoveDescriptor& move,
                           std::optional<ReusableRelocateSource> reusable,
                           std::size_t physical_row,
                           bool operator_allowed) {
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
    if (!operator_allowed) {
      ++statistics_.granular_rejections;
      return;
    }

    const auto table = physical_table(memory, support);
    PhysicalSummaryRow* row = table ? &table->rows.at(physical_row) : nullptr;
    // Literal existing-route Relocate does not alter projected opening.
    const bool projected_opening =
        project_opening_station_visits_ && support.opens_route;
    const auto evaluate_physical = [&](std::size_t index,
                                       std::span<const IncumbentRouteSlice> slices) {
      const unsigned char bit = static_cast<unsigned char>(1U << index);
      if (row != nullptr && (row->valid & bit) != 0) {
        ++statistics_.physical_summary_hits;
        return row->summaries[index];
      }
      ++statistics_.physical_summary_misses;
      const RouteCostSummary evaluated = [&]() {
        if (move.physical_reverse_two_opt || move.physical_forward_two_opt)
          return physical_two_opt_score(move, index);
        if (move.physical_relocate) {
          auto& source_profile = swap_profile(move.relocate_source_route);
          if (move.relocate_source_route == move.relocate_target_route) {
            return source_profile.scorer->evaluate_relocate(
                move.relocate_source_position, move.relocate_target_slot);
          }
          if (move.changes[index].route == move.relocate_source_route) {
            return source_profile.scorer->evaluate_removal(move.relocate_source_position);
          }
          return swap_profile(move.relocate_target_route).scorer->evaluate_insertion(
              move.relocate_target_slot,
              source_profile.visits.at(move.relocate_source_position));
        }
        if (physical_swap) {
          if (!support.second.has_value()) {
            auto& profile = swap_profile(support.first);
            const std::size_t size = current[support.first].size();
            return profile.scorer->evaluate_swap(
                profile.customer_slots.at(physical_row / size),
                profile.customer_slots.at(physical_row % size));
          }
          const std::size_t second = *support.second;
          const std::size_t width = current[second].size();
          const std::size_t first_position = physical_row / width;
          const std::size_t second_position = physical_row % width;
          const bool first_child = move.changes[index].route == support.first;
          const std::size_t route = first_child ? support.first : second;
          auto& profile = swap_profile(route);
          return profile.scorer->evaluate_replacement(
              profile.customer_slots.at(first_child ? first_position : second_position),
              first_child ? current[second][second_position]
                          : current[support.first][first_position]);
        }
        return projected_opening
            ? opening_evaluator.evaluate(opening_evaluator.route_label(Route{
                  project_opening_visits(
                      instance_, solution.plans[support.first],
                      current[support.first][physical_row],
                      move.changes[index].route == route_count)}))
            : ensure_sequence_scorer().evaluate(slices);
      }();
      const PhysicalScoreSummary summary(evaluated);
      if (row != nullptr) {
        row->summaries[index] = summary;
        row->valid |= bit;
      }
      return summary;
    };
    double candidate_cost = 0.0;
    for (std::size_t index = 0; index < move.change_count; ++index) {
      const DescribedRouteChange& change = move.changes[index];
      if (change.slice_count == 0) continue;
      // Child 0 has already been scored and cached exactly in the original
      // order. Never skip a child-1 row hit or publish a bound as an exact row.
      if (index == 1 && row != nullptr && (row->valid & 2U) == 0 &&
          exact_child_rejection_supported && support.second.has_value() &&
          !support.opens_route && move.change_count == 2 &&
          move.changes[0].slice_count != 0 &&
          physical_ids[support.first] != 0 &&
          physical_ids[*support.second] != 0 &&
          std::isfinite(candidate_cost) && candidate_cost >= 0.0 &&
          std::isfinite(memory.entry.candidate_cost) &&
          memory.entry.candidate_cost >= 0.0 &&
          candidate_cost > memory.entry.candidate_cost &&
          candidate_cost - memory.entry.candidate_cost > kCostTolerance) {
        // 2OPT's support bit certifies both children using the combined
        // parent domain. Keep its sequence scorer separate from Swap.
        if (paper_two_opt &&
            (customer_open_support_flags_[support.first * route_count +
                                           *support.second] & 2U) != 0) {
          return;
        }
        const auto& first = solution.plans[support.first];
        const auto& second = solution.plans[*support.second];
        // Canonical physical IDs bound each source to 256 visits. Replacement
        // preserves its station count; keep both kernels in the certified
        // inline domain and skip only a second child with charging work.
        if (physical_swap &&
            first.route.visits.size() >= first.customer_ids.size() &&
            first.route.visits.size() - first.customer_ids.size() <= 62 &&
            second.route.visits.size() > second.customer_ids.size() &&
            second.route.visits.size() - second.customer_ids.size() <= 62) {
          // Certified C1 >= 0 implies fl(C0 + C1) >= C0. This exact positive
          // C0-B difference excludes both improvement and tolerance ties.
          return;
        }
      }
      const auto slices = std::span<const IncumbentRouteSlice>{
          change.slices.data(), change.slice_count};
      // The legacy deletion memo drops newly joined station paths. It cannot
      // supply the projected source, whose visits retain those same stations.
      // Its separate opening-support physical row still caches the new result.
      if (!projected_opening && reusable.has_value() &&
          change.route == reusable->route) {
        auto& summary =
            relocate_source_summaries.at(reusable->customer_slot);
        if (!summary.has_value()) {
          summary = evaluate_physical(index, slices);
          ++statistics_.slice_routes_evaluated;
          ++statistics_.relocate_source_summary_misses;
        } else {
          ++statistics_.relocate_source_summary_hits;
          if (row != nullptr) {
            row->summaries[index] = *summary;
            row->valid |= static_cast<unsigned char>(1U << index);
          }
        }
        candidate_cost += summary->generalized_cost(instance_, penalty_lambda);
      } else {
        ++statistics_.slice_routes_evaluated;
        candidate_cost += evaluate_physical(index, slices).generalized_cost(
            instance_, penalty_lambda);
      }
    }
    ++statistics_.moves_evaluated;
    if (candidate_cost < memory.entry.candidate_cost - kCostTolerance ||
        (std::abs(candidate_cost - memory.entry.candidate_cost) <=
             kCostTolerance &&
         move_less(move, memory.entry.move))) {
      memory.entry = MoveMemoryEntry{.candidate_cost = candidate_cost,
                                     .move = move};
    }
  };
  const auto consider = [&](std::size_t physical_row,
                             bool operator_allowed) {
    consider_move(candidate_move, std::nullopt, physical_row, operator_allowed);
  };

  if (neighborhood == Neighborhood::Swap) {
    std::vector<std::pair<int, int>>& customer_pairs = swap_customer_pairs_;
    customer_pairs.clear();
    customer_pairs.reserve(granular_.promising_pairs().size() * 4);
    if (route_count != 0 &&
        route_count > customer_open_support_flags_.max_size() / route_count) {
      throw std::length_error("Swap support workspace is too large");
    }
    customer_open_support_flags_.assign(route_count * route_count, 0);
    // Keep the original support-memory lookup order, including closed hits.
    // A promising neighbor already identifies its target route; do not scan
    // the same source list again for every open target route.
    for (std::size_t first_route = 0; first_route < route_count;
         ++first_route) {
      if (!support_is_closed(first_route, std::nullopt, false)) {
        customer_open_support_flags_[first_route * route_count + first_route] = 1;
      }
      for (std::size_t second_route = first_route + 1;
           second_route < route_count; ++second_route) {
        if (!support_is_closed(first_route, second_route, false)) {
          customer_open_support_flags_[first_route * route_count + second_route] = 1;
          customer_open_support_flags_[second_route * route_count + first_route] = 1;
        }
      }
    }
    const auto& customers = instance_.customer_ids();
    bool numeric_order = true;
    for (std::size_t i = 0; i < customers.size(); ++i) {
      if (customers[i] != static_cast<int>(i + 1)) {
        numeric_order = false;
        break;
      }
    }
    if (numeric_order) {
      if (!customers.empty() &&
          customers.size() > swap_order_flags_.max_size() / customers.size()) {
        throw std::length_error("Swap ordering workspace is too large");
      }
      swap_order_flags_.assign(customers.size() * customers.size(), 0);
    }
    for (std::size_t source_route = 0; source_route < route_count;
         ++source_route) {
      const auto open = customer_open_support_flags_.begin() +
                        source_route * route_count;
      if (std::none_of(open, open + route_count,
                       [](unsigned char flag) { return flag != 0; })) continue;
      for (std::size_t position = 0;
           position < current[source_route].size(); ++position) {
        const int candidate_customer = current[source_route][position];
        const auto add_from_adjacent = [&](int adjacent_customer) {
          for (const int promising_customer :
               granular_.promising_neighbors_of(adjacent_customer)) {
            const auto target_route = customer_location(promising_customer).route;
            if (!open[target_route] || candidate_customer == promising_customer) {
              continue;
            }
            const auto [left, right] =
                std::minmax(candidate_customer, promising_customer);
            if (numeric_order) {
              swap_order_flags_[(static_cast<std::size_t>(left) - 1) *
                                    customers.size() +
                                static_cast<std::size_t>(right) - 1] = 1;
            } else {
              customer_pairs.emplace_back(left, right);
            }
          }
        };
        if (position > 0) {
          add_from_adjacent(current[source_route][position - 1]);
        }
        if (position + 1 < current[source_route].size()) {
          add_from_adjacent(current[source_route][position + 1]);
        }
      }
    }
    if (numeric_order) {
      // The upper-triangle scan is exactly sort+unique on numeric ID pairs.
      for (std::size_t left = 0; left < customers.size(); ++left) {
        for (std::size_t right = left + 1; right < customers.size(); ++right) {
          if (swap_order_flags_[left * customers.size() + right]) {
            customer_pairs.emplace_back(customers[left], customers[right]);
          }
        }
      }
    } else {
      std::sort(customer_pairs.begin(), customer_pairs.end());
      customer_pairs.erase(
          std::unique(customer_pairs.begin(), customer_pairs.end()),
          customer_pairs.end());
    }
    for (const auto& [left_customer, right_customer] : customer_pairs) {
      const CustomerLocation left = customer_location(left_customer);
      const CustomerLocation right = customer_location(right_customer);
      // Correlated sets contain customers only. A depot predecessor contributes
      // no qualifying arc; the other entry can still qualify the Swap.
      const bool entry_allowed =
          (left.position > 0 && granular_.promising_by_rank(
              current_ranks[left.route][left.position - 1], right.customer_slot)) ||
          (right.position > 0 && granular_.promising_by_rank(
              current_ranks[right.route][right.position - 1], left.customer_slot));
      if (left.route == right.route) {
        const std::size_t low = std::min(left.position, right.position);
        const std::size_t high = std::max(left.position, right.position);
        if (high < low + 2) continue;
        auto& changed = one_change(left.route);
        add_slice(changed, left.route, 0, low);
        add_slice(changed, left.route, high, high + 1);
        add_slice(changed, left.route, low + 1, high);
        add_slice(changed, left.route, low, low + 1);
        add_slice(changed, left.route, high + 1,
                  current[left.route].size());
        consider(low * current[left.route].size() + high, entry_allowed);
      } else {
        auto [changed_left, changed_right] = two_changes(left.route, right.route);
        add_slice(changed_left, left.route, 0, left.position);
        add_slice(changed_left, right.route, right.position,
                  right.position + 1);
        add_slice(changed_left, left.route, left.position + 1,
                  current[left.route].size());
        add_slice(changed_right, right.route, 0, right.position);
        add_slice(changed_right, left.route, left.position,
                  left.position + 1);
        add_slice(changed_right, right.route, right.position + 1,
                  current[right.route].size());
        const auto& lower = left.route < right.route ? left : right;
        const auto& upper = left.route < right.route ? right : left;
        consider(lower.position * current[upper.route].size() + upper.position,
                 entry_allowed);
      }
    }
  } else if (neighborhood == Neighborhood::Relocate) {
    using RelocateMove = std::tuple<int, std::size_t, std::size_t>;
    std::vector<RelocateMove>& moves = relocate_moves_;
    moves.clear();
    moves.reserve(instance_.customer_ids().size() *
                  (4 * granular_.gamma_size() + 1));
    if (route_count != 0 &&
        route_count > customer_open_support_flags_.max_size() / route_count) {
      throw std::length_error("Relocate support workspace is too large");
    }
    customer_open_support_flags_.assign(route_count * route_count, 0);
    // Preserve intra/opening/pair memory queries in their original order.
    for (std::size_t first_route = 0; first_route < route_count;
         ++first_route) {
      if (!support_is_closed(first_route, std::nullopt, false)) {
        customer_open_support_flags_[first_route * route_count + first_route] = 1;
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
        customer_open_support_flags_[first_route * route_count + second_route] = 1;
        customer_open_support_flags_[second_route * route_count + first_route] = 1;
      }
    }
    for (std::size_t source_route = 0; source_route < route_count;
         ++source_route) {
      const auto open = customer_open_support_flags_.begin() +
                        source_route * route_count;
      if (std::none_of(open, open + route_count,
                       [](unsigned char flag) { return flag != 0; })) continue;
      for (const int customer : current[source_route]) {
        // N1 qualifies the incoming destination arc, not source closure.
        for (const int neighbor : granular_.promising_neighbors_of(customer)) {
          const CustomerLocation target = customer_location(neighbor);
          if (!open[target.route] ||
              !granular_.promising(neighbor, customer)) continue;
          moves.emplace_back(customer, target.route, target.position + 1);
        }
      }
    }
    order_relocate_moves(current);
    for (const auto& [customer, target_route, original_gap] : moves) {
      const CustomerLocation source = customer_location(customer);
      // HLNS permits Relocate to open a route. An empty destination has no
      // customer entry arc on which Gamma is defined; keep this structural
      // case separate from the unchanged existing-route entry restriction.
      // This boundary interpretation is recorded in audit 5.15.
      const bool entry_allowed = target_route == route_count ||
          (original_gap > 0 && granular_.promising_by_rank(
              current_ranks[target_route][original_gap - 1], source.customer_slot));
      if (target_route == source.route) {
        std::size_t insertion = original_gap;
        if (source.position < original_gap) --insertion;
        if (physical_relocate_mode) {
          const auto& slots = swap_profile(source.route).customer_slots;
          const std::size_t p = slots.at(source.position);
          const std::size_t q = original_gap == 0 ? 0 : slots.at(original_gap - 1) + 1;
          // Customer-order no-ops can move a customer across an existing CS.
          if (q - static_cast<std::size_t>(p < q) == p) continue;
        } else if (insertion == source.position) {
          continue;
        }
        auto& changed = one_change(source.route);
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
        mark_physical_relocate(source.route, source.position, target_route, original_gap);
        consider(source.position * (current[source.route].size() + 1) + original_gap,
                 entry_allowed);
      } else if (target_route == current.size()) {
        auto [changed_source, singleton] = two_changes(source.route, current.size());
        add_slice(changed_source, source.route, 0, source.position);
        add_slice(changed_source, source.route, source.position + 1,
                  current[source.route].size());
        add_slice(singleton, source.route, source.position,
                  source.position + 1);
        consider_move(
            candidate_move,
            ReusableRelocateSource{.route = source.route,
                                   .customer_slot = source.customer_slot},
            source.position, entry_allowed);
      } else {
        auto [changed_source, changed_target] = two_changes(source.route, target_route);
        add_slice(changed_source, source.route, 0, source.position);
        add_slice(changed_source, source.route, source.position + 1,
                  current[source.route].size());
        add_slice(changed_target, target_route, 0, original_gap);
        add_slice(changed_target, source.route, source.position,
                  source.position + 1);
        add_slice(changed_target, target_route, original_gap,
                  current[target_route].size());
        mark_physical_relocate(source.route, source.position, target_route, original_gap);
        const std::size_t lower_route = std::min(source.route, target_route);
        const std::size_t upper_route = std::max(source.route, target_route);
        const std::size_t physical_row = source.route == lower_route
            ? source.position * (current[upper_route].size() + 1) + original_gap
            : current[lower_route].size() * (current[upper_route].size() + 1) +
                source.position * (current[lower_route].size() + 1) + original_gap;
        consider_move(
            candidate_move,
            ReusableRelocateSource{.route = source.route,
                                   .customer_slot = source.customer_slot}, physical_row,
            entry_allowed);
      }
    }
  } else if (neighborhood == Neighborhood::TwoOptStar) {
    using TwoOptMove =
        std::tuple<std::size_t, std::size_t, std::size_t, std::size_t, bool>;
    std::vector<TwoOptMove>& moves = two_opt_moves_;
    moves.clear();
    moves.reserve(granular_.promising_pairs().size() * 4 +
                  current.size() * (current.size() - 1) / 2);
    if (route_count != 0 &&
        route_count > customer_open_support_flags_.max_size() / route_count) {
      throw std::length_error("2OPT support workspace is too large");
    }
    customer_open_support_flags_.assign(route_count * route_count, 0);
    for (std::size_t left_route = 0; left_route < current.size(); ++left_route) {
      for (std::size_t right_route = left_route + 1;
           right_route < current.size(); ++right_route) {
        if (support_is_closed(left_route, right_route, false)) continue;
        auto& flags =
            customer_open_support_flags_[left_route * route_count + right_route];
        flags = 1;
        if (!ILS_SP_FORWARD_ONLY_TWO_OPT && physical_reverse_mode &&
            !current[left_route].empty() &&
            !current[right_route].empty()) {
          // Complete the departure/return-depot boundary coverage of the
          // existing reverse inter-route reconnection. The new arcs have no
          // customer/customer pair for Gamma to qualify. Each child remains
          // nonempty: reverse(B), reverse(A), with every physical visit kept.
          // The opposite true-end case yields A,B unchanged and is omitted.
          auto [changed_left, changed_right] = two_changes(left_route, right_route);
          candidate_move.physical_reverse_two_opt = true;
          candidate_move.two_opt_cuts = {0, current[right_route].size() + 1};
          add_slice(changed_left, right_route, 0, current[right_route].size(), true);
          add_slice(changed_right, left_route, 0, current[left_route].size(), true);
          // Dedicated row: do not alias generic customer cut (0,n), whose
          // physical boundary precedes trailing stations. Role order remains
          // the support's ordered route slots in both summaries and winner.
          consider(2 * (current[left_route].size() + 1) *
                       (current[right_route].size() + 1), true);
        }
        if (!exact_child_rejection_supported || physical_ids[left_route] == 0 ||
            physical_ids[right_route] == 0) continue;
        const auto& left = solution.plans[left_route];
        const auto& right = solution.plans[right_route];
        const std::size_t left_visits = left.route.visits.size();
        const std::size_t right_visits = right.route.visits.size();
        if (left.customer_ids.empty() || right.customer_ids.empty() ||
            left_visits < left.customer_ids.size() ||
            right_visits < right.customer_ids.size() ||
            left_visits > 256 || right_visits > 256 - left_visits) continue;
        const std::size_t left_stations = left_visits - left.customer_ids.size();
        const std::size_t right_stations = right_visits - right.customer_ids.size();
        // Every child's inherited stations are a subset of these two parents;
        // reversed internal edges are direct. Bound occurrences, not IDs.
        if (left_stations <= 62 && right_stations <= 62 - left_stations &&
            (left_stations != 0 || right_stations != 0)) flags = 3;
      }
    }
    for (std::size_t left_route = 0; left_route < route_count; ++left_route) {
      const auto open = customer_open_support_flags_.begin() +
                        left_route * route_count;
      if (std::none_of(open + left_route + 1, open + route_count,
                       [](unsigned char flag) { return flag != 0; })) continue;
      for (std::size_t left_position = 0;
           left_position < current[left_route].size(); ++left_position) {
        const int left_customer = current[left_route][left_position];
        for (const int right_customer :
             granular_.promising_neighbors_of(left_customer)) {
          const CustomerLocation right = customer_location(right_customer);
          if (right.route <= left_route || !open[right.route]) continue;
          moves.emplace_back(left_route, right.route, left_position + 1,
                             right.position, false);
          moves.emplace_back(left_route, right.route, left_position,
                             right.position + 1, false);
#if !ILS_SP_FORWARD_ONLY_TWO_OPT
          moves.emplace_back(left_route, right.route, left_position + 1,
                             right.position + 1, true);
          moves.emplace_back(left_route, right.route, left_position,
                             right.position, true);
#endif
        }
      }
    }
    order_two_opt_moves(current);
    for (const auto& [left_route, right_route, left_cut, right_cut,
                      reverse_reconnection] : moves) {
      const auto& left = current_ranks[left_route];
      const auto& right = current_ranks[right_route];
      // Check the actual two directed reconnections. Reversed internal arcs
      // cannot qualify a move whose reconnection arcs are both unpromising.
      const bool entry_allowed = !reverse_reconnection
          ? ((left_cut > 0 && right_cut < right.size() &&
              granular_.promising_by_rank(left[left_cut - 1], right[right_cut])) ||
             (right_cut > 0 && left_cut < left.size() &&
              granular_.promising_by_rank(right[right_cut - 1], left[left_cut])))
          : ((left_cut > 0 && right_cut > 0 &&
              granular_.promising_by_rank(left[left_cut - 1], right[right_cut - 1])) ||
             (left_cut < left.size() && right_cut < right.size() &&
              granular_.promising_by_rank(left[left_cut], right[right_cut])));
      auto [changed_left, changed_right] = two_changes(left_route, right_route);
      if (!reverse_reconnection) {
        candidate_move.physical_forward_two_opt = physical_reverse_mode;
        candidate_move.two_opt_cuts = {left_cut, right_cut};
        // A customer-empty child is an unused vehicle: its station-only tour
        // contributes zero (empty slices) and is not materialized. The other
        // child's physical prefix/suffix, including CS visits, stays intact.
        add_slice(changed_left, left_route, 0, left_cut);
        add_slice(changed_left, right_route, right_cut,
                  current[right_route].size());
        add_slice(changed_right, right_route, 0, right_cut);
        add_slice(changed_right, left_route, left_cut,
                  current[left_route].size());
      } else {
#if ILS_SP_PHYSICAL_REVERSE_EMPTY_CHILD
        // Use the same physical cuts for reverse route closure as for every
        // other reconnection. Empty customer slices contribute zero and are
        // not materialized; the productive child keeps its physical fragments.
        candidate_move.physical_reverse_two_opt = physical_reverse_mode;
#else
        // Empty-customer children retain the existing legacy behavior. For
        // nonempty children, explicit cuts also retain station-only physical
        // fragments that are absent from the customer-slice descriptor.
        candidate_move.physical_reverse_two_opt = physical_reverse_mode &&
            left_cut + right_cut != 0 &&
            left_cut + right_cut != left.size() + right.size();
#endif
        candidate_move.two_opt_cuts = {left_cut, right_cut};
        add_slice(changed_left, left_route, 0, left_cut);
        add_slice(changed_left, right_route, 0, right_cut, true);
        add_slice(changed_right, left_route, left_cut,
                  current[left_route].size(), true);
        add_slice(changed_right, right_route, right_cut,
                  current[right_route].size());
      }
      consider(2 * (left_cut * (current[right_route].size() + 1) + right_cut) +
                   static_cast<std::size_t>(reverse_reconnection), entry_allowed);
    }
  }

  const bool neighborhood_completed = !expired(deadline);
  std::optional<MoveMemoryEntry> best_move;
  double best_delta = std::numeric_limits<double>::infinity();
  for (std::size_t support = 0; support < support_count; ++support) {
    SupportWorkingMemory& memory = support_working_memory_[support];
    if (memory.epoch != support_epoch_ || !memory.active) continue;
    if (can_remember_support(memory.key) &&
        memory.state != SupportState::Closed && neighborhood_completed) {
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
  const MoveSupport best_support = move_support(best_move->move);
  const bool projected_opening =
      project_opening_station_visits_ && best_support.opens_route;

  struct PreparedChange {
    std::size_t route{};
    std::optional<Plan> plan;
  };
  std::vector<PreparedChange> prepared;
  prepared.reserve(best_changes.size());
  for (const auto& [route, sequence] : best_changes) {
    std::optional<Plan> plan;
    if (!sequence.empty()) {
      if (best_move->move.physical_reverse_two_opt ||
          best_move->move.physical_forward_two_opt) {
        const std::size_t child = route == best_move->move.changes[0].route ? 0 : 1;
        plan = factory_.plan_visits(physical_two_opt_visits(best_move->move, child));
        if (plan->customer_ids != sequence) {
          throw std::logic_error("physical 2OPT customer descriptor differs");
        }
      } else if (best_move->move.physical_relocate) {
        const auto& move = best_move->move;
        auto visits = solution.plans.at(route).route.visits;
        const int customer = solution.plans.at(move.relocate_source_route)
            .route.visits.at(move.relocate_source_position);
        if (route == move.relocate_source_route) {
          visits.erase(visits.begin() +
                       static_cast<std::ptrdiff_t>(move.relocate_source_position));
        }
        if (route == move.relocate_target_route) {
          const std::size_t slot = move.relocate_target_slot -
              static_cast<std::size_t>(move.relocate_source_route == route &&
                  move.relocate_source_position < move.relocate_target_slot);
          visits.insert(visits.begin() + static_cast<std::ptrdiff_t>(slot), customer);
        }
        plan = factory_.plan_visits(std::move(visits));
        if (plan->customer_ids != sequence) {
          throw std::logic_error("physical Relocate customer descriptor differs");
        }
      } else if (physical_swap) {
        auto visits = solution.plans.at(route).route.visits;
        std::size_t customer = 0;
        for (int& id : visits) {
          if (instance_.node(id).kind == NodeKind::Customer) {
            id = sequence.at(customer++);
          }
        }
        if (customer != sequence.size()) throw std::logic_error("Swap materialization size");
        plan = factory_.plan_visits(std::move(visits));
      } else {
        plan = projected_opening
            ? factory_.plan_visits(project_opening_visits(
                  instance_, solution.plans[best_support.first],
                  best_changes.back().second.front(), route == route_count))
            : prepared_plan(sequence);
      }
      if (completes_vnd_winner(profile_)) {
        ++statistics_.linked_completion_attempts;
        const auto direct_gaps = plan->gaps;
        plan = path_sampler_.repair_new_gaps(
            std::move(*plan), ensure_reference_gaps(), penalty_lambda, deadline);
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
    const Solution& solution, double penalty_lambda,
    std::optional<std::chrono::steady_clock::time_point> deadline) {
  ++statistics_.replace_path_calls;
  if (solution.plans.empty()) return std::nullopt;
  const double current_cost =
      generalized_cost(solution, instance_, penalty_lambda);
  for (std::size_t attempt = 0; attempt < replace_path_attempts_; ++attempt) {
    // Preserve the original single-proposal deadline/RNG behavior. Only the
    // explicit repeated stage adds deadline checks before later proposals.
    if (attempt != 0 && expired(deadline)) break;
    ++statistics_.replace_path_proposals;
    auto proposal = path_sampler_.random_replacement_proposal(solution);
    if (!proposal.has_value()) {
      ++statistics_.replace_path_null_proposals;
      continue;
    }
    // Exactly the original generalized_cost(Solution) addition order. Avoid
    // a delta or subtract/add update, which can change floating-point ties.
    double candidate_cost = 0.0;
    for (std::size_t route = 0; route < solution.plans.size(); ++route) {
      const Plan& plan = route == proposal->route_index
          ? proposal->plan : solution.plans[route];
      candidate_cost += generalized_cost(*plan.evaluation, instance_, penalty_lambda);
    }
    if (candidate_cost < current_cost - kCostTolerance) {
      ++statistics_.replace_path_accepted;
      Solution candidate = solution;
      candidate.plans[proposal->route_index] = std::move(proposal->plan);
      return candidate;
    }
    ++statistics_.replace_path_nonimproving_proposals;
  }
  // The explicit four-proposal A experiment failed its G1 quality screen.
  // The active contract remains one proposal, including null proposals.
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
