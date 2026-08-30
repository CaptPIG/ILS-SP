#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>

#include "ils_sp/lru_cache.hpp"
#include "ils_sp/pool.hpp"
#include "ils_sp/search.hpp"
#include "ils_sp/xml.hpp"

namespace ils_sp {

struct PerturbationTestAccess {
  static Solution repair(Perturbation& perturbation,
                         DestroyedSolution destroyed,
                         RepairOperator operation, double penalty_lambda) {
    return perturbation.repair(std::move(destroyed), operation,
                               penalty_lambda);
  }
};

}  // namespace ils_sp

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

ils_sp::DestroyedSolution make_destroyed_fixture(
    const ils_sp::Solution& source, ils_sp::PlanFactory& factory,
    const std::vector<int>& removed) {
  const auto inherited = factory.gap_map(source);
  ils_sp::Solution remaining;
  for (const ils_sp::Plan& source_plan : source.plans) {
    std::vector<int> sequence = source_plan.customer_ids;
    std::erase_if(sequence, [&](int customer_id) {
      return std::find(removed.begin(), removed.end(), customer_id) !=
             removed.end();
    });
    if (!sequence.empty()) {
      remaining.plans.push_back(factory.plan_sequence(sequence, inherited));
    }
  }
  return ils_sp::DestroyedSolution{.remaining = std::move(remaining),
                                   .removed_customers = removed};
}

ils_sp::Solution exhaustive_linked_repair(
    ils_sp::DestroyedSolution destroyed, ils_sp::RepairOperator operation,
    const ils_sp::Instance& instance, ils_sp::PlanFactory& factory,
    ils_sp::PathSampler& path_sampler, double penalty_lambda) {
  struct Placement {
    double delta{};
    int customer{};
    std::size_t route{};
    std::size_t position{};
    ils_sp::Plan plan;
  };
  const auto better = [](const Placement& left, const Placement& right) {
    return std::tie(left.delta, left.customer, left.route, left.position) <
           std::tie(right.delta, right.customer, right.route,
                    right.position);
  };

  ils_sp::Solution current = std::move(destroyed.remaining);
  std::vector<int> unassigned = std::move(destroyed.removed_customers);
  while (!unassigned.empty()) {
    const auto inherited = factory.gap_map(current);
    std::optional<Placement> chosen;
    double best_regret = -1.0;
    for (const int customer_id : unassigned) {
      std::vector<Placement> placements;
      for (std::size_t route = 0; route <= current.plans.size(); ++route) {
        const std::size_t positions =
            route == current.plans.size()
                ? 1
                : current.plans[route].customer_ids.size() + 1;
        for (std::size_t position = 0; position < positions; ++position) {
          std::vector<int> sequence;
          double incumbent_cost = 0.0;
          if (route == current.plans.size()) {
            sequence = {customer_id};
          } else {
            sequence = current.plans[route].customer_ids;
            sequence.insert(
                sequence.begin() + static_cast<std::ptrdiff_t>(position),
                customer_id);
            incumbent_cost = ils_sp::generalized_cost(
                *current.plans[route].evaluation, instance, penalty_lambda);
          }
          ils_sp::Plan plan = factory.plan_sequence(sequence, inherited);
          plan = path_sampler.repair_new_gaps(
              std::move(plan), inherited, penalty_lambda);
          placements.push_back(Placement{
              .delta = ils_sp::generalized_cost(
                           *plan.evaluation, instance, penalty_lambda) -
                       incumbent_cost,
              .customer = customer_id,
              .route = route,
              .position = position,
              .plan = std::move(plan)});
        }
      }
      std::sort(placements.begin(), placements.end(), better);
      if (operation == ils_sp::RepairOperator::GreedyInsertion) {
        if (!chosen.has_value() || better(placements.front(), *chosen)) {
          chosen = std::move(placements.front());
        }
      } else {
        const double regret =
            placements.size() == 1
                ? std::numeric_limits<double>::infinity()
                : placements[1].delta - placements[0].delta;
        if (!chosen.has_value() ||
            regret > best_regret + ils_sp::kCostTolerance ||
            (std::abs(regret - best_regret) <= ils_sp::kCostTolerance &&
             better(placements.front(), *chosen))) {
          best_regret = regret;
          chosen = std::move(placements.front());
        }
      }
    }
    const int inserted = chosen->customer;
    if (chosen->route == current.plans.size()) {
      current.plans.push_back(std::move(chosen->plan));
    } else {
      current.plans[chosen->route] = std::move(chosen->plan);
    }
    std::erase(unassigned, inserted);
  }
  current.validate_partition(instance);
  return current;
}

ils_sp::Solution exhaustive_virtual_repair(
    ils_sp::DestroyedSolution destroyed, ils_sp::RepairOperator operation,
    const ils_sp::Instance& instance, ils_sp::PlanFactory& factory,
    ils_sp::PathSampler& path_sampler, double penalty_lambda) {
  struct Placement {
    double delta{};
    int customer{};
    std::size_t route{};
    std::size_t position{};
  };
  const auto better = [](const Placement& left, const Placement& right) {
    return std::tie(left.delta, left.customer, left.route, left.position) <
           std::tie(right.delta, right.customer, right.route,
                    right.position);
  };

  ils_sp::Solution current = std::move(destroyed.remaining);
  std::vector<int> unassigned = std::move(destroyed.removed_customers);
  while (!unassigned.empty()) {
    const auto inherited = factory.gap_map(current);
    std::optional<Placement> chosen;
    double best_regret = -1.0;
    for (const int customer_id : unassigned) {
      std::vector<Placement> placements;
      for (std::size_t route = 0; route <= current.plans.size(); ++route) {
        const std::size_t positions =
            route == current.plans.size()
                ? 1
                : current.plans[route].customer_ids.size() + 1;
        for (std::size_t position = 0; position < positions; ++position) {
          const ils_sp::Plan* source =
              route == current.plans.size() ? nullptr : &current.plans[route];
          const double incumbent_cost =
              source == nullptr
                  ? 0.0
                  : ils_sp::generalized_cost(*source->evaluation, instance,
                                             penalty_lambda);
          const auto completion = path_sampler.virtual_complete_insertion(
              source, position, customer_id, penalty_lambda);
          placements.push_back(Placement{
              .delta = completion.generalized_cost_h - incumbent_cost,
              .customer = customer_id,
              .route = route,
              .position = position});
        }
      }
      std::sort(placements.begin(), placements.end(), better);
      if (operation == ils_sp::RepairOperator::GreedyInsertion) {
        if (!chosen.has_value() || better(placements.front(), *chosen)) {
          chosen = placements.front();
        }
      } else {
        const double regret =
            placements.size() == 1
                ? std::numeric_limits<double>::infinity()
                : placements[1].delta - placements[0].delta;
        if (!chosen.has_value() ||
            regret > best_regret + ils_sp::kCostTolerance ||
            (std::abs(regret - best_regret) <= ils_sp::kCostTolerance &&
             better(placements.front(), *chosen))) {
          best_regret = regret;
          chosen = placements.front();
        }
      }
    }
    std::vector<int> sequence;
    if (chosen->route == current.plans.size()) {
      sequence = {chosen->customer};
    } else {
      sequence = current.plans[chosen->route].customer_ids;
      sequence.insert(
          sequence.begin() + static_cast<std::ptrdiff_t>(chosen->position),
          chosen->customer);
    }
    ils_sp::Plan direct = factory.plan_sequence(sequence, inherited);
    if (chosen->route == current.plans.size()) {
      current.plans.push_back(std::move(direct));
    } else {
      current.plans[chosen->route] = std::move(direct);
    }
    std::erase(unassigned, chosen->customer);
  }
  current.validate_partition(instance);
  return current;
}

std::vector<std::vector<int>> exhaustive_station_catalog(
    const ils_sp::Instance& instance, int origin_id, int target_id) {
  struct Candidate {
    std::vector<int> stations;
    double first_leg{};
    double last_leg{};
    double total{};
  };
  const auto& stations = instance.station_ids();
  require(stations.size() <= 8,
          "exhaustive station catalog oracle is only for small fixtures");
  std::vector<Candidate> feasible;
  std::vector<int> path;
  std::vector<bool> used(stations.size(), false);
  const auto enumerate = [&](auto&& self, int previous,
                             double prefix_distance) -> void {
    for (std::size_t index = 0; index < stations.size(); ++index) {
      if (used[index] ||
          instance.energy_wh(previous, stations[index]) >
              instance.vehicle().battery_capacity_wh +
                  ils_sp::kEnergyToleranceWh) {
        continue;
      }
      used[index] = true;
      path.push_back(stations[index]);
      const double distance =
          prefix_distance + instance.distance_km(previous, stations[index]);
      if (instance.energy_wh(stations[index], target_id) <=
          instance.vehicle().battery_capacity_wh +
              ils_sp::kEnergyToleranceWh) {
        const double last_leg =
            instance.distance_km(stations[index], target_id);
        feasible.push_back(Candidate{
            .stations = path,
            .first_leg = instance.distance_km(origin_id, path.front()),
            .last_leg = last_leg,
            .total = distance + last_leg});
      }
      self(self, stations[index], distance);
      path.pop_back();
      used[index] = false;
    }
  };
  enumerate(enumerate, origin_id, 0.0);
  std::sort(feasible.begin(), feasible.end(), [](const Candidate& left,
                                                  const Candidate& right) {
    return std::tie(left.first_leg, left.last_leg, left.total, left.stations) <
           std::tie(right.first_leg, right.last_leg, right.total,
                    right.stations);
  });
  std::vector<Candidate> nondominated;
  for (Candidate candidate : feasible) {
    const auto dominates = [&](const Candidate& left,
                               const Candidate& right) {
      const bool curves_dominate =
          std::all_of(left.stations.begin(), left.stations.end(),
                      [&](int left_station) {
                        return std::all_of(
                            right.stations.begin(), right.stations.end(),
                            [&](int right_station) {
                              return charging_curve_no_slower(
                                  instance.curve_for_station(left_station),
                                  instance.curve_for_station(right_station));
                            });
                      });
      return curves_dominate &&
             left.first_leg <= right.first_leg + ils_sp::kCostTolerance &&
             left.last_leg <= right.last_leg + ils_sp::kCostTolerance &&
             left.total <= right.total + ils_sp::kCostTolerance;
    };
    if (std::any_of(nondominated.begin(), nondominated.end(),
                    [&](const Candidate& incumbent) {
                      return dominates(incumbent, candidate);
                    })) {
      continue;
    }
    std::erase_if(nondominated, [&](const Candidate& incumbent) {
      return dominates(candidate, incumbent);
    });
    nondominated.push_back(std::move(candidate));
  }
  std::vector<std::vector<int>> result{{}};
  for (Candidate& candidate : nondominated) {
    result.push_back(std::move(candidate.stations));
  }
  std::sort(result.begin(), result.end());
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    require(argc == 2, "test executable needs the repository root");
    const std::filesystem::path root = argv[1];
    const auto instance_path =
        root / "data/evrp-2017/tc0c10s2cf1.xml";
    const ils_sp::Instance instance =
        ils_sp::load_montoya_instance(instance_path);
    require(instance.name() == "tc0c10s2cf1", "parser lost instance name");
    require(instance.customer_ids().size() == 10,
            "parser customer count differs from XML");
    require(instance.station_ids().size() == 2,
            "parser station count differs from XML");
    require_close(instance.vehicle().range_km(), 128.0, 1e-12,
                  "vehicle range conversion is wrong");
    require_close(instance.curve_for_station(11).capacity_wh(), 16'000.0,
                  1e-9, "charging curve capacity is wrong");

    ils_sp::RouteEvaluator evaluator(instance, 8);
    const ils_sp::Route direct{{7, 8}};
    const auto first = evaluator.evaluate(direct);
    const std::size_t misses = evaluator.cache_misses();
    const auto second = evaluator.evaluate(direct);
    require(first == second, "evaluation cache did not retain shared result");
    require(evaluator.cache_hits() == 1 && evaluator.cache_misses() == misses,
            "evaluation cache hit/miss counters are wrong");
    require_close(
        ils_sp::generalized_cost(*first, instance, 100.0),
        first->raw_cost_h +
            100.0 * (first->energy_shortfall_wh /
                         instance.vehicle().consumption_wh_per_km +
                     first->duration_excess_h),
        1e-12, "generalized cost does not implement Eq. (39)");

    const ils_sp::NonlinearSequenceEvaluator sequence_evaluator(instance);
    const auto require_sequence_equivalent = [&](const ils_sp::Route& route) {
      const auto full = evaluator.evaluate(route);
      const ils_sp::ChargingSequenceLabel label =
          sequence_evaluator.route_label(route);
      const ils_sp::RouteCostSummary compressed =
          sequence_evaluator.evaluate(label);
      require(compressed.feasible == full->feasible,
              "compressed nonlinear evaluation changed feasibility");
      require_close(compressed.raw_cost_h, full->raw_cost_h, 1e-12,
                    "compressed nonlinear evaluation changed raw cost");
      require_close(compressed.travel_time_h, full->travel_time_h, 1e-12,
                    "compressed nonlinear evaluation changed travel time");
      require_close(compressed.charging_time_h, full->charging_time_h, 1e-12,
                    "compressed nonlinear evaluation changed charging time");
      require_close(compressed.service_time_h, full->service_time_h, 1e-12,
                    "compressed nonlinear evaluation changed service time");
      require_close(compressed.route_duration_h, full->route_duration_h,
                    1e-12,
                    "compressed nonlinear evaluation changed duration");
      require_close(compressed.energy_shortfall_wh,
                    full->energy_shortfall_wh, 1e-8,
                    "compressed nonlinear evaluation changed shortfall");
      require_close(compressed.duration_excess_h, full->duration_excess_h,
                    1e-12,
                    "compressed nonlinear evaluation changed excess");
    };
    require_sequence_equivalent(direct);
    require_sequence_equivalent(ils_sp::Route{{7, 11, 8}});
    require_sequence_equivalent(ils_sp::Route{{1, 11, 6, 12, 10}});
    require_sequence_equivalent(
        ils_sp::Route{{1, 2, 3, 4, 5, 11, 6, 7, 8, 9, 10}});

    ils_sp::ChargingSequenceLabel concatenated;
    sequence_evaluator.append(
        concatenated,
        sequence_evaluator.edge_label(instance.depot().id, {}, 1));
    const std::vector<int> nonlinear_gap{11, 12};
    sequence_evaluator.append(
        concatenated,
        sequence_evaluator.edge_label(1, nonlinear_gap, 6));
    sequence_evaluator.append(
        concatenated,
        sequence_evaluator.edge_label(6, {}, instance.depot().id));
    const auto concatenated_cost = sequence_evaluator.evaluate(concatenated);
    const auto concatenated_full =
        evaluator.evaluate(ils_sp::Route{{1, 11, 12, 6}});
    require_close(concatenated_cost.raw_cost_h,
                  concatenated_full->raw_cost_h, 1e-12,
                  "associative sequence concatenation changed route cost");
    require_close(concatenated_cost.energy_shortfall_wh,
                  concatenated_full->energy_shortfall_wh, 1e-8,
                  "associative sequence concatenation changed shortfall");

    ils_sp::GranularNeighborhood granular(instance);
    require(granular.gamma_size() == 9,
            "Vidal gamma must contain at most 40 other customers");
    require(granular.neighbors_of(1).size() == 9,
            "directed promising-arc list has wrong size");

    ils_sp::GranularNeighborhood directed_granular(instance, 1);
    require(directed_granular.promising(1, 6) &&
                !directed_granular.promising(6, 1),
            "promising arcs were incorrectly symmetrized");
    const std::array route_six{6};
    const std::array route_one{1};
    require(directed_granular.allows_insertion(route_six, 0, 1) &&
                directed_granular.allows_insertion(route_one, 1, 6) &&
                !directed_granular.allows_insertion(route_one, 0, 6),
            "granular repair insertion ignored directed promising arcs");
    const std::vector<std::vector<int>> separated_customers{{1}, {6}};
    require(directed_granular.allows(separated_customers, {{1, 6}}) &&
                !directed_granular.allows(separated_customers, {{6, 1}}),
            "granular filter ignored promising-arc orientation");
    require(!directed_granular.allows(separated_customers,
                                      separated_customers),
            "depot-only move bypassed the promising-arc requirement");

    const ils_sp::Instance large_granular_instance =
        ils_sp::load_montoya_instance(
            root / "data/evrp-2017/tc1c320s24cf2.xml");
    require(ils_sp::GranularNeighborhood(large_granular_instance).gamma_size() ==
                40,
            "Vidal gamma must stay fixed at 40 on C320");

    ils_sp::PlanFactory factory(instance, evaluator, 16);
    auto oriented_gaps = ils_sp::direct_gaps(3);
    oriented_gaps[1] = {12};
    oriented_gaps[3] = {11};
    const ils_sp::Solution oriented_reference{{
        factory.make_plan({10, 6, 1}, std::move(oriented_gaps)),
        factory.make_plan({2, 3, 4, 5, 7, 8, 9},
                          ils_sp::direct_gaps(7)),
    }};
    oriented_reference.validate_partition(instance);
    const auto reversible_gaps = factory.gap_map(oriented_reference);
    require(reversible_gaps.at(ils_sp::AnchorPair{1, 0}) ==
                std::vector<int>({11}),
            "observed directed gap was not retained");
    require(!reversible_gaps.contains(ils_sp::AnchorPair{0, 1}),
            "gap map synthesized a reverse charging-station path");
    require(factory.plan_sequence({10, 6}, reversible_gaps).gaps[1] ==
                std::vector<int>({12}) &&
                factory.plan_sequence({6, 10}, reversible_gaps).gaps[1].empty(),
            "newly reversed gap was not materialized as direct");
    const ils_sp::RouteSequenceScorer route_scorer(
        instance, oriented_reference, reversible_gaps);
    const auto require_scored_sequence = [&](std::vector<int> sequence) {
      const ils_sp::Plan materialized =
          factory.plan_sequence(sequence, reversible_gaps);
      require_close(
          route_scorer.generalized_cost(sequence, 100.0),
          ils_sp::generalized_cost(*materialized.evaluation, instance, 100.0),
          1e-10,
          "incumbent-subsequence score differs from materialized plan");
    };
    require_scored_sequence({10, 6, 1});
    require_scored_sequence({1, 6, 10});
    require_scored_sequence({9, 8, 7, 5, 4, 3, 2, 10, 6, 1});
    require_scored_sequence({10, 2, 3, 6, 1, 4, 5, 7, 8, 9});

    // Insertion scoring must preserve the generic sequence label exactly;
    // otherwise a last-bit difference can change a deterministic tie.
    const ils_sp::Solution insertion_reference{{
        oriented_reference.plans[0],
        factory.plan_sequence({2, 3, 4, 5, 7, 8}, reversible_gaps),
    }};
    const auto insertion_reference_gaps = factory.gap_map(insertion_reference);
    const ils_sp::RouteSequenceScorer insertion_scorer(
        instance, insertion_reference, insertion_reference_gaps);
    const auto require_same_summary = [](const ils_sp::RouteCostSummary& left,
                                         const ils_sp::RouteCostSummary& right) {
      require(left.feasible == right.feasible,
              "incremental insertion changed feasibility");
      const std::array left_values{
          left.raw_cost_h,          left.travel_time_h,
          left.charging_time_h,     left.service_time_h,
          left.route_duration_h,    left.energy_shortfall_wh,
          left.duration_excess_h};
      const std::array right_values{
          right.raw_cost_h,          right.travel_time_h,
          right.charging_time_h,     right.service_time_h,
          right.route_duration_h,    right.energy_shortfall_wh,
          right.duration_excess_h};
      for (std::size_t index = 0; index < left_values.size(); ++index) {
        require(std::bit_cast<std::uint64_t>(left_values[index]) ==
                    std::bit_cast<std::uint64_t>(right_values[index]),
                "incremental insertion changed a summary bit");
      }
    };
    for (std::size_t route = 0; route < insertion_reference.plans.size();
         ++route) {
      for (std::size_t position = 0;
           position <= insertion_reference.plans[route].customer_ids.size();
           ++position) {
        std::vector<int> materialized =
            insertion_reference.plans[route].customer_ids;
        materialized.insert(
            materialized.begin() + static_cast<std::ptrdiff_t>(position), 9);
        const auto generic = insertion_scorer.evaluate(materialized);
        const auto incremental =
            insertion_scorer.evaluate_insertion(route, position, 9);
        require_same_summary(incremental, generic);
        require(std::bit_cast<std::uint64_t>(ils_sp::generalized_cost(
                    incremental, instance, 100.0)) ==
                    std::bit_cast<std::uint64_t>(ils_sp::generalized_cost(
                        generic, instance, 100.0)),
                "incremental insertion changed generalized-cost bits");
      }
    }

    const auto require_scored_slices =
        [&](std::span<const ils_sp::IncumbentRouteSlice> slices,
            const std::vector<int>& sequence) {
          require_close(
              route_scorer.generalized_cost(slices, 100.0),
              route_scorer.generalized_cost(sequence, 100.0), 1e-10,
              "descriptor slice score differs from sequence score");
        };
    const std::array reverse_first{
        ils_sp::IncumbentRouteSlice{0, 0, 3, true}};
    require_scored_slices(reverse_first, {1, 6, 10});
    const std::array reverse_then_forward{
        ils_sp::IncumbentRouteSlice{1, 0, 7, true},
        ils_sp::IncumbentRouteSlice{0, 0, 3, false}};
    require_scored_slices(reverse_then_forward,
                          {9, 8, 7, 5, 4, 3, 2, 10, 6, 1});
    const std::array four_fragments{
        ils_sp::IncumbentRouteSlice{0, 0, 1, false},
        ils_sp::IncumbentRouteSlice{1, 0, 2, false},
        ils_sp::IncumbentRouteSlice{0, 1, 3, false},
        ils_sp::IncumbentRouteSlice{1, 2, 7, false}};
    require_scored_slices(four_fragments,
                          {10, 2, 3, 6, 1, 4, 5, 7, 8, 9});

    ils_sp::RouteEvaluator exact_scope_evaluator(instance, 16);
    ils_sp::PlanFactory exact_scope_factory(instance, exact_scope_evaluator,
                                            16);
    ils_sp::Plan cached_exact = exact_scope_factory.make_plan(
        {7, 8}, ils_sp::direct_gaps(2));
    require(cached_exact.evaluation->feasible,
            "exact-cache scope fixture must be feasible");
    const auto cached_heuristic_evaluation = cached_exact.evaluation;
    cached_exact.evaluation = exact_scope_evaluator.evaluate_schedule(
        cached_exact.route,
        std::vector<std::optional<double>>(cached_exact.route.visits.size()));
    cached_exact.exact_charging = true;
    exact_scope_factory.publish_exact(cached_exact);
    const ils_sp::Plan same_physical_heuristic =
        exact_scope_factory.plan_sequence({7, 8});
    require(same_physical_heuristic.evaluation == cached_heuristic_evaluation &&
                same_physical_heuristic.evaluation != cached_exact.evaluation,
            "FRVCP replay polluted the Algorithm-3 route cache");
    auto heuristic_gaps = ils_sp::direct_gaps(2);
    heuristic_gaps[0] = {11};
    const ils_sp::Plan heuristic =
        exact_scope_factory.plan_sequence({7, 8},
                                          ils_sp::PlanFactory::GapMap{
                                              {ils_sp::AnchorPair{0, 7},
                                               heuristic_gaps[0]}});
    require(!heuristic.exact_charging && heuristic.gaps[0] ==
                                             std::vector<int>({11}),
            "FRVCP cache leaked into ordinary ILS plan construction");
    require(exact_scope_factory.exact_plan({7, 8}).has_value(),
            "explicit FRVCP cache lookup lost the cached exact plan");

    const ils_sp::Plan forward =
        factory.make_plan({7, 8}, ils_sp::direct_gaps(2));
    const ils_sp::Plan reverse =
        factory.make_plan({8, 7}, ils_sp::direct_gaps(2));
    require(forward.evaluation->feasible && reverse.evaluation->feasible,
            "pool dominance fixture routes must be feasible");
    ils_sp::RoutePool pool(4);
    (void)pool.admit(ils_sp::RouteColumn::from_plan(forward));
    (void)pool.admit(ils_sp::RouteColumn::from_plan(reverse));
    require(pool.main_size() == 2,
            "different ordered customer sequences were incorrectly dominated");
    require(pool.distinct_main_coverages() == 1,
            "route-pool distinct coverage statistic is incorrect");
    (void)pool.admit(ils_sp::RouteColumn::from_plan(forward));
    require(pool.main_size() == 2,
            "the same ordered customer sequence was duplicated");
    ils_sp::Plan exact_forward = forward;
    exact_forward.exact_charging = true;
    (void)pool.admit(ils_sp::RouteColumn::from_plan(exact_forward));
    const auto ordered_columns = pool.main_columns();
    const auto upgraded = std::find_if(
        ordered_columns.begin(), ordered_columns.end(),
        [](const ils_sp::RouteColumn& column) {
          return column.plan.customer_ids == std::vector<int>({7, 8});
        });
    require(upgraded != ordered_columns.end() &&
                upgraded->plan.exact_charging,
            "equal-cost exact route did not upgrade the same ordered sequence");

    ils_sp::RoutePool coverage_pool(4);
    ils_sp::RouteColumn coverage_forward =
        ils_sp::RouteColumn::from_plan(forward);
    ils_sp::RouteColumn coverage_reverse =
        ils_sp::RouteColumn::from_plan(reverse);
    coverage_forward.cost_h = 2.0;
    coverage_reverse.cost_h = 3.0;
    require(coverage_pool.admit(coverage_forward)
                .effective_coverage_improvement,
            "first customer coverage was not reported as effective");
    require(!coverage_pool.admit(coverage_reverse)
                 .effective_coverage_improvement,
            "a more expensive ordering incorrectly refreshed coverage");
    coverage_reverse.cost_h = 1.0;
    require(coverage_pool.admit(std::move(coverage_reverse))
                .effective_coverage_improvement,
            "a cheaper ordering did not refresh coverage quality");
    require(coverage_pool.distinct_main_coverages() == 1,
            "coverage quality tracking changed coverage identity");

    ils_sp::RouteColumn bounded_forward =
        ils_sp::RouteColumn::from_plan(forward);
    ils_sp::RouteColumn bounded_reverse =
        ils_sp::RouteColumn::from_plan(reverse);
    // Control the rank explicitly: the test exercises ordered-sequence map
    // erasure during repricing, not floating-point symmetry of the fixture.
    bounded_forward.cost_h = 1.0;
    bounded_reverse.cost_h = 2.0;
    ils_sp::RoutePool bounded_pool(1);
    (void)bounded_pool.admit(std::move(bounded_forward));
    (void)bounded_pool.admit(std::move(bounded_reverse));
    require(!bounded_pool.duals_valid() && bounded_pool.main_size() == 1 &&
                bounded_pool.backup_size() == 1,
            "pre-dual raw-cost bootstrap did not keep the bounded pool");
    bounded_pool.update_duals({{7, 10.0}, {8, 10.0}});
    const ils_sp::PoolRepriceUpdate bounded_reprice =
        bounded_pool.reprice_backup();
    require(bounded_reprice.promoted == 0 && bounded_pool.main_size() == 1 &&
                bounded_pool.backup_size() == 1 &&
                bounded_pool.main_columns().front().plan.customer_ids ==
                    forward.customer_ids,
            "a self-evicted repricing candidate was not recoverable");

    ils_sp::RoutePool recovering_pool(1);
    ils_sp::RouteColumn recovering_single =
        ils_sp::RouteColumn::from_plan(
            factory.make_plan({7}, ils_sp::direct_gaps(1)));
    ils_sp::RouteColumn recovering_combined =
        ils_sp::RouteColumn::from_plan(forward);
    recovering_single.cost_h = 1.0;
    recovering_combined.cost_h = 2.0;
    (void)recovering_pool.admit(recovering_single);
    (void)recovering_pool.admit(recovering_combined);
    require(recovering_pool.main_columns().front().plan.customer_ids ==
                    std::vector<int>({7}) &&
                recovering_pool.backup_size() == 1,
            "pre-dual raw-cost ordering admitted the worse column");
    recovering_pool.update_duals({{8, 10.0}});
    const ils_sp::PoolRepriceUpdate combined_promoted =
        recovering_pool.reprice_backup();
    require(combined_promoted.promoted == 1 &&
                recovering_pool.main_columns().front().plan.customer_ids ==
                    forward.customer_ids &&
                recovering_pool.backup_size() == 1,
            "a dual change did not recover the useful backup column");
    recovering_pool.update_duals({{7, 10.0}});
    const ils_sp::PoolRepriceUpdate single_recovered =
        recovering_pool.reprice_backup();
    require(single_recovered.promoted == 1 &&
                recovering_pool.main_columns().front().plan.customer_ids ==
                    std::vector<int>({7}) &&
                recovering_pool.backup_size() == 1 &&
                recovering_pool.statistics().main_demotions >= 2,
            "an evicted main column could not recover after the duals changed");

    ils_sp::RoutePool protected_pool(2);
    const ils_sp::Plan protected_seven =
        factory.make_plan({7}, ils_sp::direct_gaps(1));
    const ils_sp::Plan protected_eight =
        factory.make_plan({8}, ils_sp::direct_gaps(1));
    protected_pool.protect_partition(
        ils_sp::Solution{{protected_seven, protected_eight}});
    protected_pool.update_duals({{7, 10.0}, {8, 10.0}});
    (void)protected_pool.admit(ils_sp::RouteColumn::from_plan(forward));
    const auto protected_columns = protected_pool.main_columns();
    require(protected_columns.size() == 2 &&
                std::any_of(protected_columns.begin(), protected_columns.end(),
                            [](const ils_sp::RouteColumn& column) {
                              return column.plan.customer_ids ==
                                     std::vector<int>({7});
                            }) &&
                std::any_of(protected_columns.begin(), protected_columns.end(),
                            [](const ils_sp::RouteColumn& column) {
                              return column.plan.customer_ids ==
                                     std::vector<int>({8});
                            }),
            "best-feasible protected partition was evicted from the pool");

    std::mt19937_64 random(36'509);
    ils_sp::PathSampler insertion_paths(instance, factory, random);
    const ils_sp::Plan& replacement_fixture =
        oriented_reference.plans.front();
    std::vector<int> replacement_anchors{instance.depot().id};
    replacement_anchors.insert(replacement_anchors.end(),
                               replacement_fixture.customer_ids.begin(),
                               replacement_fixture.customer_ids.end());
    replacement_anchors.push_back(instance.depot().id);
    const double replacement_incumbent = ils_sp::generalized_cost(
        *replacement_fixture.evaluation, instance, 100.0);
    struct RankedReplacement {
      double cost{};
      ils_sp::Plan plan;
    };
    double replacement_minimum_cost =
        std::numeric_limits<double>::infinity();
    std::vector<RankedReplacement> replacement_candidates;
    for (std::size_t gap = 0; gap < replacement_fixture.gaps.size(); ++gap) {
      for (const auto& option : insertion_paths.paths_between(
               replacement_anchors[gap], replacement_anchors[gap + 1])) {
        if (option == replacement_fixture.gaps[gap]) continue;
        auto changed_gaps = replacement_fixture.gaps;
        changed_gaps[gap] = option;
        ils_sp::Plan candidate = factory.make_plan(
            replacement_fixture.customer_ids, std::move(changed_gaps));
        const double candidate_cost = ils_sp::generalized_cost(
            *candidate.evaluation, instance, 100.0);
        if (candidate_cost <
            replacement_incumbent - ils_sp::kCostTolerance) {
          replacement_minimum_cost =
              std::min(replacement_minimum_cost, candidate_cost);
          replacement_candidates.push_back(
              RankedReplacement{candidate_cost, std::move(candidate)});
        }
      }
    }
    std::optional<ils_sp::Plan> replacement_expected;
    double replacement_best_cost = replacement_incumbent;
    for (RankedReplacement& candidate : replacement_candidates) {
      if (candidate.cost >
          replacement_minimum_cost + ils_sp::kCostTolerance) {
        continue;
      }
      if (!replacement_expected.has_value() ||
          candidate.plan.route.visits < replacement_expected->route.visits) {
        replacement_best_cost = candidate.cost;
        replacement_expected = std::move(candidate.plan);
      }
    }
    const auto replacement_actual =
        insertion_paths.best_replacement(replacement_fixture, 100.0);
    require(replacement_actual.has_value() ==
                replacement_expected.has_value(),
            "incremental ReplacePath ranking changed improvement existence");
    if (replacement_actual.has_value()) {
      require(replacement_actual->route.visits ==
                  replacement_expected->route.visits,
              "incremental ReplacePath ranking changed the best route");
      require_close(ils_sp::generalized_cost(*replacement_actual->evaluation,
                                             instance, 100.0),
                    replacement_best_cost, 1e-10,
                    "incremental ReplacePath ranking changed best cost");
    }
    std::vector<int> anchors{instance.depot().id};
    anchors.insert(anchors.end(), instance.customer_ids().begin(),
                   instance.customer_ids().end());
    for (const int origin : anchors) {
      for (const int target : anchors) {
        if (origin == target) continue;
        require(insertion_paths.paths_between(origin, target) ==
                    exhaustive_station_catalog(instance, origin, target),
                "CS-path catalog omitted or retained a dominated path");
      }
    }
    const ils_sp::Instance four_station_instance =
        ils_sp::load_montoya_instance(
            root / "data/evrp-2017/tc0c20s4cf2.xml");
    require(four_station_instance.station_ids().size() == 4,
            "four-station oracle fixture has the wrong station count");
    ils_sp::RouteEvaluator four_station_evaluator(four_station_instance, 8);
    ils_sp::PlanFactory four_station_factory(four_station_instance,
                                             four_station_evaluator, 16);
    std::mt19937_64 four_station_random(36'509);
    ils_sp::PathSampler four_station_paths(
        four_station_instance, four_station_factory, four_station_random);
    std::vector<int> four_station_anchors{four_station_instance.depot().id};
    four_station_anchors.insert(four_station_anchors.end(),
                                four_station_instance.customer_ids().begin(),
                                four_station_instance.customer_ids().end());
    for (const int origin : four_station_anchors) {
      for (const int target : four_station_anchors) {
        if (origin == target) continue;
        require(four_station_paths.paths_between(origin, target) ==
                    exhaustive_station_catalog(four_station_instance, origin,
                                               target),
                "shortest-path reduction changed the nondominated catalog");
      }
    }
    const ils_sp::Instance mixed_curve_instance =
        ils_sp::load_montoya_instance(
            root / "data/evrp-2017/tc1c10s2ct2.xml");
    require(mixed_curve_instance.node(11).station_type !=
                mixed_curve_instance.node(12).station_type,
            "mixed-curve path fixture lost its charging technologies");
    ils_sp::RouteEvaluator mixed_curve_evaluator(mixed_curve_instance, 8);
    ils_sp::PlanFactory mixed_curve_factory(mixed_curve_instance,
                                            mixed_curve_evaluator, 16);
    std::mt19937_64 mixed_curve_random(36'509);
    ils_sp::PathSampler mixed_curve_paths(
        mixed_curve_instance, mixed_curve_factory, mixed_curve_random);
    const ils_sp::Route mixed_curve_route{{1, 11, 6, 12, 10}};
    const auto mixed_curve_full =
        mixed_curve_evaluator.evaluate(mixed_curve_route);
    const ils_sp::NonlinearSequenceEvaluator mixed_sequence_evaluator(
        mixed_curve_instance);
    const auto mixed_curve_compressed = mixed_sequence_evaluator.evaluate(
        mixed_sequence_evaluator.route_label(mixed_curve_route));
    require_close(mixed_curve_compressed.raw_cost_h,
                  mixed_curve_full->raw_cost_h, 1e-12,
                  "mixed-technology sequence label changed raw cost");
    require_close(mixed_curve_compressed.energy_shortfall_wh,
                  mixed_curve_full->energy_shortfall_wh, 1e-8,
                  "mixed-technology sequence label changed shortfall");
    std::mt19937_64 sequence_equivalence_random(91'337);
    for (std::size_t trial = 0; trial < 64; ++trial) {
      std::vector<int> customers = mixed_curve_instance.customer_ids();
      std::shuffle(customers.begin(), customers.end(),
                   sequence_equivalence_random);
      customers.resize(1 + trial % customers.size());
      std::vector<int> visits;
      for (const int customer_id : customers) {
        if (sequence_equivalence_random() % 4 == 0) {
          visits.push_back(mixed_curve_instance.station_ids()[
              sequence_equivalence_random() %
              mixed_curve_instance.station_ids().size()]);
        }
        visits.push_back(customer_id);
        if (sequence_equivalence_random() % 4 == 0) {
          visits.push_back(mixed_curve_instance.station_ids()[
              sequence_equivalence_random() %
              mixed_curve_instance.station_ids().size()]);
        }
      }
      const ils_sp::Route sampled_route{std::move(visits)};
      const auto sampled_full =
          mixed_curve_evaluator.evaluate(sampled_route);
      const auto sampled_compressed = mixed_sequence_evaluator.evaluate(
          mixed_sequence_evaluator.route_label(sampled_route));
      require(sampled_compressed.feasible == sampled_full->feasible,
              "sampled sequence label changed feasibility");
      require_close(sampled_compressed.raw_cost_h, sampled_full->raw_cost_h,
                    1e-11, "sampled sequence label changed raw cost");
      require_close(sampled_compressed.energy_shortfall_wh,
                    sampled_full->energy_shortfall_wh, 1e-7,
                    "sampled sequence label changed shortfall");
      require_close(sampled_compressed.duration_excess_h,
                    sampled_full->duration_excess_h, 1e-11,
                    "sampled sequence label changed duration excess");
    }
    require(mixed_curve_paths.paths_between(
                10, mixed_curve_instance.depot().id) ==
                exhaustive_station_catalog(mixed_curve_instance, 10,
                                           mixed_curve_instance.depot().id),
            "curve-blind dominance removed a mixed-technology path");

    // A cached move stores absolute route slots.  This deterministic ILS
    // trajectory makes a historical route reappear in a different slot; the
    // cache must not apply that route's old descriptor to the new occupant.
    std::mt19937_64 move_memory_random(36'509);
    ils_sp::RouteEvaluator move_memory_evaluator(mixed_curve_instance,
                                                 200'000);
    ils_sp::PlanFactory move_memory_factory(mixed_curve_instance,
                                            move_memory_evaluator, 20'000);
    ils_sp::PathSampler move_memory_paths(
        mixed_curve_instance, move_memory_factory, move_memory_random);
    ils_sp::InitialSolutionBuilder move_memory_initial(
        mixed_curve_instance, move_memory_factory, move_memory_paths,
        move_memory_random);
    ils_sp::Perturbation move_memory_perturbation(
        mixed_curve_instance, move_memory_factory, move_memory_paths,
        move_memory_random);
    ils_sp::VariableNeighborhoodDescent move_memory_vnd(
        mixed_curve_instance, move_memory_factory, move_memory_paths,
        move_memory_random, 40);
    ils_sp::Solution move_memory_solution = move_memory_initial.build(100.0);
    for (std::size_t iteration = 1; iteration <= 42; ++iteration) {
      move_memory_solution = move_memory_vnd.improve(
          move_memory_perturbation.apply(move_memory_solution, 100.0), 100.0);
      move_memory_solution.validate_partition(mixed_curve_instance);
    }
    require(move_memory_vnd.statistics().route_support_cache_hits > 0,
            "VND route-support move memory was not exercised");
    require(move_memory_vnd.statistics().relocate_source_summary_hits > 0,
            "Relocate source-route summary reuse was not exercised");

    ils_sp::InitialSolutionBuilder initial(instance, factory, insertion_paths,
                                           random);
    const ils_sp::Solution initial_solution = initial.build(100.0);
    initial_solution.validate_partition(instance);
    require(!initial_solution.plans.empty(), "initial solution is empty");
    require(std::all_of(
                initial_solution.plans.begin(), initial_solution.plans.end(),
                [](const ils_sp::Plan& plan) {
                  return std::all_of(
                      plan.gaps.begin(), plan.gaps.end(),
                      [](const std::vector<int>& gap) { return gap.empty(); });
                }),
            "customer insertion introduced a charging-station path");
    std::mt19937_64 default_profile_random(9'876);
    std::mt19937_64 explicit_paper_random(9'876);
    ils_sp::RouteEvaluator default_profile_evaluator(instance, 256);
    ils_sp::RouteEvaluator explicit_paper_evaluator(instance, 256);
    ils_sp::PlanFactory default_profile_factory(
        instance, default_profile_evaluator, 64);
    ils_sp::PlanFactory explicit_paper_factory(
        instance, explicit_paper_evaluator, 64);
    ils_sp::PathSampler default_profile_paths(
        instance, default_profile_factory, default_profile_random);
    ils_sp::PathSampler explicit_paper_paths(
        instance, explicit_paper_factory, explicit_paper_random);
    ils_sp::InitialSolutionBuilder default_profile_initial(
        instance, default_profile_factory, default_profile_paths,
        default_profile_random);
    ils_sp::InitialSolutionBuilder explicit_paper_initial(
        instance, explicit_paper_factory, explicit_paper_paths,
        explicit_paper_random, ils_sp::SearchProfile::Paper);
    const ils_sp::Solution default_profile_solution =
        default_profile_initial.build(100.0);
    const ils_sp::Solution explicit_paper_solution =
        explicit_paper_initial.build(100.0);
    require(default_profile_solution.customer_sequences() ==
                    explicit_paper_solution.customer_sequences() &&
                default_profile_factory.gap_map(default_profile_solution) ==
                    explicit_paper_factory.gap_map(
                        explicit_paper_solution) &&
                default_profile_random() == explicit_paper_random(),
            "default and explicit paper profiles changed solution or RNG");
    require(default_profile_initial.statistics().candidates_generated == 0 &&
                explicit_paper_initial.statistics().candidates_generated == 0,
            "paper initialization entered candidate-linked scoring");
    std::mt19937_64 linked_initial_random(36'509);
    ils_sp::InitialSolutionBuilder linked_initial(
        instance, factory, insertion_paths, linked_initial_random,
        ils_sp::SearchProfile::LinkedWinner);
    const ils_sp::Solution linked_initial_solution =
        linked_initial.build(100.0);
    linked_initial_solution.validate_partition(instance);
    require(std::any_of(
                linked_initial_solution.plans.begin(),
                linked_initial_solution.plans.end(),
                [](const ils_sp::Plan& plan) {
                  return std::any_of(
                      plan.gaps.begin(), plan.gaps.end(),
                      [](const std::vector<int>& gap) { return !gap.empty(); });
                }),
            "linked initial construction did not complete any new gap");

    // The experimental candidate-linked constructor must choose the true
    // best completed endpoint, not merely complete the direct-cost winner.
    // Compare its branch-and-bound result with a small exhaustive oracle.
    std::mt19937_64 candidate_initial_random(36'509);
    std::mt19937_64 exhaustive_initial_random(36'509);
    ils_sp::RouteEvaluator candidate_initial_evaluator(instance, 4'096);
    ils_sp::RouteEvaluator exhaustive_initial_evaluator(instance, 4'096);
    ils_sp::PlanFactory candidate_initial_factory(
        instance, candidate_initial_evaluator, 256);
    ils_sp::PlanFactory exhaustive_initial_factory(
        instance, exhaustive_initial_evaluator, 256);
    ils_sp::PathSampler candidate_initial_paths(
        instance, candidate_initial_factory, candidate_initial_random);
    ils_sp::PathSampler exhaustive_initial_paths(
        instance, exhaustive_initial_factory, exhaustive_initial_random);
    ils_sp::InitialSolutionBuilder candidate_initial(
        instance, candidate_initial_factory, candidate_initial_paths,
        candidate_initial_random,
        ils_sp::SearchProfile::LinkedInsertionCandidates);
    const ils_sp::Solution candidate_initial_solution =
        candidate_initial.build(100.0);

    std::vector<int> exhaustive_unassigned = instance.customer_ids();
    const std::size_t exhaustive_first_index =
        std::uniform_int_distribution<std::size_t>(
            0, exhaustive_unassigned.size() - 1)(exhaustive_initial_random);
    const int exhaustive_first =
        exhaustive_unassigned[exhaustive_first_index];
    exhaustive_unassigned.erase(
        exhaustive_unassigned.begin() +
        static_cast<std::ptrdiff_t>(exhaustive_first_index));
    ils_sp::Plan exhaustive_first_plan =
        exhaustive_initial_factory.plan_sequence({exhaustive_first});
    exhaustive_first_plan = exhaustive_initial_paths.repair_new_gaps(
        std::move(exhaustive_first_plan), ils_sp::PlanFactory::GapMap{},
        100.0);
    ils_sp::Solution exhaustive_initial_solution{{
        std::move(exhaustive_first_plan),
    }};
    while (!exhaustive_unassigned.empty()) {
      const auto inherited =
          exhaustive_initial_factory.gap_map(exhaustive_initial_solution);
      struct ExhaustiveInsertion {
        double delta{};
        int customer{};
        std::size_t route{};
        std::size_t position{};
        ils_sp::Plan plan;
      };
      std::optional<ExhaustiveInsertion> best_insertion;
      for (const int customer_id : exhaustive_unassigned) {
        for (std::size_t route = 0;
             route <= exhaustive_initial_solution.plans.size(); ++route) {
          const std::size_t positions =
              route == exhaustive_initial_solution.plans.size()
                  ? 1
                  : exhaustive_initial_solution.plans[route]
                            .customer_ids.size() +
                        1;
          for (std::size_t position = 0; position < positions; ++position) {
            std::vector<int> sequence;
            double incumbent_cost = 0.0;
            if (route == exhaustive_initial_solution.plans.size()) {
              sequence = {customer_id};
            } else {
              sequence = exhaustive_initial_solution.plans[route].customer_ids;
              sequence.insert(
                  sequence.begin() + static_cast<std::ptrdiff_t>(position),
                  customer_id);
              incumbent_cost = ils_sp::generalized_cost(
                  *exhaustive_initial_solution.plans[route].evaluation,
                  instance, 100.0);
            }
            ils_sp::Plan plan =
                exhaustive_initial_factory.plan_sequence(sequence, inherited);
            plan = exhaustive_initial_paths.repair_new_gaps(
                std::move(plan), inherited, 100.0);
            const double delta = ils_sp::generalized_cost(
                                     *plan.evaluation, instance, 100.0) -
                                 incumbent_cost;
            ExhaustiveInsertion candidate{.delta = delta,
                                           .customer = customer_id,
                                           .route = route,
                                           .position = position,
                                           .plan = std::move(plan)};
            if (!best_insertion.has_value() ||
                std::tie(candidate.delta, candidate.customer,
                         candidate.route, candidate.position) <
                    std::tie(best_insertion->delta,
                             best_insertion->customer,
                             best_insertion->route,
                             best_insertion->position)) {
              best_insertion = std::move(candidate);
            }
          }
        }
      }
      const int inserted = best_insertion->customer;
      if (best_insertion->route == exhaustive_initial_solution.plans.size()) {
        exhaustive_initial_solution.plans.push_back(
            std::move(best_insertion->plan));
      } else {
        exhaustive_initial_solution.plans[best_insertion->route] =
            std::move(best_insertion->plan);
      }
      std::erase(exhaustive_unassigned, inserted);
    }
    candidate_initial_solution.validate_partition(instance);
    exhaustive_initial_solution.validate_partition(instance);
    const auto& candidate_initial_stats = candidate_initial.statistics();
    require(candidate_initial_solution.customer_sequences() ==
                    exhaustive_initial_solution.customer_sequences() &&
                candidate_initial_factory.gap_map(candidate_initial_solution) ==
                    exhaustive_initial_factory.gap_map(
                        exhaustive_initial_solution) &&
                candidate_initial_random() == exhaustive_initial_random(),
            "candidate-linked construction differs from exhaustive completion");
    require(candidate_initial_stats.decisions == 9 &&
                candidate_initial_stats.candidates_evaluated > 0 &&
                candidate_initial_stats.completion_memo_hits +
                        candidate_initial_stats.completion_memo_misses ==
                    candidate_initial_stats.candidates_evaluated &&
                candidate_initial_stats.winner_flips > 0 &&
                candidate_initial_stats.strict_winner_flips > 0 &&
                candidate_initial_stats.candidates_generated ==
                    candidate_initial_stats.candidates_evaluated +
                        candidate_initial_stats.candidates_pruned_by_bound +
                        candidate_initial_stats.candidates_pruned_by_regret,
            "candidate-linked initialization statistics are inconsistent");

    const std::vector<int> repair_removed{1, 4, 7};
    for (const ils_sp::RepairOperator operation :
         {ils_sp::RepairOperator::GreedyInsertion,
          ils_sp::RepairOperator::TwoRegretInsertion}) {
      std::mt19937_64 candidate_repair_random(
          operation == ils_sp::RepairOperator::GreedyInsertion ? 12'345
                                                                : 12'346);
      std::mt19937_64 exhaustive_repair_random(
          operation == ils_sp::RepairOperator::GreedyInsertion ? 12'345
                                                                : 12'346);
      ils_sp::RouteEvaluator candidate_repair_evaluator(instance, 4'096);
      ils_sp::RouteEvaluator exhaustive_repair_evaluator(instance, 4'096);
      ils_sp::PlanFactory candidate_repair_factory(
          instance, candidate_repair_evaluator, 256);
      ils_sp::PlanFactory exhaustive_repair_factory(
          instance, exhaustive_repair_evaluator, 256);
      ils_sp::PathSampler candidate_repair_paths(
          instance, candidate_repair_factory, candidate_repair_random);
      ils_sp::PathSampler exhaustive_repair_paths(
          instance, exhaustive_repair_factory, exhaustive_repair_random);
      ils_sp::Perturbation candidate_repair(
          instance, candidate_repair_factory, candidate_repair_paths,
          candidate_repair_random,
          ils_sp::SearchProfile::LinkedInsertionCandidates);
      const ils_sp::Solution candidate_repaired =
          ils_sp::PerturbationTestAccess::repair(
              candidate_repair,
              make_destroyed_fixture(candidate_initial_solution,
                                     candidate_repair_factory,
                                     repair_removed),
              operation, 100.0);
      const ils_sp::Solution exhaustive_repaired = exhaustive_linked_repair(
          make_destroyed_fixture(candidate_initial_solution,
                                 exhaustive_repair_factory, repair_removed),
          operation, instance, exhaustive_repair_factory,
          exhaustive_repair_paths, 100.0);
      require(candidate_repaired.customer_sequences() ==
                      exhaustive_repaired.customer_sequences() &&
                  candidate_repair_factory.gap_map(candidate_repaired) ==
                      exhaustive_repair_factory.gap_map(
                          exhaustive_repaired),
              "candidate-linked repair differs from exhaustive completion");
      const auto& statistics = candidate_repair.statistics();
      require(statistics.decisions == repair_removed.size() &&
                  statistics.candidates_evaluated > 0 &&
                  statistics.completion_memo_hits +
                          statistics.completion_memo_misses ==
                      statistics.candidates_evaluated &&
                  statistics.candidates_generated ==
                      statistics.candidates_evaluated +
                          statistics.candidates_pruned_by_bound +
                          statistics.candidates_pruned_by_regret &&
                  (operation == ils_sp::RepairOperator::GreedyInsertion
                       ? statistics.greedy_decisions == repair_removed.size()
                       : statistics.two_regret_decisions ==
                             repair_removed.size()),
              "candidate-linked repair statistics are inconsistent");
      if (operation == ils_sp::RepairOperator::TwoRegretInsertion) {
        require(statistics.completion_memo_hits > 0 &&
                    statistics.customers_pruned_by_regret > 0 &&
                    statistics.candidates_pruned_by_regret > 0,
                "2-Regret insertion pruning was not exercised");
      }
    }

    std::mt19937_64 direct_policy_random(12'345);
    ils_sp::Perturbation direct_policy_perturbation(
        instance, factory, insertion_paths, direct_policy_random);
    const ils_sp::Solution direct_perturbed =
        direct_policy_perturbation.apply(initial_solution, 100.0);
    require(std::all_of(
                direct_perturbed.plans.begin(), direct_perturbed.plans.end(),
                [](const ils_sp::Plan& plan) {
                  return std::all_of(
                      plan.gaps.begin(), plan.gaps.end(),
                      [](const std::vector<int>& gap) { return gap.empty(); });
                }),
            "destroy/repair insertion introduced a charging-station path");

    std::mt19937_64 linked_policy_random(12'345);
    ils_sp::Perturbation linked_policy_perturbation(
        instance, factory, insertion_paths, linked_policy_random,
        ils_sp::SearchProfile::LinkedWinner);
    const ils_sp::Solution linked_perturbed =
        linked_policy_perturbation.apply(initial_solution, 100.0);
    linked_perturbed.validate_partition(instance);
    require(std::any_of(
                linked_perturbed.plans.begin(), linked_perturbed.plans.end(),
                [](const ils_sp::Plan& plan) {
                  return std::any_of(
                      plan.gaps.begin(), plan.gaps.end(),
                      [](const std::vector<int>& gap) { return !gap.empty(); });
                }),
            "linked repair did not complete any selected insertion gap");

    std::mt19937_64 granular_random(36'509);
    ils_sp::RouteEvaluator granular_evaluator(instance, 4'096);
    ils_sp::PlanFactory granular_factory(instance, granular_evaluator, 256);
    ils_sp::PathSampler granular_paths(instance, granular_factory,
                                       granular_random);
    ils_sp::InitialSolutionBuilder granular_initial(
        instance, granular_factory, granular_paths, granular_random,
        ils_sp::SearchProfile::GranularAdaptive);
    const ils_sp::Solution granular_initial_solution =
        granular_initial.build(100.0);
    require(
        granular_initial.statistics().candidates_rejected_by_granular == 0 &&
            granular_initial.statistics().virtual_candidates_evaluated > 0 &&
            granular_initial.statistics().virtual_strict_winner_flips > 0 &&
            std::all_of(
                granular_initial_solution.plans.begin(),
                granular_initial_solution.plans.end(),
                [](const ils_sp::Plan& plan) {
                  return std::all_of(
                      plan.gaps.begin(), plan.gaps.end(),
                      [](const std::vector<int>& gap) { return gap.empty(); });
                }),
            "lazy path-aware initialization did not reverse a candidate "
            "ranking or wrote its virtual path into the repair result");

    require(granular_initial_solution.plans.size() >= 2,
            "virtual completion cache fixture needs two routes");
    std::mt19937_64 virtual_cache_random(91'001);
    ils_sp::RouteEvaluator virtual_cache_evaluator(instance, 512);
    ils_sp::PlanFactory virtual_cache_factory(instance,
                                              virtual_cache_evaluator, 64);
    ils_sp::PathSampler virtual_cache_paths(instance, virtual_cache_factory,
                                            virtual_cache_random);
    ils_sp::RepairTimingStatistics virtual_cache_timing;
    const ils_sp::Plan& virtual_cache_source =
        granular_initial_solution.plans.front();
    const int virtual_cache_customer =
        granular_initial_solution.plans.back().customer_ids.front();
    const auto virtual_cache_first =
        virtual_cache_paths.virtual_complete_insertion(
            &virtual_cache_source, 0, virtual_cache_customer, 100.0,
            &virtual_cache_timing);
    const auto virtual_cache_second =
        virtual_cache_paths.virtual_complete_insertion(
            &virtual_cache_source, 0, virtual_cache_customer, 100.0,
            &virtual_cache_timing);
    (void)virtual_cache_paths.virtual_complete_insertion(
        &virtual_cache_source, 0, virtual_cache_customer, 10'000.0,
        &virtual_cache_timing);
    require(virtual_cache_first.generalized_cost_h ==
                    virtual_cache_second.generalized_cost_h &&
                virtual_cache_first.paths == virtual_cache_second.paths &&
                virtual_cache_timing.virtual_completion_cross_cache_hits == 1 &&
                virtual_cache_timing.virtual_completion_cross_cache_misses ==
                    2 &&
                virtual_cache_paths.virtual_completion_cache_size() == 2,
            "cross-repair virtual completion cache ignored physical-route or "
            "exact-lambda identity");

    for (const ils_sp::RepairOperator operation :
         {ils_sp::RepairOperator::GreedyInsertion,
          ils_sp::RepairOperator::TwoRegretInsertion}) {
      std::mt19937_64 certified_random(71'000 +
                                       static_cast<unsigned>(operation));
      std::mt19937_64 virtual_oracle_random(
          71'000 + static_cast<unsigned>(operation));
      ils_sp::RouteEvaluator certified_evaluator(instance, 4'096);
      ils_sp::RouteEvaluator virtual_oracle_evaluator(instance, 4'096);
      ils_sp::PlanFactory certified_factory(instance, certified_evaluator,
                                            256);
      ils_sp::PlanFactory virtual_oracle_factory(
          instance, virtual_oracle_evaluator, 256);
      ils_sp::PathSampler certified_paths(instance, certified_factory,
                                          certified_random);
      ils_sp::PathSampler virtual_oracle_paths(
          instance, virtual_oracle_factory, virtual_oracle_random);
      ils_sp::Perturbation certified_repair(
          instance, certified_factory, certified_paths, certified_random,
          ils_sp::SearchProfile::GranularAdaptive);
      const std::vector<int> removed{instance.customer_ids()[0],
                                     instance.customer_ids()[1],
                                     instance.customer_ids()[2]};
      const ils_sp::Solution certified =
          ils_sp::PerturbationTestAccess::repair(
              certified_repair,
              make_destroyed_fixture(granular_initial_solution,
                                     certified_factory, removed),
              operation, 100.0);
      const ils_sp::Solution oracle = exhaustive_virtual_repair(
          make_destroyed_fixture(granular_initial_solution,
                                 virtual_oracle_factory, removed),
          operation, instance, virtual_oracle_factory,
          virtual_oracle_paths, 100.0);
      require(certified.customer_sequences() == oracle.customer_sequences() &&
                  certified_factory.gap_map(certified) ==
                      virtual_oracle_factory.gap_map(oracle),
              "lower-bound certified virtual repair differs from exhaustive "
              "virtual completion");
    }
    std::mt19937_64 forced_granular_initial_random(36'509);
    ils_sp::InitialSolutionBuilder forced_granular_initial(
        instance, granular_factory, granular_paths,
        forced_granular_initial_random,
        ils_sp::SearchProfile::GranularAdaptive, 1, 0);
    (void)forced_granular_initial.build(100.0);
    require(forced_granular_initial.statistics()
                    .candidates_rejected_by_granular > 0 &&
                forced_granular_initial.statistics()
                        .granular_escape_candidates_selected > 0,
            "granular construction did not use the shared promising-arc "
            "contract and its reservoir escape above the threshold");
    ils_sp::Perturbation granular_perturbation(
        instance, granular_factory, granular_paths, granular_random,
        ils_sp::SearchProfile::GranularAdaptive, 1, 0);
    const std::vector<int> granular_removed{
        instance.customer_ids()[0], instance.customer_ids()[1],
        instance.customer_ids()[2]};
    const auto preserves_repair_paths = [&](const ils_sp::Solution& solution,
                                            const auto& inherited) {
      for (const ils_sp::Plan& plan : solution.plans) {
        std::vector<int> anchors;
        anchors.reserve(plan.customer_ids.size() + 2);
        anchors.push_back(instance.depot().id);
        anchors.insert(anchors.end(), plan.customer_ids.begin(),
                       plan.customer_ids.end());
        anchors.push_back(instance.depot().id);
        for (std::size_t gap = 0; gap < plan.gaps.size(); ++gap) {
          const auto known = inherited.find(
              ils_sp::AnchorPair{anchors[gap], anchors[gap + 1]});
          if (known == inherited.end()) {
            if (!plan.gaps[gap].empty()) return false;
          } else if (plan.gaps[gap] != known->second) {
            return false;
          }
        }
      }
      return true;
    };
    const ils_sp::DestroyedSolution granular_destroyed =
        make_destroyed_fixture(granular_initial_solution, granular_factory,
                               granular_removed);
    const auto granular_inherited =
        granular_factory.gap_map(granular_destroyed.remaining);
    const ils_sp::Solution granular_repaired =
        ils_sp::PerturbationTestAccess::repair(
            granular_perturbation, granular_destroyed,
            ils_sp::RepairOperator::GreedyInsertion, 100.0);
    granular_repaired.validate_partition(instance);
    require(granular_perturbation.statistics()
                    .candidates_rejected_by_granular > 0 &&
                preserves_repair_paths(granular_repaired,
                                       granular_inherited),
            "granular repair changed charging paths or ignored its "
            "promising-arc contract");
    const auto& granular_repair_statistics =
        granular_perturbation.statistics();
    require(
        granular_repair_statistics.decisions == granular_removed.size() &&
            granular_repair_statistics.greedy_decisions ==
                granular_removed.size() &&
            granular_repair_statistics.candidates_generated ==
                granular_repair_statistics.candidates_evaluated &&
            granular_repair_statistics.virtual_candidates_evaluated +
                    granular_repair_statistics
                        .virtual_candidates_pruned_by_bound ==
                granular_repair_statistics.candidates_evaluated &&
            granular_repair_statistics.virtual_completion_memo_hits +
                    granular_repair_statistics
                        .virtual_completion_memo_misses ==
                granular_repair_statistics.virtual_candidates_evaluated &&
            granular_repair_statistics.candidates_ranked_by_rough_bound == 0 &&
            granular_repair_statistics.completion_memo_hits == 0 &&
            granular_repair_statistics.completion_memo_misses == 0 &&
            granular_repair_statistics.affected_path_completion.calls == 0,
        "certified lazy Greedy repair did not close its lower-bound ranking "
        "without materializing charging paths");

    std::mt19937_64 granular_regret_random(36'509);
    ils_sp::Perturbation granular_regret_perturbation(
        instance, granular_factory, granular_paths, granular_regret_random,
        ils_sp::SearchProfile::GranularAdaptive, 1, 0);
    const ils_sp::Solution granular_regret_repaired =
        ils_sp::PerturbationTestAccess::repair(
            granular_regret_perturbation, granular_destroyed,
            ils_sp::RepairOperator::TwoRegretInsertion, 100.0);
    granular_regret_repaired.validate_partition(instance);
    const auto& granular_regret_statistics =
        granular_regret_perturbation.statistics();
    require(
        granular_regret_statistics.decisions == granular_removed.size() &&
            granular_regret_statistics.two_regret_decisions ==
                granular_removed.size() &&
            granular_regret_statistics.candidates_generated ==
                granular_regret_statistics.candidates_evaluated &&
            granular_regret_statistics.virtual_candidates_evaluated +
                    granular_regret_statistics
                        .virtual_candidates_pruned_by_bound +
                    granular_regret_statistics
                        .virtual_candidates_pruned_by_regret ==
                granular_regret_statistics.candidates_evaluated &&
            granular_regret_statistics.virtual_completion_memo_hits +
                    granular_regret_statistics
                        .virtual_completion_memo_misses ==
                granular_regret_statistics.virtual_candidates_evaluated &&
            granular_regret_statistics.candidates_ranked_by_rough_bound == 0 &&
            granular_regret_statistics.completion_memo_hits == 0 &&
            granular_regret_statistics.completion_memo_misses == 0 &&
            granular_regret_statistics.affected_path_completion.calls == 0 &&
            preserves_repair_paths(granular_regret_repaired,
                                   granular_inherited),
        "certified lazy 2-Regret did not close every customer ranking while "
        "preserving direct repair output");

    std::mt19937_64 hinted_repair_random(44'001);
    ils_sp::Perturbation hinted_repair(
        instance, granular_factory, granular_paths, hinted_repair_random,
        ils_sp::SearchProfile::GranularAdaptive);
    const ils_sp::Solution hinted_solution =
        hinted_repair.apply(granular_initial_solution, 100.0);
    hinted_solution.validate_partition(instance);
    require(hinted_repair.statistics().post_repair_hint_attempts > 0 &&
                hinted_repair.statistics().post_repair_hint_accepted > 0 &&
                std::any_of(hinted_solution.plans.begin(),
                            hinted_solution.plans.end(),
                            [](const ils_sp::Plan& plan) {
                              return std::any_of(
                                  plan.gaps.begin(), plan.gaps.end(),
                                  [](const std::vector<int>& gap) {
                                    return !gap.empty();
                                  });
                            }),
            "post-repair ReplacePath did not independently accept a virtual "
            "insertion hint");

    ils_sp::VariableNeighborhoodDescent granular_vnd(
        instance, granular_factory, granular_paths, granular_random, 1,
        ils_sp::SearchProfile::GranularAdaptive);
    const double granular_repaired_cost = ils_sp::generalized_cost(
        granular_repaired, instance, 100.0);
    const ils_sp::Solution granular_vnd_solution =
        granular_vnd.improve(granular_repaired, 100.0);
    const auto& granular_vnd_statistics = granular_vnd.statistics();
    std::size_t granular_replaceable_gap_count = 0;
    for (const ils_sp::Plan& plan : granular_vnd_solution.plans) {
      std::vector<int> anchors;
      anchors.reserve(plan.customer_ids.size() + 2);
      anchors.push_back(instance.depot().id);
      anchors.insert(anchors.end(), plan.customer_ids.begin(),
                     plan.customer_ids.end());
      anchors.push_back(instance.depot().id);
      for (std::size_t gap = 0; gap < plan.gaps.size(); ++gap) {
        const auto options =
            granular_paths.paths_between(anchors[gap], anchors[gap + 1]);
        if (std::any_of(options.begin(), options.end(),
                        [&](const std::vector<int>& option) {
                          return option != plan.gaps[gap];
                        })) {
          ++granular_replaceable_gap_count;
        }
      }
    }
    require(ils_sp::generalized_cost(granular_vnd_solution, instance, 100.0) <=
                granular_repaired_cost + ils_sp::kCostTolerance &&
                granular_vnd_statistics.adaptive_customer_advances == 3 &&
                granular_vnd_statistics.adaptive_replace_path_selections >
                    0 &&
                granular_vnd_statistics.adaptive_replace_path_selections <=
                    3 &&
                granular_vnd_statistics.adaptive_final_replace_path_calls ==
                    1 &&
                granular_vnd_statistics
                        .adaptive_final_replace_path_gaps_scanned ==
                    granular_replaceable_gap_count &&
                granular_vnd_statistics.replace_path_calls >=
                    granular_vnd_statistics
                        .adaptive_final_replace_path_gaps_scanned &&
                granular_vnd_statistics.linked_completion_attempts == 0,
            "granular-adaptive VND violated its stochastic transition "
            "contract");

    ils_sp::RouteEvaluator linked_evaluator(instance, 256);
    ils_sp::PlanFactory linked_factory(instance, linked_evaluator, 64);
    auto one_gaps = ils_sp::direct_gaps(1);
    one_gaps[1] = {12};
    auto seven_four_nine_gaps = ils_sp::direct_gaps(3);
    seven_four_nine_gaps[1] = {11};
    auto ten_six_gaps = ils_sp::direct_gaps(2);
    ten_six_gaps[0] = {12};
    ten_six_gaps[2] = {12};
    ils_sp::Solution linked_basin{{
        linked_factory.make_plan({1}, std::move(one_gaps)),
        linked_factory.make_plan({2, 5}, ils_sp::direct_gaps(2)),
        linked_factory.make_plan({8, 3}, ils_sp::direct_gaps(2)),
        linked_factory.make_plan({7, 4, 9}, std::move(seven_four_nine_gaps)),
        linked_factory.make_plan({10, 6}, std::move(ten_six_gaps)),
    }};
    linked_basin.validate_partition(instance);
    require(linked_basin.feasible(),
            "linked-path regression basin must be feasible");
    std::mt19937_64 linked_random(36'509);
    ils_sp::PathSampler linked_paths(instance, linked_factory, linked_random);
    ils_sp::VariableNeighborhoodDescent linked_vnd(
        instance, linked_factory, linked_paths, linked_random);
    const double linked_basin_cost =
        ils_sp::generalized_cost(linked_basin, instance, 100.0);
    const ils_sp::Solution linked =
        linked_vnd.improve(linked_basin, 100.0);
    require(ils_sp::generalized_cost(linked, instance, 100.0) <=
                linked_basin_cost + ils_sp::kCostTolerance,
            "VND accepted a non-improving materialized customer move");
    ils_sp::VariableNeighborhoodDescent linked_profile_vnd(
        instance, factory, insertion_paths, random, 40,
        ils_sp::SearchProfile::LinkedWinner);
    const double linked_perturbed_cost =
        ils_sp::generalized_cost(linked_perturbed, instance, 100.0);
    const ils_sp::Solution linked_profile_solution =
        linked_profile_vnd.improve(linked_perturbed, 100.0);
    require(ils_sp::generalized_cost(linked_profile_solution, instance,
                                     100.0) <=
                linked_perturbed_cost + ils_sp::kCostTolerance &&
                linked_profile_vnd.statistics().linked_completion_attempts >
                    0,
            "linked-winner VND did not use its isolated completion profile");
    const int depot_id = instance.depot().id;
    const auto has_cs_option = [&](int origin, int target) {
      const auto options = insertion_paths.paths_between(origin, target);
      return std::any_of(options.begin(), options.end(),
                         [](const std::vector<int>& path) {
                           return !path.empty();
                         });
    };
    const auto bootstrap_customer = std::find_if(
        instance.customer_ids().begin(), instance.customer_ids().end(),
        [&](int customer_id) {
          return has_cs_option(depot_id, customer_id) &&
                 has_cs_option(customer_id, depot_id);
        });
    require(bootstrap_customer != instance.customer_ids().end(),
            "ReplacePath bootstrap fixture has no two-sided CS alternative");
    const ils_sp::Solution direct_only{{factory.make_plan(
        {*bootstrap_customer}, ils_sp::direct_gaps(1))}};
    const auto limited_customer = std::find_if(
        instance.customer_ids().begin(), instance.customer_ids().end(),
        [&](int customer_id) {
          const ils_sp::Plan direct =
              factory.make_plan({customer_id}, ils_sp::direct_gaps(1));
          return direct.evaluation->feasible &&
                 has_cs_option(depot_id, customer_id);
        });
    require(limited_customer != instance.customer_ids().end(),
            "limited ReplacePath fixture has no feasible direct gap");
    std::mt19937_64 limited_random(72'901);
    ils_sp::PathSampler limited_paths(instance, factory, limited_random);
    const ils_sp::Solution limited_direct{{factory.make_plan(
        {*limited_customer}, ils_sp::direct_gaps(1))}};
    const std::array<ils_sp::AnchorPair, 1> limited_dirty{
        ils_sp::AnchorPair{depot_id, *limited_customer}};
    const ils_sp::LimitedRandomPathResult limited_result =
        limited_paths.limited_random_replacement(limited_direct, 0.1,
                                                 limited_dirty);
    require(limited_result.eligible_gaps == 1 &&
                limited_result.failure_limit == 1 &&
                limited_result.attempts == 1 &&
                limited_result.failed_attempts == 1 &&
                limited_result.accepted == 0 &&
                limited_result.solution.plans[0].gaps ==
                    limited_direct.plans[0].gaps,
            "limited ReplacePath repeated a small dirty gap or accepted a "
            "non-improving atomic move");
    const auto direct_replaced =
        insertion_paths.random_replacement(direct_only);
    require(direct_replaced.has_value() &&
                std::count_if(
                    direct_replaced->plans[0].gaps.begin(),
                    direct_replaced->plans[0].gaps.end(),
                    [](const std::vector<int>& gap) { return !gap.empty(); }) ==
                    1,
            "ReplacePath could not introduce a CS path into a direct route");

    const auto repair_customer = std::find_if(
        instance.customer_ids().begin(), instance.customer_ids().end(),
        [&](int customer_id) {
          const ils_sp::Plan direct = factory.make_plan(
              {customer_id}, ils_sp::direct_gaps(1));
          return !direct.evaluation->feasible &&
                 has_cs_option(depot_id, customer_id) &&
                 has_cs_option(customer_id, depot_id);
        });
    require(repair_customer != instance.customer_ids().end(),
            "path-repair fixture has no energy-infeasible two-sided route");
    const ils_sp::Solution repair_direct{{factory.make_plan(
        {*repair_customer}, ils_sp::direct_gaps(1))}};
    const ils_sp::Plan expired_completion = insertion_paths.repair_new_gaps(
        repair_direct.plans[0], ils_sp::PlanFactory::GapMap{},
        10'000.0, std::chrono::steady_clock::now() -
                      std::chrono::milliseconds(1));
    require(expired_completion.gaps == repair_direct.plans[0].gaps,
            "expired new-gap completion changed its incumbent");
    const ils_sp::Plan repaired_singleton = insertion_paths.repair_new_gaps(
        repair_direct.plans[0], ils_sp::PlanFactory::GapMap{},
        10'000.0);
    require(ils_sp::generalized_cost(*repaired_singleton.evaluation, instance,
                                     10'000.0) <
                ils_sp::generalized_cost(
                    *repair_direct.plans[0].evaluation, instance,
                    10'000.0) -
                    ils_sp::kCostTolerance,
            "introduced-gap completion did not strictly improve its plan");
    ils_sp::PlanFactory::GapMap inherited_outbound;
    inherited_outbound.emplace(
        ils_sp::AnchorPair{depot_id, *repair_customer},
        repair_direct.plans[0].gaps[0]);
    const ils_sp::Plan directionally_repaired =
        insertion_paths.repair_new_gaps(repair_direct.plans[0],
                                        inherited_outbound, 10'000.0);
    require(directionally_repaired.gaps[0] ==
                repair_direct.plans[0].gaps[0],
            "new-gap completion modified an inherited directed gap");
    auto replaceable_gaps = ils_sp::direct_gaps(1);
    const auto outbound_options =
        insertion_paths.paths_between(depot_id, *bootstrap_customer);
    replaceable_gaps[0] = *std::find_if(
        outbound_options.begin(), outbound_options.end(),
        [](const std::vector<int>& path) { return !path.empty(); });
    const ils_sp::Solution one_cs_path{{
        factory.make_plan({*bootstrap_customer}, std::move(replaceable_gaps)),
    }};
    const auto replaced = insertion_paths.random_replacement(one_cs_path);
    require(replaced.has_value() &&
                replaced->plans[0].customer_ids ==
                    one_cs_path.plans[0].customer_ids &&
                std::inner_product(
                    replaced->plans[0].gaps.begin(),
                    replaced->plans[0].gaps.end(),
                    one_cs_path.plans[0].gaps.begin(), std::size_t{0},
                    std::plus<>(), std::not_equal_to<>()) == 1,
            "ReplacePath did not change exactly one sampled path");

    // Rough insertion completion must also be able to remove an inherited
    // charging path when direct travel is the better endpoint.  This is the
    // part that repair_new_gaps intentionally does not cover.
    std::optional<ils_sp::Plan> removable_cs_plan;
    for (const int customer_id : instance.customer_ids()) {
      const ils_sp::Plan direct =
          factory.make_plan({customer_id}, ils_sp::direct_gaps(1));
      for (std::size_t gap = 0; gap < direct.gaps.size(); ++gap) {
        const int origin = gap == 0 ? depot_id : customer_id;
        const int target = gap == 0 ? customer_id : depot_id;
        for (const auto& option : insertion_paths.paths_between(origin,
                                                                 target)) {
          if (option.empty()) continue;
          auto candidate_gaps = direct.gaps;
          candidate_gaps[gap] = option;
          ils_sp::Plan candidate = factory.make_plan(
              {customer_id}, std::move(candidate_gaps));
          if (direct.evaluation->feasible &&
              ils_sp::generalized_cost(*direct.evaluation, instance, 100.0) <
              ils_sp::generalized_cost(*candidate.evaluation, instance,
                                       100.0) -
                  ils_sp::kCostTolerance) {
            removable_cs_plan = std::move(candidate);
            break;
          }
        }
        if (removable_cs_plan.has_value()) break;
      }
      if (removable_cs_plan.has_value()) break;
    }
    require(removable_cs_plan.has_value(),
            "rough repair fixture has no removable charging path");
    const ils_sp::Solution removable_solution{{*removable_cs_plan}};
    const auto removable_inherited = factory.gap_map(removable_solution);
    const double removable_incumbent = ils_sp::generalized_cost(
        *removable_cs_plan->evaluation, instance, 100.0);
    const ils_sp::Plan removed_cs = insertion_paths.repair_affected_paths(
        *removable_cs_plan, removable_inherited, 100.0);
    const double removed_cost = ils_sp::generalized_cost(
        *removed_cs.evaluation, instance, 100.0);
    const auto removed_count = std::count_if(
        removed_cs.gaps.begin(), removed_cs.gaps.end(),
        [](const std::vector<int>& gap) { return !gap.empty(); });
    const auto incumbent_count = std::count_if(
        removable_cs_plan->gaps.begin(), removable_cs_plan->gaps.end(),
        [](const std::vector<int>& gap) { return !gap.empty(); });
    require(removed_cost < removable_incumbent - ils_sp::kCostTolerance &&
                removed_count < incumbent_count,
            "rough repair could not delete an inherited charging path: "
            "incumbent=" +
                std::to_string(removable_incumbent) +
                " repaired=" + std::to_string(removed_cost) +
                " old_cs=" + std::to_string(incumbent_count) +
                " new_cs=" + std::to_string(removed_count));

    const std::array one_customer_anchors{
        ils_sp::AnchorPair{depot_id, *bootstrap_customer},
        ils_sp::AnchorPair{*bootstrap_customer, depot_id}};
    const double quasi_incumbent =
        ils_sp::generalized_cost(one_cs_path, instance, 100.0);
    const ils_sp::QuasiExhaustivePathResult quasi =
        insertion_paths.quasi_exhaustive_replacement(
            one_cs_path, 100.0, one_customer_anchors);
    require(quasi.gaps_scanned == one_customer_anchors.size() &&
                ils_sp::generalized_cost(quasi.solution, instance, 100.0) <=
                    quasi_incumbent + ils_sp::kCostTolerance,
            "final ReplacePath did not scan each eligible gap exactly once");

    ils_sp::LruCache<int, int> lru(2);
    lru.put(1, 10);
    lru.put(2, 20);
    require(lru.get(1).value() == 10, "LRU lookup failed");
    lru.put(3, 30);
    require(!lru.get(2).has_value() && lru.get(3).value() == 30,
            "LRU did not evict the least recently used entry");

    std::cout << "all C++ tests passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "test failure: " << error.what() << '\n';
    return 1;
  }
}
