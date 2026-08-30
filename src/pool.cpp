#include "ils_sp/pool.hpp"

#include <cmath>
#include <stdexcept>
#include <tuple>

namespace ils_sp {
namespace {

[[nodiscard]] bool column_rank_less(
    const RouteColumn& left, const RouteColumn& right,
    const std::unordered_map<int, double>& duals) {
  const double left_reduced_cost = left.reduced_cost(duals);
  const double right_reduced_cost = right.reduced_cost(duals);
  if (left_reduced_cost != right_reduced_cost) {
    return left_reduced_cost < right_reduced_cost;
  }
  if (left.cost_h != right.cost_h) return left.cost_h < right.cost_h;
  if (left.customer_set != right.customer_set) {
    return left.customer_set < right.customer_set;
  }
  return left.signature < right.signature;
}

[[nodiscard]] bool incumbent_dominates(
    const RouteColumn& incumbent, const RouteColumn& candidate,
    double tolerance) {
  if (incumbent.cost_h < candidate.cost_h - tolerance) return true;
  if (std::abs(incumbent.cost_h - candidate.cost_h) > tolerance) return false;
  return incumbent.plan.exact_charging || !candidate.plan.exact_charging;
}

}  // namespace

double RouteColumn::reduced_cost(
    const std::unordered_map<int, double>& duals) const {
  double result = cost_h;
  for (const int customer_id : customer_set) {
    const auto found = duals.find(customer_id);
    if (found != duals.end()) result -= found->second;
  }
  return result;
}

RouteColumn RouteColumn::from_plan(const Plan& plan) {
  if (plan.evaluation == nullptr || !plan.evaluation->feasible ||
      plan.customer_ids.empty()) {
    throw std::invalid_argument("SP columns must be nonempty feasible routes");
  }
  const std::vector<int> customer_set = sorted_customer_set(plan.customer_ids);
  if (customer_set.size() != plan.customer_ids.size()) {
    throw std::invalid_argument("SP columns cannot repeat customers");
  }
  return RouteColumn{.plan = plan,
                     .customer_set = customer_set,
                     .cost_h = plan.evaluation->raw_cost_h,
                     .signature = plan.route.visits};
}

RoutePool::RoutePool(std::size_t capacity, double tolerance)
    : capacity_(capacity), tolerance_(tolerance) {
  if (capacity_ == 0 || tolerance_ < 0.0) {
    throw std::invalid_argument("invalid route-pool configuration");
  }
  main_.reserve(capacity_);
}

void RoutePool::put_main(RouteColumn column) {
  const std::vector<int> customer_key = column.plan.customer_ids;
  backup_.erase(customer_key);
  const auto incumbent = main_.find(customer_key);
  if (incumbent != main_.end()) index_main_erase(incumbent->second);
  const auto stored =
      main_.insert_or_assign(customer_key, std::move(column)).first;
  index_main_insert(stored->second);
  ++statistics_.main_insertions;
}

std::optional<double> RoutePool::best_coverage_cost(
    const std::vector<int>& coverage) const {
  const auto found = coverage_costs_.find(coverage);
  if (found == coverage_costs_.end() || found->second.empty()) {
    return std::nullopt;
  }
  return *found->second.begin();
}

bool RoutePool::coverage_improved(std::optional<double> before,
                                  std::optional<double> after) const {
  return after.has_value() &&
         (!before.has_value() || *after < *before - tolerance_);
}

void RoutePool::index_main_insert(const RouteColumn& column) {
  coverage_costs_[column.customer_set].insert(column.cost_h);
}

void RoutePool::index_main_erase(const RouteColumn& column) {
  const auto coverage = coverage_costs_.find(column.customer_set);
  if (coverage == coverage_costs_.end()) {
    throw std::logic_error("main-pool coverage index is missing");
  }
  const auto cost = coverage->second.find(column.cost_h);
  if (cost == coverage->second.end()) {
    throw std::logic_error("main-pool coverage cost is missing");
  }
  coverage->second.erase(cost);
  if (coverage->second.empty()) coverage_costs_.erase(coverage);
}

void RoutePool::erase_main(ColumnMap::iterator column) {
  index_main_erase(column->second);
  main_.erase(column);
}

void RoutePool::demote_main(ColumnMap::iterator column) {
  RouteColumn demoted = std::move(column->second);
  index_main_erase(demoted);
  main_.erase(column);
  put_backup(std::move(demoted));
  ++statistics_.main_demotions;
}

void RoutePool::put_backup(RouteColumn column) {
  const auto found = backup_.find(column.plan.customer_ids);
  if (found == backup_.end() ||
      column.cost_h < found->second.cost_h - tolerance_ ||
      (std::abs(column.cost_h - found->second.cost_h) <= tolerance_ &&
       column.plan.exact_charging && !found->second.plan.exact_charging)) {
    const std::vector<int> customer_key = column.plan.customer_ids;
    backup_.insert_or_assign(customer_key, std::move(column));
    ++statistics_.backup_insertions;
  }
}

void RoutePool::record_effective_coverage_improvement(bool improved) {
  if (!improved) return;
  ++effective_generation_;
  ++statistics_.effective_coverage_improvements;
}

RoutePool::ColumnMap::iterator RoutePool::worst_main() {
  if (main_.empty()) {
    throw std::logic_error("cannot evict from an empty main pool");
  }
  auto worst = std::find_if(main_.begin(), main_.end(), [&](const auto& item) {
    return !protected_sequences_.contains(item.first);
  });
  if (worst == main_.end()) {
    throw std::logic_error("route pool has no evictable main column");
  }
  for (auto current = std::next(worst); current != main_.end(); ++current) {
    if (protected_sequences_.contains(current->first)) continue;
    if (column_rank_less(worst->second, current->second, duals_)) {
      worst = current;
    }
  }
  return worst;
}

PoolUpdate RoutePool::admit(RouteColumn column) {
  ++statistics_.admission_attempts;
  const auto incumbent = main_.find(column.plan.customer_ids);
  if (incumbent != main_.end()) {
    if (incumbent_dominates(incumbent->second, column, tolerance_)) {
      ++statistics_.dominated_rejections;
      return PoolUpdate{};
    }
    const std::optional<double> old_coverage_cost =
        best_coverage_cost(column.customer_set);
    const std::vector<int> replaced = incumbent->first;
    put_main(std::move(column));
    ++statistics_.same_sequence_replacements;
    const bool effective = coverage_improved(
        old_coverage_cost, best_coverage_cost(main_.at(replaced).customer_set));
    record_effective_coverage_improvement(effective);
    return PoolUpdate{.main = true,
                      .replaced_dominated = true,
                      .effective_coverage_improvement = effective,
                      .evicted_customer_sequence = replaced};
  }
  const auto backup_incumbent = backup_.find(column.plan.customer_ids);
  if (backup_incumbent != backup_.end() &&
      incumbent_dominates(backup_incumbent->second, column, tolerance_)) {
    ++statistics_.dominated_rejections;
    return PoolUpdate{};
  }
  if (main_.size() < capacity_) {
    const std::vector<int> coverage = column.customer_set;
    const std::optional<double> old_coverage_cost =
        best_coverage_cost(coverage);
    put_main(std::move(column));
    const bool effective =
        coverage_improved(old_coverage_cost, best_coverage_cost(coverage));
    record_effective_coverage_improvement(effective);
    return PoolUpdate{.main = true,
                      .backup = false,
                      .replaced_dominated = false,
                      .effective_coverage_improvement = effective,
                      .evicted_customer_sequence = std::nullopt};
  }
  if (duals_valid_ && column.reduced_cost(duals_) >= -tolerance_) {
    put_backup(std::move(column));
    return PoolUpdate{.main = false,
                      .backup = true,
                      .replaced_dominated = false,
                      .evicted_customer_sequence = std::nullopt};
  }

  const std::vector<int> new_sequence = column.plan.customer_ids;
  const std::vector<int> new_coverage = column.customer_set;
  const std::optional<double> old_coverage_cost =
      best_coverage_cost(new_coverage);
  put_main(std::move(column));
  auto worst = worst_main();
  const std::vector<int> evicted = worst->first;
  const bool new_survived = evicted != new_sequence;
  demote_main(worst);
  const bool effective =
      new_survived &&
      coverage_improved(old_coverage_cost, best_coverage_cost(new_coverage));
  record_effective_coverage_improvement(effective);
  return PoolUpdate{.main = new_survived,
                    .backup = !new_survived,
                    .replaced_dominated = false,
                    .effective_coverage_improvement = effective,
                    .evicted_customer_sequence = evicted};
}

PoolRepriceUpdate RoutePool::reprice_backup() {
  if (!duals_valid_) return PoolRepriceUpdate{};
  std::vector<RouteColumn> candidates;
  for (const auto& [_, column] : backup_) {
    if (column.reduced_cost(duals_) < -tolerance_) {
      candidates.push_back(column);
    }
  }
  std::sort(candidates.begin(), candidates.end(), [&](const RouteColumn& left,
                                                       const RouteColumn& right) {
    return column_rank_less(left, right, duals_);
  });
  PoolRepriceUpdate update{.candidates = candidates.size()};
  statistics_.repricing_candidates += candidates.size();
  for (RouteColumn& candidate : candidates) {
    backup_.erase(candidate.plan.customer_ids);
    const PoolUpdate admitted = admit(std::move(candidate));
    if (admitted.main) ++update.promoted;
    if (admitted.effective_coverage_improvement) {
      ++update.effective_coverage_improvements;
    }
  }
  statistics_.repricing_promotions += update.promoted;
  return update;
}

void RoutePool::update_duals(std::unordered_map<int, double> duals) {
  duals_ = std::move(duals);
  duals_valid_ = true;
}

void RoutePool::protect_partition(const Solution& solution) {
  if (!solution.feasible() || solution.plans.size() > capacity_) {
    throw std::invalid_argument(
        "protected route partition must be feasible and fit the main pool");
  }
  std::unordered_set<std::vector<int>, IntVectorHash> protected_sequences;
  protected_sequences.reserve(solution.plans.size());
  for (const Plan& plan : solution.plans) {
    if (!protected_sequences.insert(plan.customer_ids).second) {
      throw std::invalid_argument(
          "protected route partition repeats an ordered customer sequence");
    }
  }
  protected_sequences_ = std::move(protected_sequences);

  for (const Plan& plan : solution.plans) {
    RouteColumn candidate = RouteColumn::from_plan(plan);
    const std::vector<int> coverage = candidate.customer_set;
    const std::optional<double> old_coverage_cost =
        best_coverage_cost(coverage);
    const auto incumbent = main_.find(plan.customer_ids);
    if (incumbent == main_.end() ||
        !incumbent_dominates(incumbent->second, candidate, tolerance_)) {
      put_main(std::move(candidate));
      record_effective_coverage_improvement(
          coverage_improved(old_coverage_cost, best_coverage_cost(coverage)));
    }
  }
  while (main_.size() > capacity_) {
    demote_main(worst_main());
  }
}

std::vector<RouteColumn> RoutePool::main_columns() const {
  std::vector<RouteColumn> result;
  result.reserve(main_.size());
  for (const auto& [_, column] : main_) result.push_back(column);
  std::sort(result.begin(), result.end(), [](const RouteColumn& left,
                                             const RouteColumn& right) {
    return std::tie(left.customer_set, left.cost_h, left.signature) <
           std::tie(right.customer_set, right.cost_h, right.signature);
  });
  return result;
}

std::vector<RouteColumn> RoutePool::backup_columns() const {
  std::vector<RouteColumn> result;
  result.reserve(backup_.size());
  for (const auto& [_, column] : backup_) result.push_back(column);
  std::sort(result.begin(), result.end(), [&](const RouteColumn& left,
                                              const RouteColumn& right) {
    return column_rank_less(left, right, duals_);
  });
  return result;
}

std::size_t RoutePool::distinct_main_coverages() const {
  return coverage_costs_.size();
}

}  // namespace ils_sp
