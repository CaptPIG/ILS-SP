#pragma once

#include <set>
#include <unordered_map>
#include <unordered_set>

#include "ils_sp/core.hpp"

namespace ils_sp {

struct RouteColumn {
  Plan plan;
  std::vector<int> customer_set;
  double cost_h{};
  std::vector<int> signature;

  [[nodiscard]] double reduced_cost(
      const std::unordered_map<int, double>& duals) const;
  [[nodiscard]] static RouteColumn from_plan(const Plan& plan);
};

struct PoolUpdate {
  bool main{};
  bool backup{};
  bool replaced_dominated{};
  // True only when the admitted column creates a new customer coverage or
  // strictly lowers the cheapest main-pool cost for that coverage.  Different
  // ordered sequences with the same, no-better coverage do not count.
  bool effective_coverage_improvement{};
  std::optional<std::vector<int>> evicted_customer_sequence;
};

struct PoolRepriceUpdate {
  std::size_t candidates{};
  std::size_t promoted{};
  std::size_t effective_coverage_improvements{};
};

struct PoolStatistics {
  std::uint64_t admission_attempts{};
  std::uint64_t dominated_rejections{};
  std::uint64_t same_sequence_replacements{};
  std::uint64_t main_insertions{};
  std::uint64_t backup_insertions{};
  std::uint64_t main_demotions{};
  std::uint64_t repricing_candidates{};
  std::uint64_t repricing_promotions{};
  std::uint64_t effective_coverage_improvements{};
};

class RoutePool {
 public:
  explicit RoutePool(std::size_t capacity = 5'000,
                     double tolerance = kCostTolerance);

  [[nodiscard]] PoolUpdate admit(RouteColumn column);
  [[nodiscard]] PoolRepriceUpdate reprice_backup();
  void update_duals(std::unordered_map<int, double> duals);
  // Keep one representable column for every route of the incumbent best
  // feasible partition.  A cheaper/equally exact route with the same ordered
  // customer sequence may still replace the protected column.
  void protect_partition(const Solution& solution);

  [[nodiscard]] std::vector<RouteColumn> main_columns() const;
  [[nodiscard]] std::vector<RouteColumn> backup_columns() const;
  [[nodiscard]] std::size_t distinct_main_coverages() const;
  [[nodiscard]] bool duals_valid() const noexcept { return duals_valid_; }
  [[nodiscard]] std::uint64_t effective_generation() const noexcept {
    return effective_generation_;
  }
  [[nodiscard]] const PoolStatistics& statistics() const noexcept {
    return statistics_;
  }
  [[nodiscard]] const std::unordered_map<int, double>& duals() const noexcept {
    return duals_;
  }
  [[nodiscard]] std::size_t main_size() const noexcept { return main_.size(); }
  [[nodiscard]] std::size_t backup_size() const noexcept {
    return backup_.size();
  }
  [[nodiscard]] std::size_t capacity() const noexcept { return capacity_; }

 private:
  using ColumnMap =
      std::unordered_map<std::vector<int>, RouteColumn, IntVectorHash>;
  using CoverageCostMap = std::unordered_map<std::vector<int>,
                                             std::multiset<double>,
                                             IntVectorHash>;

  [[nodiscard]] ColumnMap::iterator worst_main();
  [[nodiscard]] std::optional<double> best_coverage_cost(
      const std::vector<int>& coverage) const;
  [[nodiscard]] bool coverage_improved(
      std::optional<double> before, std::optional<double> after) const;
  void index_main_insert(const RouteColumn& column);
  void index_main_erase(const RouteColumn& column);
  void erase_main(ColumnMap::iterator column);
  void demote_main(ColumnMap::iterator column);
  void put_main(RouteColumn column);
  void put_backup(RouteColumn column);
  void record_effective_coverage_improvement(bool improved);

  std::size_t capacity_;
  double tolerance_;
  ColumnMap main_;
  CoverageCostMap coverage_costs_;
  ColumnMap backup_;
  std::unordered_set<std::vector<int>, IntVectorHash> protected_sequences_;
  std::unordered_map<int, double> duals_;
  bool duals_valid_{false};
  std::uint64_t effective_generation_{};
  PoolStatistics statistics_;
};

}  // namespace ils_sp
