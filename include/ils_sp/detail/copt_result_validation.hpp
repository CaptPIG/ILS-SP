#pragma once

#include <cmath>

#include "copt.h"
#include "ils_sp/evaluator.hpp"

namespace ils_sp::detail {

[[nodiscard]] inline bool finite_copt_value(double value) {
  return std::isfinite(value) && std::abs(value) < COPT_INFINITY;
}

// Keep physical feasibility stricter than the backend's feasibility tolerance.
[[nodiscard]] inline bool valid_frvcp_replay(
    const RouteEvaluation& evaluation, double model_objective) {
  return evaluation.feasible && finite_copt_value(evaluation.raw_cost_h) &&
         finite_copt_value(model_objective) &&
         std::abs(evaluation.raw_cost_h - model_objective) <= 1e-5;
}

}  // namespace ils_sp::detail
