#pragma once

// Output-only diagnostic. This consumes the retained winning evaluations;
// it never re-evaluates, optimizes, mutates, or calls RNG/solver code.
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string_view>
#include "ils_sp/core.hpp"

namespace ils_sp_validation {
inline void quoted(std::ostream& out, std::string_view value) {
  out << '"';
  for (unsigned char c : value) {
    switch (c) {
      case '"': out << "\\\""; break;
      case '\\': out << "\\\\"; break;
      case '\n': out << "\\n"; break;
      case '\r': out << "\\r"; break;
      case '\t': out << "\\t"; break;
      default:
        if (c < 0x20) throw std::runtime_error("unsupported certificate string control character");
        out << static_cast<char>(c);
    }
  }
  out << '"';
}
template <class Range> inline void array(std::ostream& out, const Range& values) {
  out << '[';
  bool first = true;
  for (const auto& value : values) {
    if (!first) out << ',';
    first = false;
    out << value;
  }
  out << ']';
}
inline void write_certificate(const ils_sp::Instance& instance,
                              const ils_sp::Solution& solution,
                              std::string_view origin = "search_record") {
  const char* target = std::getenv("ILS_SP_VALIDATION_JSON");
  if (target == nullptr || *target == '\0') return;
  for (const auto& plan : solution.plans) {
    if (!plan.evaluation) throw std::runtime_error("missing retained route evaluation in schedule certificate");
  }
  std::ofstream out(target, std::ios::binary | std::ios::trunc);
  if (!out) throw std::runtime_error("cannot open schedule certificate output");
  out << std::setprecision(std::numeric_limits<double>::max_digits10);
  out << "{\"schema\":\"ils-sp-retained-schedule-v1\",\"instance\":";
  quoted(out, instance.name());
  out << ",\"solution_origin\":"; quoted(out, origin);
  out << ",\"objective_h\":" << solution.raw_cost_h() << ",\"routes\":[";
  for (std::size_t r = 0; r < solution.plans.size(); ++r) {
    if (r) out << ',';
    const auto& plan = solution.plans[r];
    const auto& e = *plan.evaluation;
    out << "{\"route_index\":" << r << ",\"visits\":";
    array(out, plan.route.visits);
    out << ",\"evaluation_visits\":"; array(out, e.route.visits);
    out << ",\"customers\":"; array(out, plan.customer_ids);
    out << ",\"arrival_energy_wh\":"; array(out, e.arrival_energy_wh);
    out << ",\"departure_energy_wh\":"; array(out, e.departure_energy_wh);
    out << ",\"feasible\":" << (e.feasible ? "true" : "false")
        << ",\"exact_charging\":" << (plan.exact_charging ? "true" : "false")
        << ",\"raw_cost_h\":" << e.raw_cost_h
        << ",\"travel_time_h\":" << e.travel_time_h
        << ",\"charging_time_h\":" << e.charging_time_h
        << ",\"service_time_h\":" << e.service_time_h
        << ",\"route_duration_h\":" << e.route_duration_h
        << ",\"energy_shortfall_wh\":" << e.energy_shortfall_wh
        << ",\"duration_excess_h\":" << e.duration_excess_h
        << ",\"charging_events\":[";
    for (std::size_t k = 0; k < e.charging_events.size(); ++k) {
      if (k) out << ',';
      const auto& event = e.charging_events[k];
      out << "{\"visit_index\":" << event.visit_index
          << ",\"station_id\":" << event.station_id
          << ",\"arrival_energy_wh\":" << event.arrival_energy_wh
          << ",\"departure_energy_wh\":" << event.departure_energy_wh
          << ",\"duration_h\":" << event.duration_h << '}';
    }
    out << "]}";
  }
  out << "]}\n";
  out.flush();
  if (!out) throw std::runtime_error("failed writing schedule certificate output");
}
}  // namespace ils_sp_validation
