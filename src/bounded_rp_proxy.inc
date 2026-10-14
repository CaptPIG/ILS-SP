PathSampler::BoundedProxyScan PathSampler::bounded_proxy_choice(
    const Solution& solution,
    std::span<const std::pair<std::size_t, std::size_t>> gaps,
    std::size_t first, double penalty_lambda,
    std::size_t descriptor_cap, std::size_t work_cap) {
  BoundedProxyScan result;
  detail::RpProxyWork work{work_cap};
  const auto fallback = [&]() {
    result.choice.reset();
    result.work = work.used;
    result.completed = false;
    return result;
  };
  if (first >= gaps.size() || gaps.size() - first > descriptor_cap ||
      !std::isfinite(penalty_lambda) || penalty_lambda < 0) return fallback();
  // Temporary shared handles only. No persistent cache, path generation,
  // flattened catalogue, or first-N lexical approximation is introduced.
  std::vector<std::shared_ptr<const PathCatalog>> catalogs;
  catalogs.reserve(gaps.size() - first);
  for (std::size_t position = first; position < gaps.size(); ++position) {
    if (!work.take()) return fallback();
    const auto [route, gap] = gaps[position];
    if (route >= solution.plans.size()) return fallback();
    const Plan& plan = solution.plans[route];
    if (gap >= plan.gaps.size()) return fallback();
    const AnchorPair key{gap == 0 ? instance_.depot().id : plan.customer_ids[gap - 1],
        gap == plan.customer_ids.size() ? instance_.depot().id : plan.customer_ids[gap]};
    std::shared_ptr<const PathCatalog> catalog;
    if (catalog_cache_.contains(key)) {
      catalog = *catalog_cache_.get(key);
    } else {
      const auto record = catalog_disk_cache_.records.find(key);
      if (!catalog_disk_cache_.enabled || record == catalog_disk_cache_.records.end() ||
          !record->second.complete || record->second.header[6] == 0) return fallback();
      // Only packed existing disk records can supply a small indexed handle.
      catalog = disk_cache_load(key);
    }
    if (!catalog || (catalog->disk_backing && !catalog->packed) ||
        catalog->size() > descriptor_cap - result.raw_descriptors) return fallback();
    result.raw_descriptors += catalog->size();
    catalogs.push_back(std::move(catalog));
  }
  // Preflight counted EVERY raw descriptor before any winner can be chosen.
  // Ephemeral candidates are bounded by the same raw-descriptor/work caps.
  std::vector<BoundedProxyChoice> candidates;
  std::vector<double> scores;
  candidates.reserve(result.raw_descriptors);
  scores.reserve(result.raw_descriptors);
  for (std::size_t position = first; position < gaps.size(); ++position) {
    const auto [route, gap] = gaps[position];
    const Plan& plan = solution.plans[route];
    const auto& catalog = *catalogs[position - first];
    CatalogReader reader(catalog);
    for (std::size_t index = 0; index < catalog.size(); ++index) {
      if (!work.take()) return fallback();
      const auto descriptor = reader.descriptor(index);
      // Charge decode and equality comparison before either loops over IDs.
      const std::size_t count = catalog.path_size(descriptor);
      if (count > (work.cap - work.used) / 2 || !work.take(count * 2)) return fallback();
      PathCatalog::DecodeBuffer decoded{};
      const auto path = catalog.stations(descriptor, decoded);
      if (std::equal(path.begin(), path.end(), plan.gaps[gap].begin(), plan.gaps[gap].end()))
        continue;  // Exclude every physical copy of the current path.
      const auto route_proxy = detail::rp_route_proxy(
          instance_, plan, gap, path, penalty_lambda, work);
      if (!route_proxy) return fallback();
      double whole_proxy = 0;
      for (std::size_t r = 0; r < solution.plans.size(); ++r) {
        if (!work.take()) return fallback();
        whole_proxy += r == route ? *route_proxy
            : generalized_cost(*solution.plans[r].evaluation, instance_, penalty_lambda);
      }
      if (!std::isfinite(whole_proxy)) return fallback();
      if (!work.take(path.size())) return fallback();
      candidates.push_back(BoundedProxyChoice{
          position, std::vector<int>(path.begin(), path.end()), whole_proxy});
      scores.push_back(whole_proxy);
    }
  }
  // Charge both final minimum and band passes before any tie RNG draw.
  if (!work.take(scores.size()) || !work.take(scores.size())) return fallback();
  const auto ties = detail::rp_proxy_minimum_band(scores);
  if (!ties.empty()) {
    const std::size_t chosen = ties.size() == 1 ? ties.front()
        : ties[std::uniform_int_distribution<std::size_t>(0, ties.size() - 1)(random_)];
    result.choice = std::move(candidates[chosen]);
  }
  result.work = work.used;
  result.completed = true;
  return result;
}
