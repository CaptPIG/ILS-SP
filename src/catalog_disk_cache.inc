// Included inside namespace ils_sp, after PathCatalog's representation helpers.
// This optional local cache does not define or approximate any path semantics.
namespace {
constexpr std::uint64_t kCatalogDiskMagic = 0x494c535043415431ULL;
constexpr std::uint64_t kCatalogDiskVersion = 2;
constexpr std::size_t kCatalogDiskChunk = 64 * 1024;
constexpr std::size_t kCatalogDiskMaximumRecords = 65536;
constexpr std::uint64_t kCatalogDiskReserve = 64ULL * 1024 * 1024;
constexpr std::uint64_t kCatalogHashOffset = 14695981039346656037ULL;

std::uint64_t catalog_budget(const char* name, std::uint64_t fallback) {
  const char* value = std::getenv(name);
  if (value == nullptr) return fallback;
  std::uint64_t result{};
  const char* end = value + std::strlen(value);
  const auto parsed = std::from_chars(value, end, result);
  // Explicit budgets, never host-memory autodetection. One TiB also bounds
  // arithmetic on this experimental format; larger requests fail clearly.
  if (value == end || parsed.ec != std::errc{} || parsed.ptr != end ||
      result > (1ULL << 40)) {
    throw std::invalid_argument(std::string("invalid catalogue byte budget: ") + name);
  }
  return result;
}

bool catalog_pread(int fd, void* destination, std::size_t size,
                   std::uint64_t offset, std::uint64_t& bytes) noexcept {
  auto* cursor = static_cast<unsigned char*>(destination);
  while (size != 0) {
    const auto n = ::pread(fd, cursor, size, static_cast<off_t>(offset));
    if (n < 0 && errno == EINTR) continue;
    if (n <= 0) return false;
    const auto amount = static_cast<std::size_t>(n);
    bytes += amount;
    cursor += amount;
    offset += amount;
    size -= amount;
  }
  return true;
}

std::uint64_t catalog_fp_signature() noexcept {
#if defined(__x86_64__)
  unsigned short control;
  asm volatile("fnstcw %0" : "=m"(control));
  // Control bits only: arithmetic status flags are intentionally irrelevant.
  return (static_cast<std::uint64_t>(control) & 0x1f3fULL) |
      ((static_cast<std::uint64_t>(_mm_getcsr()) & 0xffc0ULL) << 16);
#else
  return 0;
#endif
}

void catalog_hash_bytes(std::uint64_t& hash, const unsigned char* data,
                        std::size_t size) noexcept {
  for (std::size_t index = 0; index < size; ++index) {
    hash ^= data[index];
    hash *= 1099511628211ULL;
  }
}

struct CatalogDiskFile {
  std::FILE* file{};
  ~CatalogDiskFile() { if (file != nullptr) std::fclose(file); }
  bool close() noexcept {
    std::FILE* previous = file;
    file = nullptr;
    return std::fclose(previous) == 0;
  }
};

struct CatalogDiskTimer {
  double& seconds;
  std::chrono::steady_clock::time_point started{std::chrono::steady_clock::now()};
  ~CatalogDiskTimer() { seconds += elapsed_seconds(started); }
};

bool catalog_transfer(std::FILE* file, void* data, std::size_t size, bool writing,
                      std::uint64_t& byte_counter, std::uint64_t* hash,
                      std::vector<std::array<std::uint64_t, 2>>* blocks = nullptr) {
  auto* cursor = static_cast<unsigned char*>(data);
  while (size != 0) {
    const auto amount = std::min(size, kCatalogDiskChunk);
    const auto transferred = writing
        ? std::fwrite(cursor, 1, amount, file)
        : std::fread(cursor, 1, amount, file);
    byte_counter += transferred;
    if (hash != nullptr) catalog_hash_bytes(*hash, cursor, transferred);
    if (transferred != amount) return false;
    if (blocks != nullptr) {
      const auto block_hash = XXH3_128bits(cursor, amount);
      blocks->push_back({block_hash.low64, block_hash.high64});
    }
    cursor += amount;
    size -= amount;
  }
  return true;
}

std::filesystem::path catalog_disk_file(const std::filesystem::path& root,
                                       AnchorPair key) {
  return root / (std::to_string(key.origin) + "_" +
                 std::to_string(key.target) + ".cat");
}
}  // namespace

[[noreturn]] void PathSampler::CatalogReader::fail(const char* reason) const {
  ++backing->state->failures;
  throw std::runtime_error(std::string("indexed catalogue read failed: ") + reason +
                           " file=" + backing->file);
}

PathSampler::CatalogReader::CatalogReader(const PathCatalog& value)
    : catalog(value), backing(value.disk_backing) {
  if (!backing) return;  // No disk work or page bookkeeping on resident vectors.
  auto& state = *backing->state;
  if (state.active_readers >= 8) fail("active reader limit");
  CatalogDiskTimer timer{state.read_s};
  struct OpenFile {
    int value{-1};
    ~OpenFile() { if (value >= 0) ::close(value); }
  } file{::open(backing->file.c_str(), O_RDONLY | O_CLOEXEC | O_NOFOLLOW)};
  if (file.value < 0) fail("open");
  struct stat status{};
  if (::fstat(file.value, &status) != 0 || !S_ISREG(status.st_mode) ||
      status.st_size < 0 || static_cast<std::uint64_t>(status.st_size) != backing->header[2])
    fail("file length/type");
  CatalogDiskCache::Header header{};
  if (!catalog_pread(file.value, header.data(), sizeof(header), 0, state.bytes_read) ||
      header != backing->header) fail("header");
  fd = file.value;
  file.value = -1;
  ++state.active_readers;
  state.peak_readers = std::max(state.peak_readers, state.active_readers);
}

PathSampler::CatalogReader::~CatalogReader() {
  if (fd >= 0) {
    ::close(fd);
    --backing->state->active_readers;
  }
}

void PathSampler::CatalogReader::copy_value(
    void* destination, std::size_t width, std::uint64_t section,
    std::uint64_t section_bytes, std::uint64_t relative,
    std::size_t hash_base) const {
  if (relative > section_bytes || width > section_bytes - relative)
    fail("scalar outside section");
  const std::uint64_t block_number = relative / kCatalogDiskChunk;
  const std::uint64_t block_relative = block_number * kCatalogDiskChunk;
  const auto inside = static_cast<std::size_t>(relative - block_relative);
  const auto amount = static_cast<std::size_t>(
      std::min<std::uint64_t>(kCatalogDiskChunk, section_bytes - block_relative));
  if (inside > amount || width > amount - inside) fail("scalar crosses chunk");
  const auto hash_index = hash_base + static_cast<std::size_t>(block_number);
  if (hash_index >= backing->block_hashes.size()) fail("missing trusted block hash");
  auto& state = *backing->state;
  const CatalogBlockKey key{backing->record_id, section + block_relative};
  if (const auto cached = state.pages.get(key)) {
    ++state.hits;
    std::memcpy(destination, (*cached)->data() + inside, width);
    return;
  }
  ++state.misses;
  CatalogDiskTimer timer{state.read_s};
  // At most one transient 64-KiB page exists per active reader. Cache charges
  // include a conservative fixed allowance for the page and LRU/index nodes.
  auto page = std::make_shared<CatalogPageState::Page>(amount);
  if (!catalog_pread(fd, page->data(), amount, key.offset, state.bytes_read))
    fail("short block read");
  const auto hash = XXH3_128bits(page->data(), page->size());
  if (CatalogDiskBacking::BlockChecksum{hash.low64, hash.high64} !=
      backing->block_hashes[hash_index]) fail("block checksum");
  const std::size_t charge = page->capacity() + 256;
  if (charge <= state.budget) {
    while (state.bytes > state.budget - charge) {
      const auto old = state.pages.discard_oldest();
      if (!old) fail("page accounting");
      state.bytes -= (*old)->capacity() + 256;
    }
    try {
      state.pages.put(key, page);
    } catch (...) {
      // The legacy LRU appends its list before allocating its index node.
      // A failed insertion must not leave an unindexed, uncharged hot page.
      state.pages.clear();
      state.bytes = 0;
      throw;
    }
    state.bytes += charge;
    state.peak_bytes = std::max(state.peak_bytes, state.bytes);
  }
  std::memcpy(destination, page->data() + inside, width);
}

PathSampler::PathDescriptor PathSampler::CatalogReader::descriptor(std::size_t index) const {
  if (!backing) return catalog.paths.at(index);
  if (index >= catalog.size()) fail("descriptor index");
  PathDescriptor result;
  copy_value(&result, sizeof(result), sizeof(CatalogDiskCache::Header),
             backing->header[8] * sizeof(result), index * sizeof(result), 0);
  return result;
}

std::uint32_t PathSampler::CatalogReader::lexical_index(std::size_t ordinal) const {
  if (!backing) return catalog.lexical_order.at(ordinal);
  if (ordinal >= catalog.size()) fail("lexical ordinal");
  const std::uint64_t descriptor_bytes = backing->header[8] * sizeof(PathDescriptor);
  const std::size_t hash_base = static_cast<std::size_t>(
      (descriptor_bytes + kCatalogDiskChunk - 1) / kCatalogDiskChunk);
  std::uint32_t result;
  copy_value(&result, sizeof(result), sizeof(CatalogDiskCache::Header) + descriptor_bytes,
             backing->header[10] * sizeof(result), ordinal * sizeof(result), hash_base);
  if (result >= catalog.size()) fail("lexical descriptor index");
  return result;
}

std::pair<std::size_t, std::size_t> PathSampler::CatalogReader::equal_range(
    const std::vector<int>& path) const {
  PathCatalog::DecodeBuffer decoded;
  const auto less = [&](std::uint32_t index, const std::vector<int>& key) {
    const auto stations = catalog.stations(descriptor(index), decoded);
    return std::lexicographical_compare(stations.begin(), stations.end(), key.begin(), key.end());
  };
  const auto greater = [&](const std::vector<int>& key, std::uint32_t index) {
    const auto stations = catalog.stations(descriptor(index), decoded);
    return std::lexicographical_compare(key.begin(), key.end(), stations.begin(), stations.end());
  };
  if (!backing) {
    const auto first = std::lower_bound(catalog.lexical_order.begin(), catalog.lexical_order.end(), path, less);
    const auto last = std::upper_bound(first, catalog.lexical_order.end(), path, greater);
    return {static_cast<std::size_t>(first - catalog.lexical_order.begin()),
            static_cast<std::size_t>(last - catalog.lexical_order.begin())};
  }
  std::size_t low = 0, high = catalog.size();
  while (low < high) {
    const auto middle = low + (high - low) / 2;
    if (less(lexical_index(middle), path)) low = middle + 1;
    else high = middle;
  }
  const auto first = low;
  high = catalog.size();
  while (low < high) {
    const auto middle = low + (high - low) / 2;
    if (!greater(path, lexical_index(middle))) low = middle + 1;
    else high = middle;
  }
  return {first, low};
}

void PathSampler::CatalogReader::verify_all() const {
  if (!backing) return;
  auto& state = *backing->state;
  CatalogDiskTimer timer{state.read_s};
  // A deliberate full audit bypasses verified hot pages. Untouched disk
  // corruption is otherwise detected only when that chunk is first needed.
  struct stat status{};
  CatalogDiskCache::Header header{};
  if (::fstat(fd, &status) != 0 || status.st_size < 0 ||
      static_cast<std::uint64_t>(status.st_size) != backing->header[2] ||
      !catalog_pread(fd, header.data(), sizeof(header), 0, state.bytes_read) ||
      header != backing->header) fail("audit header/length");
  std::array<unsigned char, kCatalogDiskChunk> scratch;
  std::uint64_t hash = kCatalogHashOffset;
  for (std::uint64_t offset = sizeof(CatalogDiskCache::Header); offset < backing->header[2];) {
    const auto amount = static_cast<std::size_t>(
        std::min<std::uint64_t>(scratch.size(), backing->header[2] - offset));
    if (!catalog_pread(fd, scratch.data(), amount, offset, state.bytes_read)) fail("audit short read");
    catalog_hash_bytes(hash, scratch.data(), amount);
    offset += amount;
  }
  if (hash != backing->header[16]) fail("whole-payload checksum");
  if (::fstat(fd, &status) != 0 || status.st_size < 0 ||
      static_cast<std::uint64_t>(status.st_size) != backing->header[2])
    fail("audit final length");
  ++state.audits;
}

PathSampler::CatalogDiskCache::~CatalogDiskCache() {
  if (!requested) return;
  std::fprintf(stderr,
      "catalog_disk_cache enabled=%d hits=%llu misses=%llu stores=%llu "
      "failures=%llu quota_skips=%llu fp_skips=%llu bytes_read=%llu "
      "bytes_written=%llu charged_bytes=%llu byte_limit=%llu records=%zu "
      "read_s=%.9f write_s=%.9f\n",
      enabled ? 1 : 0, static_cast<unsigned long long>(hits),
      static_cast<unsigned long long>(misses), static_cast<unsigned long long>(stores),
      static_cast<unsigned long long>(failures),
      static_cast<unsigned long long>(quota_skips),
      static_cast<unsigned long long>(fp_skips),
      static_cast<unsigned long long>(bytes_read),
      static_cast<unsigned long long>(bytes_written),
      static_cast<unsigned long long>(charged_bytes),
      static_cast<unsigned long long>(byte_limit), records.size(), read_s, write_s);
  if (page_state) {
    const auto& p = *page_state;
    std::fprintf(stderr,
        "catalog_indexed_blocks hits=%llu misses=%llu bytes_read=%llu failures=%llu "
        "audits=%llu bytes=%zu peak_bytes=%zu byte_limit=%zu metadata_bytes=%llu "
        "metadata_limit=%llu active_readers=%zu peak_readers=%zu read_s=%.9f\n",
        static_cast<unsigned long long>(p.hits), static_cast<unsigned long long>(p.misses),
        static_cast<unsigned long long>(p.bytes_read), static_cast<unsigned long long>(p.failures),
        static_cast<unsigned long long>(p.audits), p.bytes, p.peak_bytes, p.budget,
        static_cast<unsigned long long>(metadata_bytes), static_cast<unsigned long long>(metadata_limit),
        p.active_readers, p.peak_readers, p.read_s);
  }
}

void PathSampler::initialize_catalog_disk_cache() {
  auto& disk = catalog_disk_cache_;
  catalog_cache_budget_ = static_cast<std::size_t>(
      catalog_budget("ILS_SP_CATALOG_RAM_BYTES", catalog_cache_budget_));
  const auto block_budget = catalog_budget("ILS_SP_CATALOG_BLOCK_BYTES", 64ULL * 1024 * 1024);
  disk.byte_limit = catalog_budget("ILS_SP_CATALOG_DISK_BYTES", disk.byte_limit);
  disk.metadata_limit = catalog_budget("ILS_SP_CATALOG_METADATA_BYTES", disk.metadata_limit);
  const char* requested = std::getenv("ILS_SP_CATALOG_DISK_DIR");
  if (requested == nullptr || *requested == '\0') return;
  disk.requested = true;
  if (XXH_versionNumber() != 803)
    throw std::runtime_error("indexed XXH3 checksum requires libxxhash 0.8.3");
  if (disk.byte_limit < 4096) { ++disk.quota_skips; return; }
  try {
    const std::filesystem::path input{requested};
    if (!input.is_absolute() || input.filename().empty()) {
      ++disk.failures; return;
    }
    std::error_code error;
    const auto parent = std::filesystem::canonical(input.parent_path(), error);
    if (error || !parent.string().starts_with("/workspace/")) {
      ++disk.failures; return;
    }
    struct statfs filesystem{};
    if (::statfs(parent.c_str(), &filesystem) != 0 ||
        filesystem.f_type == 0x01021994 || filesystem.f_type == 0x858458f6 ||
        filesystem.f_bsize != 4096 || filesystem.f_frsize != 4096) {
      ++disk.failures; return;  // tmpfs and ramfs are never a disk backing tier.
    }
    const auto space = std::filesystem::space(parent, error);
    if (error || space.available < disk.byte_limit + kCatalogDiskReserve) {
      ++disk.quota_skips; return;
    }
    const auto root = parent / input.filename();
    // Existing files/directories, including prior runs, are never reused,
    // overwritten or deleted. A second sampler with the same dir falls back.
    if (!std::filesystem::create_directory(root, error) || error) {
      ++disk.failures; return;
    }
    disk.root = root;
    disk.charged_bytes = 4096;  // Reserve the initial directory block too.
    disk.page_state = std::make_shared<CatalogPageState>();
    disk.page_state->budget = static_cast<std::size_t>(block_budget);
    disk.enabled = true;
  } catch (const std::exception&) {
    ++disk.failures;
  }
}

std::shared_ptr<const PathSampler::PathCatalog>
PathSampler::disk_cache_load(AnchorPair key) {
  auto& disk = catalog_disk_cache_;
  if (!disk.enabled) return {};
  if (!supported_bound_fp_environment()) { ++disk.fp_skips; return {}; }
  const auto found = disk.records.find(key);
  if (found == disk.records.end() || !found->second.complete) {
    ++disk.misses; return {};
  }
  const auto& expected = found->second.header;
  const std::uint64_t flags = static_cast<std::uint64_t>(complete_paths_) |
      (static_cast<std::uint64_t>(packed_catalogs_enabled_) << 1) |
      (static_cast<std::uint64_t>(packed_construction_enabled_) << 2) |
      (static_cast<std::uint64_t>(catalog_first_leg_index_enabled_) << 3);
  if (expected[5] != flags || expected[15] != catalog_fp_signature()) {
    ++disk.fp_skips; return {};
  }
  if (expected[6] != 0) {
    // The run-local manifest is trusted construction output. Creating a new
    // logical identity does not read or allocate the complete payload. The
    // operation-scoped reader validates the file and every used chunk.
    if (!found->second.backing) throw std::runtime_error("missing indexed catalogue backing");
    auto result = std::make_shared<PathCatalog>();
    result->packed = true;
    result->packed_station_count = static_cast<std::uint8_t>(expected[7]);
    for (std::size_t index = 0; index < result->packed_station_ids.size(); ++index)
      result->packed_station_ids[index] = std::bit_cast<int>(
          static_cast<std::uint32_t>(expected[17 + index]));
    result->disk_backing = found->second.backing;
    ++disk.hits;  // Indexed handle hits; physical page reads are separate.
    return result;
  }
  CatalogDiskTimer timer{disk.read_s};
  try {
    const auto name = catalog_disk_file(disk.root, key);
    std::error_code error;
    if (std::filesystem::file_size(name, error) != expected[2] || error) {
      ++disk.failures; return {};
    }
    CatalogDiskFile file{std::fopen(name.c_str(), "rb")};
    if (file.file == nullptr) { ++disk.failures; return {}; }
    CatalogDiskCache::Header header{};
    if (!catalog_transfer(file.file, header.data(), sizeof(header), false,
                          disk.bytes_read, nullptr) || header != expected) {
      ++disk.failures; return {};
    }
    // The in-memory manifest was created from this exact sampler's completed
    // output. Checking it before allocation validates schema, key, flags,
    // counts, capacities, packed alphabet and payload checksum together.
    auto result = std::make_shared<PathCatalog>();
    result->packed = header[6] != 0;
    result->packed_station_count = static_cast<std::uint8_t>(header[7]);
    for (std::size_t index = 0; index < result->packed_station_ids.size(); ++index)
      result->packed_station_ids[index] = std::bit_cast<int>(
          static_cast<std::uint32_t>(header[17 + index]));
    result->paths.reserve(static_cast<std::size_t>(header[11]));
    result->station_ids.reserve(static_cast<std::size_t>(header[12]));
    result->lexical_order.reserve(static_cast<std::size_t>(header[13]));
    result->paths.resize(static_cast<std::size_t>(header[8]));
    result->station_ids.resize(static_cast<std::size_t>(header[9]));
    result->lexical_order.resize(static_cast<std::size_t>(header[10]));
    std::uint64_t hash = kCatalogHashOffset;
    const bool okay = catalog_transfer(file.file, result->paths.data(),
            result->paths.size() * sizeof(PathDescriptor), false, disk.bytes_read, &hash) &&
        catalog_transfer(file.file, result->station_ids.data(),
            result->station_ids.size() * sizeof(int), false, disk.bytes_read, &hash) &&
        catalog_transfer(file.file, result->lexical_order.data(),
            result->lexical_order.size() * sizeof(std::uint32_t), false,
            disk.bytes_read, &hash);
    if (!okay || hash != header[16] || std::fgetc(file.file) != EOF ||
        std::ferror(file.file) != 0 || !file.close()) {
      ++disk.failures; return {};
    }
    ++disk.hits;
    return result;
  } catch (const std::exception&) {
    ++disk.failures;
    return {};
  }
}

void PathSampler::disk_cache_store(AnchorPair key, const PathCatalog& catalog) {
  auto& disk = catalog_disk_cache_;
  if (!disk.enabled || disk.records.contains(key)) return;
  if (!supported_bound_fp_environment()) { ++disk.fp_skips; return; }
  CatalogDiskTimer timer{disk.write_s};
  static_assert(sizeof(int) == sizeof(std::uint32_t));
  static_assert(std::is_trivially_copyable_v<PathDescriptor>);
  try {
    CatalogDiskCache::Header header{};
    header[0] = kCatalogDiskMagic;
    header[1] = kCatalogDiskVersion;
    header[3] = std::bit_cast<std::uint32_t>(key.origin);
    header[4] = std::bit_cast<std::uint32_t>(key.target);
    header[5] = static_cast<std::uint64_t>(complete_paths_) |
        (static_cast<std::uint64_t>(packed_catalogs_enabled_) << 1) |
        (static_cast<std::uint64_t>(packed_construction_enabled_) << 2) |
        (static_cast<std::uint64_t>(catalog_first_leg_index_enabled_) << 3);
    header[6] = catalog.packed;
    header[7] = catalog.packed_station_count;
    header[8] = catalog.paths.size();
    header[9] = catalog.station_ids.size();
    header[10] = catalog.lexical_order.size();
    header[11] = catalog.paths.capacity();
    header[12] = catalog.station_ids.capacity();
    header[13] = catalog.lexical_order.capacity();
    header[14] = kCatalogDiskChunk;
    header[15] = catalog_fp_signature();
    for (std::size_t index = 0; index < catalog.packed_station_ids.size(); ++index)
      header[17 + index] = std::bit_cast<std::uint32_t>(catalog.packed_station_ids[index]);
    std::uint64_t bytes = sizeof(header);
    const auto add_size = [&](std::uint64_t count, std::uint64_t width) {
      if (bytes > disk.byte_limit || count > (disk.byte_limit - bytes) / width)
        return false;
      bytes += count * width;
      return true;
    };
    if (!add_size(header[8], sizeof(PathDescriptor)) ||
        !add_size(header[9], sizeof(int)) ||
        !add_size(header[10], sizeof(std::uint32_t))) {
      ++disk.quota_skips; return;
    }
    header[2] = bytes;
    // Charge complete filesystem blocks plus per-file directory slack. Failed
    // partial writes retain their full reservation and are never retried.
    const std::uint64_t charged = ((bytes + 4095) / 4096) * 4096 + 4096;
    if (disk.records.size() >= kCatalogDiskMaximumRecords ||
        disk.charged_bytes > disk.byte_limit ||
        charged > disk.byte_limit - disk.charged_bytes) {
      ++disk.quota_skips; return;
    }
    std::error_code error;
    const auto space = std::filesystem::space(disk.root, error);
    if (error || space.available < charged + kCatalogDiskReserve) {
      ++disk.quota_skips; return;
    }
    const auto destination = catalog_disk_file(disk.root, key);
    const auto partial = std::filesystem::path{destination.string() + ".part"};
    if (std::filesystem::exists(destination, error) || error ||
        std::filesystem::exists(partial, error) || error) {
      ++disk.failures; return;
    }
    const std::uint64_t descriptor_bytes = header[8] * sizeof(PathDescriptor);
    const std::uint64_t lexical_bytes = header[10] * sizeof(std::uint32_t);
    const std::uint64_t block_count = catalog.packed
        ? (descriptor_bytes + kCatalogDiskChunk - 1) / kCatalogDiskChunk +
          (lexical_bytes + kCatalogDiskChunk - 1) / kCatalogDiskChunk : 0;
    // Fixed per-record allowance bounds manifest/container/path overhead; the
    // exact checksum vector capacity is charged below. Failed stores retain it.
    const std::uint64_t metadata_charge = 1024 + 2 * destination.native().size() +
        block_count * sizeof(CatalogDiskBacking::BlockChecksum);
    if (disk.metadata_bytes > disk.metadata_limit ||
        metadata_charge > disk.metadata_limit - disk.metadata_bytes) {
      ++disk.quota_skips; return;
    }
    auto [position, inserted] = disk.records.emplace(
        key, CatalogDiskCache::Record{header, false, {}});
    if (!inserted) return;
    disk.charged_bytes += charged;
    disk.metadata_bytes += metadata_charge;
    std::shared_ptr<CatalogDiskBacking> backing;
    if (catalog.packed) {
      backing = std::make_shared<CatalogDiskBacking>();
      backing->file = destination.native();
      if (backing->file.capacity() > 2 * destination.native().size() + 1)
        throw std::runtime_error("unexpected backing filename capacity");
      backing->state = disk.page_state;
      backing->record_id = disk.records.size();
      backing->block_hashes.reserve(static_cast<std::size_t>(block_count));
      if (backing->block_hashes.capacity() != block_count)
        throw std::runtime_error("unexpected checksum vector capacity");
    }
    auto* block_hashes = backing ? &backing->block_hashes : nullptr;
    CatalogDiskFile file{std::fopen(partial.c_str(), "wbx")};
    if (file.file == nullptr) { ++disk.failures; return; }
    std::uint64_t hash = kCatalogHashOffset;
    const bool okay = catalog_transfer(file.file, header.data(), sizeof(header), true,
            disk.bytes_written, nullptr) &&
        catalog_transfer(file.file, const_cast<PathDescriptor*>(catalog.paths.data()),
            catalog.paths.size() * sizeof(PathDescriptor), true, disk.bytes_written, &hash, block_hashes) &&
        catalog_transfer(file.file, const_cast<int*>(catalog.station_ids.data()),
            catalog.station_ids.size() * sizeof(int), true, disk.bytes_written, &hash) &&
        catalog_transfer(file.file, const_cast<std::uint32_t*>(catalog.lexical_order.data()),
            catalog.lexical_order.size() * sizeof(std::uint32_t), true,
            disk.bytes_written, &hash, block_hashes);
    header[16] = hash;
    if (!okay || std::fseek(file.file, 0, SEEK_SET) != 0 ||
        !catalog_transfer(file.file, header.data(), sizeof(header), true,
                          disk.bytes_written, nullptr) ||
        std::fflush(file.file) != 0 || !file.close()) {
      ++disk.failures; return;
    }
    // Atomic publication must also refuse to replace any unexpected existing
    // destination, even if one appeared after the earlier existence check.
    if (::syscall(SYS_renameat2, AT_FDCWD, partial.c_str(), AT_FDCWD,
                  destination.c_str(), RENAME_NOREPLACE) != 0) {
      ++disk.failures; return;
    }
    position->second.header = header;
    if (backing) {
      if (backing->block_hashes.size() != block_count)
        throw std::runtime_error("indexed checksum count differs");
      backing->header = header;
      position->second.backing = std::move(backing);
    }
    position->second.complete = true;
    ++disk.stores;
  } catch (const std::exception&) {
    ++disk.failures;
  }
}
