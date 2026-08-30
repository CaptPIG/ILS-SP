#pragma once

#include <filesystem>

#include "ils_sp/core.hpp"

namespace ils_sp {

[[nodiscard]] Instance load_montoya_instance(const std::filesystem::path& path);

}  // namespace ils_sp

