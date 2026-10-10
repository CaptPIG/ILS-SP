#pragma once

#include <cstddef>
#include <cstdint>
#include <type_traits>

// Minimal declaration of the installed xxHash 0.8.3 stable public C ABI.
// Verified against https://github.com/Cyan4973/xxHash/blob/v0.8.3/xxhash.h
// This file supplies no hashing implementation and installs no dependency.
// Digests are run-local native integers, never serialized into the v2 format.
extern "C" {
struct XXH128_hash_t {
  std::uint64_t low64;
  std::uint64_t high64;
};
unsigned XXH_versionNumber(void);
XXH128_hash_t XXH3_128bits(const void* data, std::size_t len);
}

static_assert(sizeof(XXH128_hash_t) == 16);
static_assert(offsetof(XXH128_hash_t, low64) == 0);
static_assert(offsetof(XXH128_hash_t, high64) == 8);
static_assert(alignof(XXH128_hash_t) == alignof(std::uint64_t));
static_assert(std::is_trivially_copyable_v<XXH128_hash_t>);
static_assert(std::is_standard_layout_v<XXH128_hash_t>);
