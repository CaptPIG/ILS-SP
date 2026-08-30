#pragma once

#include <cstddef>
#include <list>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace ils_sp {

template <typename Key, typename Value, typename Hash = std::hash<Key>,
          typename Equal = std::equal_to<Key>>
class LruCache {
 public:
  explicit LruCache(std::size_t capacity) : capacity_(capacity) {
    if (capacity_ == 0) {
      throw std::invalid_argument("LRU capacity must be positive");
    }
  }

  [[nodiscard]] std::optional<Value> get(const Key& key) {
    const auto found = index_.find(key);
    if (found == index_.end()) {
      ++misses_;
      return std::nullopt;
    }
    order_.splice(order_.end(), order_, found->second);
    ++hits_;
    return found->second->second;
  }

  void put(Key key, Value value) {
    const auto found = index_.find(key);
    if (found != index_.end()) {
      found->second->second = std::move(value);
      order_.splice(order_.end(), order_, found->second);
      return;
    }
    order_.emplace_back(std::move(key), std::move(value));
    const auto item = std::prev(order_.end());
    index_.emplace(item->first, item);
    if (index_.size() > capacity_) {
      index_.erase(order_.front().first);
      order_.pop_front();
      ++evictions_;
    }
  }

  [[nodiscard]] bool contains(const Key& key) const {
    return index_.contains(key);
  }

  void clear() {
    order_.clear();
    index_.clear();
  }

  [[nodiscard]] std::size_t size() const noexcept { return index_.size(); }
  [[nodiscard]] std::size_t capacity() const noexcept { return capacity_; }
  [[nodiscard]] std::size_t hits() const noexcept { return hits_; }
  [[nodiscard]] std::size_t misses() const noexcept { return misses_; }
  [[nodiscard]] std::size_t evictions() const noexcept { return evictions_; }

 private:
  using Item = std::pair<Key, Value>;
  using Iterator = typename std::list<Item>::iterator;

  std::size_t capacity_;
  std::list<Item> order_;
  std::unordered_map<Key, Iterator, Hash, Equal> index_;
  std::size_t hits_{0};
  std::size_t misses_{0};
  std::size_t evictions_{0};
};

}  // namespace ils_sp

