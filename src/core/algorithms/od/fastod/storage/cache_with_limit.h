#pragma once

#include <cstddef>
#include <queue>
#include <unordered_map>

namespace algos::fastod {

template <typename K, typename V>
class CacheWithLimit {
private:
    std::unordered_map<K, V> entries_;
    std::queue<K> keys_in_order_;
    const size_t max_size_;

public:
    explicit CacheWithLimit(size_t max_size) : max_size_(max_size){};

    void Clear() {
        entries_.clear();
        keys_in_order_ = {};
    }

    bool Contains(const K& key) const noexcept {
        return entries_.find(key) != entries_.end();
    }

    const V& Get(const K& key) const {
        return entries_.at(key);
    }

    void Set(const K& key, const V& value) {
        if (!entries_.try_emplace(key, value).second) return;

        if (keys_in_order_.size() >= max_size_) {
            entries_.erase(keys_in_order_.front());
            keys_in_order_.pop();
        }

        keys_in_order_.push(key);
    }
};

}  // namespace algos::fastod
