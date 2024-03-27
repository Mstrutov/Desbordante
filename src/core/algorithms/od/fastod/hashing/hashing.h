#pragma once

#include <cstddef>
#include <vector>

namespace algos::fastod::hashing {

inline size_t CombineHashes(size_t first, size_t second) {
    const size_t wave = second + 2654435769UL + (first << 6) + (first >> 2);
    return first ^ wave;
}

template <class T>
inline size_t CombineHashes(std::vector<T> const& items) {
    size_t result_hash = 0;

    for (auto const& item : items) {
        const size_t item_hash = std::hash<T>{}(item);
        result_hash = algos::fastod::hashing::CombineHashes(result_hash, item_hash);
    }

    return result_hash;
}

}  // namespace algos::fastod::hashing