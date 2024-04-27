#pragma once

#include <cstddef>
#include <vector>

namespace util {
template <typename T>
std::vector<T> PickMHighestBias(std::vector<T> const& initial, std::size_t const size_limit) {
    std::size_t const initial_size = initial.size();
    if (initial_size <= size_limit) return initial;
    std::vector<T> elements;
    elements.reserve(size_limit);
    std::size_t const add_index = initial_size / size_limit;
    std::size_t const add_rem = initial_size % size_limit;
    std::size_t rem = add_rem;
    std::size_t index = add_index;
    while (index <= initial_size) {
        elements.push_back(initial[index - 1]);
        index += add_index;
        rem += add_rem;
        if (rem >= size_limit) {
            ++index;
            rem -= size_limit;
        }
    }
    return elements;
}
}  // namespace util
