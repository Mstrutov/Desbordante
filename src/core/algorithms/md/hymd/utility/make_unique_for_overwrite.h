#pragma once

#include <cstddef>
#include <memory>
#include <type_traits>

namespace algos::hymd::utility {

// https://en.cppreference.com/w/cpp/memory/unique_ptr/make_unique
template <class T>
std::enable_if_t<std::is_unbounded_array_v<T>, std::unique_ptr<T>> MakeUniqueForOverwrite(
        std::size_t n) {
    return std::unique_ptr<T>(new std::remove_extent_t<T>[n]);
}

}  // namespace algos::hymd::utility
