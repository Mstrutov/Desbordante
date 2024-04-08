#pragma once

#include <algorithm>
#include <vector>

namespace util {

// The standard does not allow the std::find_if predicate to modify the element.
template <typename Iter, typename Pred>
Iter FindIf(Iter first, Iter last, Pred pred) {
    while (first != last && !pred(*first)) ++first;
    return first;
}

template <typename T, typename UnaryFunc>
void EraseIfReplace(std::vector<T>& vec, UnaryFunc pred) {
    auto it = vec.begin();
    auto end_it = vec.end();
    it = FindIf(it, end_it, pred);
    while (it != end_it) {
        *it = std::move(*--end_it);
        vec.pop_back();
        end_it = vec.end();
        it = FindIf(it, end_it, pred);
    }
}

}  // namespace util
