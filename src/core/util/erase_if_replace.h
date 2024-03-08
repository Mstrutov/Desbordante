#pragma once

#include <algorithm>
#include <vector>

namespace util {

template <typename T, typename UnaryFunc>
void EraseIfReplace(std::vector<T>& vec, UnaryFunc pred) {
    auto it = vec.begin();
    auto end_it = vec.end();
    it = std::find_if(it, end_it, pred);
    while (it != end_it) {
        *it = std::move(*--end_it);
        vec.pop_back();
        end_it = vec.end();
        it = std::find_if(it, end_it, pred);
    }
}

}  // namespace util
