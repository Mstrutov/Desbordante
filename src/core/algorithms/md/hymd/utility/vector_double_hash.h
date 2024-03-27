#pragma once

//#include <bit>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "algorithms/md/hymd/utility/bit_cast.h"
#include "algorithms/md/hymd/utility/java_hash.h"
#include "util/py_tuple_hash.h"

namespace std {
template <>
struct hash<std::vector<double>> {
    std::size_t operator()(std::vector<double> const& p) const noexcept {
        constexpr bool use_java_hash = true;
        if constexpr (use_java_hash) {
            return algos::hymd::utility::HashIterable(p, [](double element) {
                return /* TODO: replace with std::bit_cast when GCC in CI is upgraded */ algos::
                        hymd::utility::BitCast<std::int64_t>(element);
            });
        } else {
            util::PyTupleHash<double> hasher{p.size()};
            for (double el : p) {
                hasher.AddValue(el);
            }
            return hasher.GetResult();
        }
    }
};
}  // namespace std
