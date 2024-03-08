#pragma once

#include <boost/dynamic_bitset.hpp>

namespace util {

template <typename ForwardIt>
boost::dynamic_bitset<> IndicesToBitset(ForwardIt begin, ForwardIt end, size_t num_columns) {
    boost::dynamic_bitset<> bitset(num_columns);
    for (auto it = begin; it != end; ++it) {
        bitset.set(*it);
    }
    return bitset;
}

template <typename Container>
boost::dynamic_bitset<> IndicesToBitset(Container const& indices, size_t num_columns) {
    return IndicesToBitset(indices.cbegin(), indices.cend(), num_columns);
}

template <typename UnaryFunction>
void ForEachIndex(boost::dynamic_bitset<> const& bitset, UnaryFunction func) {
    for (auto index = bitset.find_first(); index != boost::dynamic_bitset<>::npos;
         index = bitset.find_next(index)) {
        func(index);
    }
}

template <typename Index>
std::vector<Index> BitsetToIndices(boost::dynamic_bitset<> const& bitset) {
    std::vector<Index> indices;
    indices.reserve(bitset.count());
    ForEachIndex(bitset, [&](auto i) {
        assert(i <= std::numeric_limits<Index>::max());
        indices.push_back(static_cast<Index>(i));
    });
    return indices;
}

}  // namespace util
