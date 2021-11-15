#pragma once

#include "RelationalSchema.h"
#include "Vertical.h"

class CustomHashing {
private:
    enum class BitsetHashingMethod {
        kTryConvertToUlong,
        kTrimAndConvertToUlong
    };

    static constexpr BitsetHashingMethod kDefaultHashingMethod =
#ifdef SAFE_VERTICAL_HASHING
        BitsetHashingMethod::kTrimAndConvertToUlong;
#else
        BitsetHashingMethod::kTryConvertToUlong;
#endif

    template <auto bitsetHashingMethod = kDefaultHashingMethod>
    static size_t BitsetHash(boost::dynamic_bitset<> const& bitset);

    friend std::hash<Vertical>;
    friend std::hash<Column>;
};

template <>
inline size_t CustomHashing::BitsetHash<CustomHashing::BitsetHashingMethod::kTryConvertToUlong>
        (boost::dynamic_bitset<> const& bitset) {

    return bitset.to_ulong();
}

template <>
inline size_t CustomHashing::BitsetHash<CustomHashing::BitsetHashingMethod::kTrimAndConvertToUlong>
        (boost::dynamic_bitset<> const& bitset) {

    boost::dynamic_bitset<> copyBitset = bitset;
    copyBitset.resize(std::numeric_limits<unsigned long>::digits);
    return copyBitset.to_ulong();
}

namespace std {
    template<>
    struct hash<Vertical> {
        size_t operator()(Vertical const& k) const {
            return CustomHashing::BitsetHash(k.getColumnIndicesRef());
        }
    };

    template<>
    struct hash<Column> {
        size_t operator()(Column const& k) const {
            return k.getIndex();
        }
    };

    template<>
    struct hash<std::shared_ptr<Vertical>> {
        size_t operator()(std::shared_ptr<Vertical> const& k) const {
            return std::hash<Vertical>()(*k);
        }
    };

    template<>
    struct hash<std::shared_ptr<Column>> {
        size_t operator()(std::shared_ptr<Column> const& k) const {
            return std::hash<Column>()(*k);
        }
    };

    template<class T>
    struct hash<std::pair<Vertical, T>> {
        size_t operator()(std::pair<Vertical, T> const& k) const {
            return std::hash<Vertical>()(k.first);
        }
    };

    template<class T>
    struct hash<std::pair<std::shared_ptr<Vertical>, T>> {
        size_t operator()(std::pair<std::shared_ptr<Vertical>, T> const& k) const {
            return std::hash<Vertical>()(*k.first);
        }
    };
}

