#pragma once

#include <bitset>
#include <functional>
#include <stdexcept>
#include <string>

#include <boost/functional/hash.hpp>

#include "model/table/column_index.h"

namespace algos::fastod {

class AttributeSet {
private:
    static constexpr model::ColumnIndex kBitsNum = 64;

    std::bitset<kBitsNum> bitset_;

    explicit AttributeSet(std::bitset<kBitsNum> bitset) noexcept : bitset_(std::move(bitset)) {}

public:
    AttributeSet() noexcept = default;

    explicit AttributeSet([[maybe_unused]] model::ColumnIndex attribute_count) {
        if (attribute_count >= kBitsNum) {
            throw std::invalid_argument("Maximum possible number of attributes is " +
                                        std::to_string(kBitsNum - 1));
        }
    }

    explicit AttributeSet([[maybe_unused]] model::ColumnIndex attribute_count,
                          model::ColumnIndex value)
        : bitset_(value) {
        if (attribute_count >= kBitsNum) {
            throw std::invalid_argument("Maximum possible number of attributes is " +
                                        std::to_string(kBitsNum - 1));
        }
    }

    AttributeSet& operator&=(AttributeSet const& b) noexcept {
        bitset_ &= b.bitset_;
        return *this;
    }

    AttributeSet& operator|=(AttributeSet const& b) noexcept {
        bitset_ |= b.bitset_;
        return *this;
    }

    AttributeSet& operator^=(AttributeSet const& b) noexcept {
        bitset_ ^= b.bitset_;
        return *this;
    }

    AttributeSet operator~() const noexcept {
        AttributeSet as(~bitset_);
        return as;
    }

    AttributeSet& Set(model::ColumnIndex n, bool value = true) {
        bitset_.set(n, value);
        return *this;
    }

    AttributeSet& Reset(model::ColumnIndex n) {
        bitset_.reset(n);
        return *this;
    }

    bool Test(model::ColumnIndex n) const noexcept {
        return bitset_.test(n);
    }

    bool All() const noexcept {
        return bitset_.all();
    }

    bool Any() const noexcept {
        return bitset_.any();
    }

    bool None() const noexcept {
        return bitset_.none();
    }

    model::ColumnIndex Count() const noexcept {
        return bitset_.count();
    }

    model::ColumnIndex Size() const noexcept {
        return bitset_.size();
    }

    model::ColumnIndex FindFirst() const noexcept {
        return bitset_._Find_first();
    }

    model::ColumnIndex FindNext(model::ColumnIndex pos) const noexcept {
        return bitset_._Find_next(pos);
    }

    std::string ToString() const;
    void Iterate(std::function<void(model::ColumnIndex)> callback) const;

    friend AttributeSet operator&(AttributeSet const& b1, AttributeSet const& b2) noexcept;
    friend AttributeSet operator|(AttributeSet const& b1, AttributeSet const& b2) noexcept;
    friend AttributeSet operator^(AttributeSet const& b1, AttributeSet const& b2) noexcept;
    friend bool operator==(AttributeSet const& b1, AttributeSet const& b2) noexcept;
    friend bool operator!=(AttributeSet const& b1, AttributeSet const& b2) noexcept;
    friend bool operator<(AttributeSet const& b1, AttributeSet const& b2) noexcept;

    friend struct std::hash<AttributeSet>;
    friend struct boost::hash<AttributeSet>;
};

inline AttributeSet operator&(AttributeSet const& b1, AttributeSet const& b2) noexcept {
    AttributeSet as(b1.bitset_ & b2.bitset_);
    return as;
}

inline AttributeSet operator|(AttributeSet const& b1, AttributeSet const& b2) noexcept {
    AttributeSet as(b1.bitset_ | b2.bitset_);
    return as;
}

inline AttributeSet operator^(AttributeSet const& b1, AttributeSet const& b2) noexcept {
    AttributeSet as(b1.bitset_ ^ b2.bitset_);
    return as;
}

inline bool operator==(AttributeSet const& b1, AttributeSet const& b2) noexcept {
    return b1.bitset_ == b2.bitset_;
}

inline bool operator!=(AttributeSet const& b1, AttributeSet const& b2) noexcept {
    return !(b1 == b2);
}

inline bool operator<(AttributeSet const& b1, AttributeSet const& b2) noexcept {
    return b1.bitset_.to_ullong() < b2.bitset_.to_ullong();
}

}  // namespace algos::fastod

template <>
struct std::hash<algos::fastod::AttributeSet> {
    size_t operator()(algos::fastod::AttributeSet const& x) const noexcept {
        return x.bitset_.to_ullong();
    }
};

template <>
struct boost::hash<algos::fastod::AttributeSet> {
    size_t operator()(algos::fastod::AttributeSet const& x) const noexcept {
        return x.bitset_.to_ullong();
    }
};

namespace algos::fastod {

inline AttributeSet CreateAttributeSet(std::initializer_list<model::ColumnIndex> attributes,
                                       model::ColumnIndex size) {
    AttributeSet attr_set(size);

    for (auto const attr : attributes) {
        attr_set.Set(attr);
    }

    return attr_set;
}

inline bool ContainsAttribute(AttributeSet const& value, model::ColumnIndex attribute) noexcept {
    return value.Test(attribute);
}

inline AttributeSet AddAttribute(AttributeSet const& value, model::ColumnIndex attribute) {
    auto value_copy = value;
    return value_copy.Set(attribute);
}

inline AttributeSet DeleteAttribute(AttributeSet const& value, model::ColumnIndex attribute) {
    auto value_copy = value;
    return value_copy.Reset(attribute);
}

inline AttributeSet Intersect(AttributeSet const& value1, AttributeSet const& value2) noexcept {
    return value1 & value2;
}

inline AttributeSet Difference(AttributeSet const& value1, AttributeSet const& value2) noexcept {
    return value1 & (~value2);
}

inline bool IsEmptySet(AttributeSet const& value) noexcept {
    return value.None();
}

inline model::ColumnIndex GetAttributeCount(AttributeSet const& value) noexcept {
    return value.Count();
}

}  // namespace algos::fastod
