// This part of code not currently used

#pragma once

#include <string>
#include <vector>

#include "range_based_stripped_partition.h"

namespace algos::fastod {

class RangeBasedStrippedPartition;

class StrippedPartition {
private:
    std::vector<size_t> indexes_;
    std::vector<size_t> begins_;
    DataFrame const& data_;

    StrippedPartition(DataFrame const& data, std::vector<size_t> const& indexes,
                      std::vector<size_t> const& begins);

    friend class RangeBasedStrippedPartition;

public:
    StrippedPartition() = delete;
    explicit StrippedPartition(DataFrame const& data);
    StrippedPartition(StrippedPartition const& origin) = default;

    StrippedPartition& operator=(StrippedPartition const& other);

    std::string ToString() const;

    void Product(model::ColumnIndex attribute);
    bool Split(model::ColumnIndex right) const;
    bool Swap(model::ColumnIndex left, model::ColumnIndex right, bool ascending) const;
};

}  // namespace algos::fastod
