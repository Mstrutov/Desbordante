// This part of code not currently used

#pragma once

#include <string>
#include <vector>

#include "algorithms/od/fastod/storage/data_frame.h"
#include "stripped_partition.h"

namespace algos::fastod {

class StrippedPartition;

class RangeBasedStrippedPartition {
private:
    std::vector<DataFrame::Range> indexes_;
    std::vector<size_t> begins_;
    DataFrame const& data_;
    bool should_be_converted_to_sp_;

    static constexpr inline double kSmallRangesRatioToConvert = 0.5;
    static constexpr inline size_t kMinMeaningfulRangeSize = static_cast<size_t>(10);

    RangeBasedStrippedPartition(DataFrame const& data, std::vector<DataFrame::Range> const& indexes,
                                std::vector<size_t> const& begins);

    std::vector<DataFrame::ValueIndices> IntersectWithAttribute(model::ColumnIndex attribute,
                                                                size_t group_start,
                                                                size_t group_end);

public:
    RangeBasedStrippedPartition() = delete;
    explicit RangeBasedStrippedPartition(DataFrame const& data);
    RangeBasedStrippedPartition(RangeBasedStrippedPartition const& origin) = default;

    bool ShouldBeConvertedToStrippedPartition() const;

    std::string ToString() const;
    StrippedPartition ToStrippedPartition() const;

    RangeBasedStrippedPartition& operator=(RangeBasedStrippedPartition const& other);

    void Product(model::ColumnIndex attribute);
    bool Split(model::ColumnIndex right) const;
    bool Swap(model::ColumnIndex left, model::ColumnIndex right, bool ascending) const;
};

}  // namespace algos::fastod
