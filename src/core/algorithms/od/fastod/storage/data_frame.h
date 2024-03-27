#pragma once

#include <filesystem>
#include <optional>
#include <vector>

#include "algorithms/od/fastod/model/attribute_set.h"
#include "config/equal_nulls/type.h"
#include "config/tabular_data/input_table_type.h"
#include "table/column_layout_typed_relation_data.h"

namespace algos::fastod {

class DataFrame {
public:
    using Range = std::pair<size_t, size_t>;
    using ValueIndices = std::pair<int, Range>;
    using DataAndIndex = std::pair<std::byte const*, int>;

private:
    std::vector<std::vector<int>> data_;
    std::vector<std::vector<DataFrame::ValueIndices>> data_ranges_;
    std::vector<std::vector<size_t>> range_item_placement_;

    AttributeSet attrs_with_ranges_;

    void RecognizeAttributesWithRanges();

    static std::vector<std::pair<std::byte const*, int>> CreateIndexedColumnData(
            model::TypedColumnData const& column);

    static std::vector<int> ConvertColumnDataToIntegers(model::TypedColumnData const& column);

    static std::vector<DataFrame::ValueIndices> ExtractRangesFromColumn(
            std::vector<int> const& column);

    static std::optional<size_t> FindRangeIndexByItem(
            size_t item, std::vector<DataFrame::ValueIndices> const& ranges);

public:
    DataFrame(DataFrame const&) = delete;
    DataFrame& operator=(DataFrame const&) = delete;

    DataFrame() = default;
    DataFrame(DataFrame&&) noexcept = default;
    DataFrame& operator=(DataFrame&&) = default;
    ~DataFrame() noexcept = default;

    explicit DataFrame(std::vector<model::TypedColumnData> const& columns_data);

    int GetValue(int tuple_index, model::ColumnIndex attribute_index) const;
    std::vector<std::vector<DataFrame::ValueIndices>> const& GetDataRanges() const;
    size_t GetRangeIndexByItem(size_t item, model::ColumnIndex attribute) const;

    model::ColumnIndex GetColumnCount() const;
    size_t GetTupleCount() const;

    bool IsAttributesMostlyRangeBased(AttributeSet attributes) const;

    static DataFrame FromCsv(std::filesystem::path const& path, char separator = ',',
                             bool has_header = true, config::EqNullsType is_null_equal_null = true);

    static DataFrame FromInputTable(config::InputTable input_table,
                                    config::EqNullsType is_null_equal_null = true);
};

inline size_t RangeSize(DataFrame::Range const& range) {
    return range.second - range.first + 1;
}

}  // namespace algos::fastod
