#include "data_frame.h"

#include <memory>
#include <stdexcept>

#include "algorithms/od/fastod/util/type_util.h"
#include "csv_parser/csv_parser.h"

namespace algos::fastod {

DataFrame::DataFrame(std::vector<model::TypedColumnData> const& columns_data) {
    model::ColumnIndex cols_num = columns_data.size();

    if (cols_num == 0) {
        throw std::invalid_argument("Number of columns must be greater than zero");
    }

    data_.reserve(cols_num);
    data_ranges_.reserve(cols_num);
    range_item_placement_.reserve(cols_num);

    std::transform(columns_data.cbegin(), columns_data.cend(), std::back_inserter(data_),
                   ConvertColumnDataToIntegers);

    std::transform(data_.cbegin(), data_.cend(), std::back_inserter(data_ranges_),
                   ExtractRangesFromColumn);

    for (size_t column = 0; column < cols_num; ++column) {
        const size_t tuple_count = data_[column].size();

        std::vector<size_t> curr_column_item_placement;
        curr_column_item_placement.reserve(tuple_count);

        for (size_t i = 0; i < tuple_count; ++i) {
            const std::optional<size_t> range_index = FindRangeIndexByItem(i, data_ranges_[column]);

            if (!range_index.has_value()) {
                throw std::logic_error("Range index not found");
            }

            curr_column_item_placement.push_back(range_index.value());
        }

        range_item_placement_.push_back(std::move(curr_column_item_placement));
    }

    RecognizeAttributesWithRanges();
}

int DataFrame::GetValue(int tuple_index, model::ColumnIndex attribute_index) const {
    return data_[attribute_index][tuple_index];
}

std::vector<std::vector<DataFrame::ValueIndices>> const& DataFrame::GetDataRanges() const {
    return data_ranges_;
}

size_t DataFrame::GetRangeIndexByItem(size_t item, model::ColumnIndex attribute) const {
    return range_item_placement_[attribute][item];
}

model::ColumnIndex DataFrame::GetColumnCount() const {
    return data_.size();
}

size_t DataFrame::GetTupleCount() const {
    return data_.size() > 0 ? data_[0].size() : 0;
}

bool DataFrame::IsAttributesMostlyRangeBased(AttributeSet attributes) const {
    if (!attributes.Any()) {
        return false;
    }

    AttributeSet remaining_attrs = Intersect(attrs_with_ranges_, attributes);

    model::ColumnIndex attrs_count = attributes.Count();
    model::ColumnIndex remaining_attrs_count = remaining_attrs.Count();

    double constexpr accept_range_based_partition_factor = 0.5;

    return static_cast<double>(remaining_attrs_count) / attrs_count >=
           accept_range_based_partition_factor;
}

DataFrame DataFrame::FromCsv(std::filesystem::path const& path, char separator, bool has_header,
                             config::EqNullsType is_null_equal_null) {
    std::shared_ptr<CSVParser> parser = std::make_shared<CSVParser>(path, separator, has_header);
    return FromInputTable(parser, is_null_equal_null);
}

DataFrame DataFrame::FromInputTable(config::InputTable input_table,
                                    config::EqNullsType is_null_equal_null) {
    std::vector<model::TypedColumnData> columns_data =
            model::CreateTypedColumnData(*input_table, is_null_equal_null);

    return DataFrame(std::move(columns_data));
}

void DataFrame::RecognizeAttributesWithRanges() {
    double constexpr accept_factor = 0.001;

    for (size_t i = 0; i < data_ranges_.size(); ++i) {
        const size_t items_count = data_[i].size();
        const size_t ranges_count = data_ranges_[i].size();

        if (static_cast<double>(ranges_count) / items_count >= accept_factor) {
            attrs_with_ranges_.Set(i, true);
        }
    }
}

std::vector<std::pair<std::byte const*, int>> DataFrame::CreateIndexedColumnData(
        model::TypedColumnData const& column) {
    const std::vector<std::byte const*> data = column.GetData();
    std::vector<std::pair<std::byte const*, int>> indexed_column_data(data.size());

    for (size_t i = 0; i < data.size(); ++i) {
        indexed_column_data[i] = std::make_pair(data[i], i);
    }

    return indexed_column_data;
}

std::vector<int> DataFrame::ConvertColumnDataToIntegers(model::TypedColumnData const& column) {
    auto less_mixed = [&column](DataAndIndex l, DataAndIndex r) {
        return CompareData<true>(l, r, column) == model::CompareResult::kLess;
    };

    auto less_not_mixed = [&column](DataAndIndex l, DataAndIndex r) {
        return CompareData<false>(l, r, column) == model::CompareResult::kLess;
    };

    auto equal_mixed = [&column](DataAndIndex l, DataAndIndex r) {
        return CompareData<true>(l, r, column) == model::CompareResult::kEqual;
    };

    auto equal_not_mixed = [&column](DataAndIndex l, DataAndIndex r) {
        return CompareData<false>(l, r, column) == model::CompareResult::kEqual;
    };

    std::vector<std::pair<std::byte const*, int>> indexed_column_data =
            CreateIndexedColumnData(column);

    bool const is_column_mixed = column.IsMixed();

    if (is_column_mixed) {
        std::sort(indexed_column_data.begin(), indexed_column_data.end(), less_mixed);
    } else {
        std::sort(indexed_column_data.begin(), indexed_column_data.end(), less_not_mixed);
    }

    std::vector<int> converted_column(indexed_column_data.size());

    int current_value = 0;
    converted_column[indexed_column_data[0].second] = current_value;

    for (size_t i = 1; i < indexed_column_data.size(); ++i) {
        std::pair<std::byte const*, int> const& prev = indexed_column_data[i - 1];
        std::pair<std::byte const*, int> const& curr = indexed_column_data[i];

        if (is_column_mixed) {
            converted_column[curr.second] =
                    equal_mixed(prev, curr) ? current_value : ++current_value;
        } else {
            converted_column[curr.second] =
                    equal_not_mixed(prev, curr) ? current_value : ++current_value;
        }
    }

    return converted_column;
}

std::vector<DataFrame::ValueIndices> DataFrame::ExtractRangesFromColumn(
        std::vector<int> const& column) {
    std::vector<ValueIndices> ranges;

    size_t start = 0;

    for (size_t i = 1; i < column.size(); ++i) {
        int const curr_value = column[i];
        int const prev_value = column[i - 1];

        if (curr_value != prev_value) {
            ranges.emplace_back(prev_value, std::pair(start, i - 1));
            start = i;
        }
    }

    ranges.emplace_back(column[column.size() - 1], std::pair(start, column.size() - 1));

    return ranges;
}

std::optional<size_t> DataFrame::FindRangeIndexByItem(
        size_t item, std::vector<DataFrame::ValueIndices> const& ranges) {
    auto iter = std::find_if(ranges.cbegin(), ranges.cend(), [item](auto const& p) {
        Range const& range = p.second;
        return item >= range.first && item <= range.second;
    });

    return (iter != ranges.cend()) ? std::optional<size_t>{iter - ranges.cbegin()} : std::nullopt;
}

}  // namespace algos::fastod
