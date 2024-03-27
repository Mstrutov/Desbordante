#pragma once

#include "algorithms/od/fastod/storage/data_frame.h"
#include "model/types/builtin.h"

namespace algos::fastod {

bool IsUnorderedType(model::TypeId type_id);

model::CompareResult CompareDataAsStrings(std::byte const* left, std::byte const* right,
                                          model::MixedType const* mixed_type);

template <bool IsColumnMixed>
model::CompareResult CompareData(DataFrame::DataAndIndex const& left,
                                 DataFrame::DataAndIndex const& right,
                                 model::TypedColumnData const& column) {
    const model::TypeId left_type_id = column.GetValueTypeId(left.second);
    const model::TypeId right_type_id = column.GetValueTypeId(right.second);

    bool const is_both_types_unordered =
            IsUnorderedType(left_type_id) && IsUnorderedType(right_type_id);

    bool const is_one_type_unordered =
            IsUnorderedType(left_type_id) || IsUnorderedType(right_type_id);

    if (is_both_types_unordered) {
        return model::CompareResult::kEqual;
    }

    if (is_one_type_unordered) {
        return IsUnorderedType(left_type_id) ? model::CompareResult::kLess
                                             : model::CompareResult::kGreater;
    }

    if constexpr (IsColumnMixed) {
        model::MixedType const* mixed_type = column.GetIfMixed();

        return left_type_id == right_type_id
                       ? mixed_type->Compare(left.first, right.first)
                       : CompareDataAsStrings(left.first, right.first, mixed_type);
    } else {
        model::Type const& type = column.GetType();
        return type.Compare(left.first, right.first);
    }
}

}  // namespace algos::fastod