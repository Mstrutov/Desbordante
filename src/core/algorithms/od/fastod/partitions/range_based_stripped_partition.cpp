// This part of code not currently used

#include "range_based_stripped_partition.h"

#include <sstream>

namespace algos::fastod {

RangeBasedStrippedPartition::RangeBasedStrippedPartition(DataFrame const& data)
    : data_(std::move(data)), should_be_converted_to_sp_(false) {
    const size_t tuple_count = data.GetTupleCount();
    begins_.push_back(0);

    if (tuple_count != 0) {
        indexes_.push_back({0, tuple_count - 1});
        begins_.push_back(1);
    }
}

RangeBasedStrippedPartition::RangeBasedStrippedPartition(
        DataFrame const& data, std::vector<DataFrame::Range> const& indexes,
        std::vector<size_t> const& begins)
    : indexes_(indexes), begins_(begins), data_(data), should_be_converted_to_sp_(false) {}

bool RangeBasedStrippedPartition::ShouldBeConvertedToStrippedPartition() const {
    return should_be_converted_to_sp_;
}

std::string RangeBasedStrippedPartition::ToString() const {
    std::stringstream result;
    std::string indexes_string;

    for (size_t i = 0; i < indexes_.size(); i++) {
        if (i != 0) {
            indexes_string += ", ";
        }

        indexes_string += "(" + std::to_string(indexes_[i].first) + ";" +
                          std::to_string(indexes_[i].second) + ")";
    }

    std::string begins_string;

    for (size_t i = 0; i < begins_.size(); i++) {
        if (i != 0) {
            begins_string += ", ";
        }

        begins_string += std::to_string(begins_[i]);
    }

    result << "RangeBasedStrippedPartition { indexes = [ " << indexes_string << " ]; begins = [ "
           << begins_string << " ] }";

    return result.str();
}

StrippedPartition RangeBasedStrippedPartition::ToStrippedPartition() const {
    std::vector<size_t> sp_indexes;
    std::vector<size_t> sp_begins;

    size_t sp_begin = 0;
    sp_begins.push_back(sp_begin);

    for (size_t begin_pointer = 0; begin_pointer < begins_.size() - 1; begin_pointer++) {
        const size_t group_begin = begins_[begin_pointer];
        const size_t group_end = begins_[begin_pointer + 1];

        for (size_t i = group_begin; i < group_end; ++i) {
            const DataFrame::Range range = indexes_[i];
            sp_begin += range.second - range.first + 1;

            for (size_t sp_index = range.first; sp_index <= range.second; ++sp_index) {
                sp_indexes.push_back(sp_index);
            }
        }

        sp_begins.push_back(sp_begin);
    }

    return StrippedPartition(data_, std::move(sp_indexes), std::move(sp_begins));
}

RangeBasedStrippedPartition& RangeBasedStrippedPartition::operator=(
        RangeBasedStrippedPartition const& other) {
    if (this == &other) {
        return *this;
    }

    indexes_ = other.indexes_;
    begins_ = other.begins_;

    return *this;
}

void RangeBasedStrippedPartition::Product(model::ColumnIndex attribute) {
    std::vector<size_t> new_begins;
    std::vector<DataFrame::Range> new_indexes;

    size_t curr_begin = 0;

    for (size_t begin_pointer = 0; begin_pointer < begins_.size() - 1; ++begin_pointer) {
        const size_t group_begin = begins_[begin_pointer];
        const size_t group_end = begins_[begin_pointer + 1];

        std::vector<DataFrame::ValueIndices> intersection =
                IntersectWithAttribute(attribute, group_begin, group_end - 1);

        const size_t intersection_size = intersection.size();
        size_t small_ranges_count = 0;

        auto add_group = [&new_indexes, &new_begins, &intersection, &curr_begin,
                          &small_ranges_count](size_t start_index, size_t end_index) {
            if (start_index == end_index) {
                DataFrame::Range range = intersection[start_index].second;

                if (range.second == range.first) {
                    return;
                }
            }

            for (size_t i = start_index; i <= end_index; ++i) {
                DataFrame::Range const& range = intersection[i].second;

                if (RangeSize(range) < kMinMeaningfulRangeSize) {
                    small_ranges_count++;
                }

                new_indexes.push_back(std::move(range));
            }

            new_begins.push_back(curr_begin);
            curr_begin += end_index - start_index + 1;
        };

        size_t group_start = 0;

        for (size_t i = 1; i < intersection_size; ++i) {
            if (intersection[i].first != intersection[i - 1].first) {
                add_group(group_start, i - 1);
                group_start = i;
            }
        }

        add_group(group_start, intersection_size - 1);

        if (intersection_size > 0 && small_ranges_count / static_cast<double>(intersection_size) >=
                                             kSmallRangesRatioToConvert) {
            should_be_converted_to_sp_ = true;
        }
    }

    new_begins.push_back(new_indexes.size());

    indexes_ = std::move(new_indexes);
    begins_ = std::move(new_begins);
}

bool RangeBasedStrippedPartition::Split(model::ColumnIndex right) const {
    for (size_t begin_pointer = 0; begin_pointer < begins_.size() - 1; ++begin_pointer) {
        const size_t group_begin = begins_[begin_pointer];
        const size_t group_end = begins_[begin_pointer + 1];

        int const group_value = data_.GetValue(indexes_[group_begin].first, right);

        for (size_t i = group_begin; i < group_end; ++i) {
            const DataFrame::Range range = indexes_[i];

            for (size_t j = range.first; j <= range.second; ++j) {
                if (data_.GetValue(j, right) != group_value) {
                    return true;
                }
            }
        }
    }

    return false;
}

bool RangeBasedStrippedPartition::Swap(model::ColumnIndex left, model::ColumnIndex right,
                                       bool ascending) const {
    for (size_t begin_pointer = 0; begin_pointer < begins_.size() - 1; begin_pointer++) {
        const size_t group_begin = begins_[begin_pointer];
        const size_t group_end = begins_[begin_pointer + 1];

        std::vector<std::pair<int, int>> values;

        for (size_t i = group_begin; i < group_end; ++i) {
            const DataFrame::Range range = indexes_[i];

            for (size_t j = range.first; j <= range.second; ++j) {
                values.push_back({data_.GetValue(j, left), data_.GetValue(j, right)});
            }
        }

        if (ascending) {
            std::sort(values.begin(), values.end(),
                      [](auto const& p1, auto const& p2) { return p1.first < p2.first; });
        } else {
            std::sort(values.begin(), values.end(),
                      [](auto const& p1, auto const& p2) { return p2.first < p1.first; });
        }

        size_t prev_group_max_index = 0;
        size_t current_group_max_index = 0;
        bool is_first_group = true;

        for (size_t i = 0; i < values.size(); i++) {
            auto const& first = values[i].first;
            auto const& second = values[i].second;

            if (i != 0 && values[i - 1].first != first) {
                is_first_group = false;
                prev_group_max_index = current_group_max_index;
                current_group_max_index = i;
            } else if (values[current_group_max_index].second <= second) {
                current_group_max_index = i;
            }

            if (!is_first_group && values[prev_group_max_index].second > second) {
                return true;
            }
        }
    }

    return false;
}

std::vector<DataFrame::ValueIndices> RangeBasedStrippedPartition::IntersectWithAttribute(
        model::ColumnIndex attribute, size_t group_start, size_t group_end) {
    std::vector<DataFrame::ValueIndices> result;
    std::vector<DataFrame::ValueIndices> const& attr_ranges = data_.GetDataRanges()[attribute];

    for (size_t i = group_start; i <= group_end; ++i) {
        DataFrame::Range const& range = indexes_[i];

        const size_t lower_bound_range_index = data_.GetRangeIndexByItem(range.first, attribute);
        const size_t upper_bound_range_index = data_.GetRangeIndexByItem(range.second, attribute);

        for (size_t j = lower_bound_range_index; j <= upper_bound_range_index; ++j) {
            DataFrame::ValueIndices const& attr_value_range = attr_ranges[j];
            DataFrame::Range const& attr_range = attr_value_range.second;

            const size_t start = std::max(range.first, attr_range.first);
            const size_t end = std::min(range.second, attr_range.second);

            if (start <= end) {
                result.push_back({attr_value_range.first, {start, end}});
            }
        }
    }

    std::sort(result.begin(), result.end(),
              [](auto const& x, auto const& y) { return x.first < y.first; });

    return result;
}

}  // namespace algos::fastod
