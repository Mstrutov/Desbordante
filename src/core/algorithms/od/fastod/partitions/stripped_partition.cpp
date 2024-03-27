// This part of code not currently used

#include "stripped_partition.h"

#include <sstream>

namespace algos::fastod {

StrippedPartition::StrippedPartition(DataFrame const& data) : data_(std::move(data)) {
    indexes_.reserve(data.GetTupleCount());

    for (size_t i = 0; i < data.GetTupleCount(); i++) {
        indexes_.push_back(i);
    }

    if (data.GetTupleCount() != 0) {
        begins_.push_back(0);
    }

    begins_.push_back(data.GetTupleCount());
}

StrippedPartition::StrippedPartition(DataFrame const& data, std::vector<size_t> const& indexes,
                                     std::vector<size_t> const& begins)
    : indexes_(indexes), begins_(begins), data_(data) {}

std::string StrippedPartition::ToString() const {
    std::stringstream result;
    std::string indexes_string;

    for (size_t i = 0; i < indexes_.size(); i++) {
        if (i != 0) {
            indexes_string += ", ";
        }

        indexes_string += std::to_string(indexes_[i]);
    }

    std::string begins_string;

    for (size_t i = 0; i < begins_.size(); i++) {
        if (i != 0) {
            begins_string += ", ";
        }

        begins_string += std::to_string(begins_[i]);
    }

    result << "StrippedPartition { indexes = [ " << indexes_string << " ]; begins = [ "
           << begins_string << " ] }";

    return result.str();
}

StrippedPartition& StrippedPartition::operator=(StrippedPartition const& other) {
    if (this == &other) {
        return *this;
    }

    indexes_ = other.indexes_;
    begins_ = other.begins_;

    return *this;
}

void StrippedPartition::Product(model::ColumnIndex attribute) {
    std::vector<size_t> new_indexes;
    new_indexes.reserve(data_.GetColumnCount());

    std::vector<size_t> new_begins;
    size_t fill_pointer = 0;

    for (size_t begin_pointer = 0; begin_pointer < begins_.size() - 1; begin_pointer++) {
        const size_t group_begin = begins_[begin_pointer];
        const size_t group_end = begins_[begin_pointer + 1];

        std::vector<std::pair<int, size_t>> values(group_end - group_begin);

        for (size_t i = group_begin; i < group_end; i++) {
            const size_t index = indexes_[i];
            values[i - group_begin] = {data_.GetValue(index, attribute), index};
        }

        std::sort(values.begin(), values.end(),
                  [](auto const& p1, auto const& p2) { return p1.first < p2.first; });

        size_t group_start = 0;
        size_t group_index = 1;

        auto add_group = [&new_begins, &new_indexes, &values, &fill_pointer, &group_start,
                          &group_index]() {
            size_t group_size = group_index - group_start;

            if (group_size > 1) {
                new_begins.push_back(fill_pointer);
                fill_pointer += group_size;

                for (size_t j = group_start; j < group_start + group_size; ++j) {
                    new_indexes.push_back(values[j].second);
                }
            }

            group_start = group_index;
        };

        for (; group_index < values.size(); ++group_index) {
            if (values[group_index - 1].first != values[group_index].first) {
                add_group();
            }
        }

        add_group();
    }

    indexes_ = std::move(new_indexes);
    begins_ = std::move(new_begins);

    begins_.push_back(indexes_.size());
}

bool StrippedPartition::Split(model::ColumnIndex right) const {
    for (size_t begin_pointer = 0; begin_pointer < begins_.size() - 1; begin_pointer++) {
        const size_t group_begin = begins_[begin_pointer];
        const size_t group_end = begins_[begin_pointer + 1];

        int const group_value = data_.GetValue(indexes_[group_begin], right);

        for (size_t i = group_begin + 1; i < group_end; i++) {
            if (data_.GetValue(indexes_[i], right) != group_value) {
                return true;
            }
        }
    }

    return false;
}

bool StrippedPartition::Swap(model::ColumnIndex left, model::ColumnIndex right,
                             bool ascending) const {
    for (size_t begin_pointer = 0; begin_pointer < begins_.size() - 1; begin_pointer++) {
        const size_t group_begin = begins_[begin_pointer];
        const size_t group_end = begins_[begin_pointer + 1];

        std::vector<std::pair<int, int>> values(group_end - group_begin);

        for (size_t i = group_begin; i < group_end; ++i) {
            const size_t index = indexes_[i];
            values[i - group_begin] = {data_.GetValue(index, left), data_.GetValue(index, right)};
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

            // values are sorted by "first"
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

}  // namespace algos::fastod
