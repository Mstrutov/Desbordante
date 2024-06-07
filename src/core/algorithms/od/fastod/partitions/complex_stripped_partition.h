#pragma once

#include <memory>
#include <string>
#include <vector>

#include "algorithms/od/fastod/storage/data_frame.h"

namespace algos::fastod {

class ComplexStrippedPartition {
protected:
    std::shared_ptr<std::vector<size_t>> sp_indexes_;
    std::shared_ptr<std::vector<size_t>> sp_begins_;
    std::shared_ptr<std::vector<DataFrame::Range>> rb_indexes_;
    std::shared_ptr<std::vector<size_t>> rb_begins_;
    std::shared_ptr<DataFrame> data_;
    bool is_stripped_partition_;
    bool should_be_converted_to_sp_;

    static constexpr inline double kSmallRangesRatioToConvert = 0.5;
    static constexpr inline size_t kMinMeaningfulRangeSize = static_cast<size_t>(40);

    std::string CommonToString() const;
    void CommonProduct(model::ColumnIndex attribute);
    bool CommonSplit(model::ColumnIndex right) const;

    std::string RangeBasedToString() const;
    void RangeBasedProduct(model::ColumnIndex attribute);
    bool RangeBasedSplit(model::ColumnIndex right) const;

    std::vector<DataFrame::ValueIndices> IntersectWithAttribute(model::ColumnIndex attribute,
                                                                size_t group_start,
                                                                size_t group_end);

    ComplexStrippedPartition(std::shared_ptr<DataFrame> data,
                             std::shared_ptr<std::vector<size_t>> indexes,
                             std::shared_ptr<std::vector<size_t>> begins);

    ComplexStrippedPartition(std::shared_ptr<DataFrame> data,
                             std::shared_ptr<std::vector<DataFrame::Range>> indexes,
                             std::shared_ptr<std::vector<size_t>> begins);

public:
    ComplexStrippedPartition();
    ComplexStrippedPartition(ComplexStrippedPartition const& origin) = default;

    ComplexStrippedPartition& operator=(ComplexStrippedPartition const& other);

    std::string ToString() const;
    void Product(model::ColumnIndex attribute);
    bool Split(model::ColumnIndex right) const;

    bool ShouldBeConvertedToStrippedPartition() const;
    void ToStrippedPartition();

    template <bool Ascending>
    bool Swap(model::ColumnIndex left, model::ColumnIndex right) const {
        const size_t group_count = is_stripped_partition_ ? sp_begins_->size() : rb_begins_->size();

        for (size_t begin_pointer = 0; begin_pointer < group_count - 1; begin_pointer++) {
            const size_t group_begin = is_stripped_partition_ ? (*sp_begins_)[begin_pointer]
                                                              : (*rb_begins_)[begin_pointer];

            const size_t group_end = is_stripped_partition_ ? (*sp_begins_)[begin_pointer + 1]
                                                            : (*rb_begins_)[begin_pointer + 1];

            std::vector<std::pair<int, int>> values;

            if (is_stripped_partition_) {
                values.reserve(group_end - group_begin);

                for (size_t i = group_begin; i < group_end; ++i) {
                    const size_t index = (*sp_indexes_)[i];

                    values.emplace_back(data_->GetValue(index, left),
                                        data_->GetValue(index, right));
                }
            } else {
                for (size_t i = group_begin; i < group_end; ++i) {
                    const DataFrame::Range range = (*rb_indexes_)[i];

                    for (size_t j = range.first; j <= range.second; ++j) {
                        values.emplace_back(data_->GetValue(j, left), data_->GetValue(j, right));
                    }
                }
            }

            if constexpr (Ascending) {
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
                auto const& [first, second] = values[i];

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

    template <bool RangeBasedMode>
    static ComplexStrippedPartition Create(std::shared_ptr<DataFrame> data) {
        if constexpr (RangeBasedMode) {
            auto rb_indexes = std::make_unique<std::vector<DataFrame::Range>>();
            auto rb_begins = std::make_unique<std::vector<size_t>>();

            const size_t tuple_count = data->GetTupleCount();
            rb_begins->push_back(0);

            if (tuple_count != 0) {
                rb_indexes->push_back({0, tuple_count - 1});
                rb_begins->push_back(1);
            }

            return ComplexStrippedPartition(std::move(data), std::move(rb_indexes),
                                            std::move(rb_begins));
        }

        auto sp_indexes = std::make_unique<std::vector<size_t>>();
        auto sp_begins = std::make_unique<std::vector<size_t>>();

        sp_indexes->reserve(data->GetTupleCount());

        for (size_t i = 0; i < data->GetTupleCount(); i++) {
            sp_indexes->push_back(i);
        }

        if (data->GetTupleCount() != 0) {
            sp_begins->push_back(0);
        }

        sp_begins->push_back(data->GetTupleCount());

        return ComplexStrippedPartition(std::move(data), std::move(sp_indexes),
                                        std::move(sp_begins));
    }
};

}  // namespace algos::fastod
