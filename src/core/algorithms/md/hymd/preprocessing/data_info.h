#pragma once

#include <algorithm>
#include <cstddef>
#include <functional>
#include <memory>
#include <unordered_set>

#include "algorithms/md/hymd/indexes/keyed_position_list_index.h"
#include "algorithms/md/hymd/table_identifiers.h"
#include "model/index.h"
#include "model/types/type.h"

namespace algos::hymd::preprocessing {

class DataInfo {
private:
    using DataDeleter = std::function<void(std::byte*)>;

    std::size_t const elements_;
    std::size_t const type_size_;
    std::unordered_set<ValueIdentifier> const nulls_;
    std::unordered_set<ValueIdentifier> const empty_;
    std::unique_ptr<std::byte[], DataDeleter> const data_;

    static std::unique_ptr<std::byte[], DataDeleter> MakeDataPtr(
            std::unique_ptr<std::byte[]> data, model::Type::Destructor destructor,
            std::size_t const elements, std::size_t const type_size,
            std::unordered_set<ValueIdentifier> const& nulls,
            std::unordered_set<ValueIdentifier> const& empty);

public:
    DataInfo(std::size_t elements, std::size_t type_size, std::unordered_set<ValueIdentifier> nulls,
             std::unordered_set<ValueIdentifier> empty, model::Type::Destructor destructor,
             std::unique_ptr<std::byte[]> data);

    std::byte const* GetAt(ValueIdentifier const value_id) const noexcept {
        return &data_[value_id * type_size_];
    }

    std::size_t GetElementNumber() const noexcept {
        return elements_;
    }

    std::unordered_set<ValueIdentifier> const& GetNulls() const noexcept {
        return nulls_;
    }

    std::unordered_set<ValueIdentifier> const& GetEmpty() const noexcept {
        return empty_;
    }

    template <typename Comparator>
    std::vector<std::pair<std::byte*, ValueIdentifier>> MakeSortedInfo(
            Comparator const& comparator) const {
        std::vector<std::pair<std::byte*, ValueIdentifier>> sorted_info;
        sorted_info.reserve(elements_ - nulls_.size() - empty_.size());
        for (ValueIdentifier value_id = 0; value_id < elements_; ++value_id) {
            if (nulls_.contains(value_id) || empty_.contains(value_id)) continue;
            sorted_info.emplace_back(&data_[value_id * type_size_], value_id);
        }
        std::sort(sorted_info.begin(), sorted_info.end(),
                  [&comparator](std::byte const* left, std::byte const* right) {
                      return comparator(left, right);
                  });
        return sorted_info;
    }

    static std::shared_ptr<DataInfo> MakeFrom(indexes::KeyedPositionListIndex const& pli,
                                              model::Type const& type);
};

}  // namespace algos::hymd::preprocessing
