#include "algorithms/md/hymd/preprocessing/data_info.h"

#include <cstddef>

namespace algos::hymd::preprocessing {

auto DataInfo::MakeDataPtr(std::unique_ptr<std::byte[]> data, model::Type::Destructor destructor,
                           std::size_t const elements, std::size_t const type_size,
                           std::unordered_set<ValueIdentifier> const& nulls,
                           std::unordered_set<ValueIdentifier> const& empty)
        -> std::unique_ptr<std::byte[], DataDeleter> {
    DataDeleter deleter = [destructor = std::move(destructor), elements, type_size, &nulls,
                           &empty](std::byte* data) {
        if (destructor) {
            for (ValueIdentifier value_id = 0; value_id < elements; ++value_id) {
                if (nulls.contains(value_id) || empty.contains(value_id)) continue;
                destructor(data + value_id * type_size);
            }
        }
        delete[] data;
    };
    std::unique_ptr<std::byte[], DataDeleter> ptr{nullptr, std::move(deleter)};
    ptr.reset(data.release());
    return ptr;
}

DataInfo::DataInfo(std::size_t elements, std::size_t type_size,
                   std::unordered_set<ValueIdentifier> nulls,
                   std::unordered_set<ValueIdentifier> empty, model::Type::Destructor destructor,
                   std::unique_ptr<std::byte[]> data)
    : elements_(elements),
      type_size_(type_size),
      nulls_(std::move(nulls)),
      empty_(std::move(empty)),
      data_(MakeDataPtr(std::move(data), std::move(destructor), elements_, type_size_, nulls_,
                        empty_)) {}

std::shared_ptr<DataInfo> DataInfo::MakeFrom(indexes::KeyedPositionListIndex const& pli,
                                             model::Type const& type) {
    std::size_t const value_number = pli.GetClusters().size();
    std::size_t const type_size = type.GetSize();
    auto data = std::unique_ptr<std::byte[]>(type.Allocate(value_number));
    std::unordered_set<ValueIdentifier> nulls;
    std::unordered_set<ValueIdentifier> empty;
    for (auto const& [string, value_id] : pli.GetMapping()) {
        if (string.empty()) {
            empty.insert(value_id);
            continue;
        }
        if (string == model::Null::kValue) {
            nulls.insert(value_id);
            continue;
        }
        type.ValueFromStr(&data[value_id * type_size], string);
    }
    return std::make_shared<DataInfo>(value_number, type_size, std::move(nulls), std::move(empty),
                                      type.GetDestructor(), std::move(data));
}

}  // namespace algos::hymd::preprocessing
