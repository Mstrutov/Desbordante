#include "order.h"

#include <algorithm>

#include "config/names_and_descriptions.h"
#include "config/tabular_data/input_table/option.h"
#include "model/types/types.h"

namespace algos::order {

Order::Order() : Algorithm({}) {
    RegisterOptions();
    MakeOptionsAvailable({config::kTableOpt.GetName()});
}

void Order::RegisterOptions() {
    using namespace config::names;
    using namespace config::descriptions;
    using config::Option;

    RegisterOption(config::kTableOpt(&input_table_));
}

void Order::LoadDataInternal() {
    typed_relation_ = model::ColumnLayoutTypedRelationData::CreateFrom(*input_table_, false);
}

void Order::ResetState() {}

void Order::CreateSortedPartitions() {
    std::vector<model::TypedColumnData> const& data = typed_relation_->GetColumnData();
    for (unsigned int i = 0; i < data.size(); ++i) {
        if (!model::Type::IsOrdered(data[i].GetTypeId())) {
            continue;
        }
        std::vector<std::pair<unsigned long, std::byte const*>> indexed_byte_data;
        indexed_byte_data.reserve(data[i].GetNumRows());
        std::vector<std::byte const*> const& byte_data = data[i].GetData();
        for (size_t k = 0; k < byte_data.size(); ++k) {
            indexed_byte_data.emplace_back(k, byte_data[k]);
        }
        model::Type const& type = data[i].GetType();
        auto less = [&type](std::pair<unsigned long, std::byte const*> l,
                            std::pair<unsigned long, std::byte const*> r) {
            return type.Compare(l.second, r.second) == model::CompareResult::kLess;
        };
        std::sort(indexed_byte_data.begin(), indexed_byte_data.end(), less);
        std::vector<std::unordered_set<unsigned long>> equivalence_classes;
        equivalence_classes.push_back({indexed_byte_data.front().first});
        auto equal = [&type](std::pair<unsigned long, std::byte const*> l,
                             std::pair<unsigned long, std::byte const*> r) {
            return type.Compare(l.second, r.second) == model::CompareResult::kEqual;
        };
        for (size_t k = 1; k < indexed_byte_data.size(); ++k) {
            if (equal(indexed_byte_data[k - 1], indexed_byte_data[k])) {
                equivalence_classes.back().insert(indexed_byte_data[k].first);
            } else {
                equivalence_classes.push_back({indexed_byte_data[k].first});
            }
        }
        sorted_partitions_[{i}] = SortedPartition(std::move(equivalence_classes));
    }
}

unsigned long long Order::ExecuteInternal() {}

}  // namespace algos::order
