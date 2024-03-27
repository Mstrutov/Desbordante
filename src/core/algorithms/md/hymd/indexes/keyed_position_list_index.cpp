#include "algorithms/md/hymd/indexes/keyed_position_list_index.h"

namespace algos::hymd::indexes {

ValueIdentifier KeyedPositionListIndex::AddNextValue(std::string value) {
    auto [it, is_value_new] = value_id_mapping_.try_emplace(std::move(value), next_value_id_);
    if (is_value_new) {
        clusters_.emplace_back();
        ++next_value_id_;
    }
    clusters_[it->second].emplace_back(cur_record_id_);
    ++cur_record_id_;
    return it->second;
}

}  // namespace algos::hymd::indexes
