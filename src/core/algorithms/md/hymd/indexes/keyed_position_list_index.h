#pragma once

#include <cassert>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

#include "algorithms/md/hymd/indexes/pli_cluster.h"
#include "algorithms/md/hymd/table_identifiers.h"

namespace algos::hymd::indexes {

class KeyedPositionListIndex {
private:
    std::unordered_map<std::string, ValueIdentifier> value_id_mapping_;
    std::vector<PliCluster> clusters_;
    RecordIdentifier cur_record_id_ = 0;
    ValueIdentifier next_value_id_ = 0;

public:
    ValueIdentifier AddNextValue(std::string value);

    std::unordered_map<std::string, ValueIdentifier> const& GetMapping() const noexcept {
        assert(value_id_mapping_.size() == clusters_.size());
        return value_id_mapping_;
    }

    std::vector<PliCluster> const& GetClusters() const noexcept {
        return clusters_;
    }
};

}  // namespace algos::hymd::indexes
