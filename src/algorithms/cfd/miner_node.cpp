#include "algorithms/cfd/miner_node.h"

// see ./LICENSE
namespace algos {

int support(const SimpleTidList& tids) {
    return tids.size();
}

int support(const PartitionTidList& tids) {
    if (tids.tids.size()) return tids.tids.size() + 1 - tids.sets_number;
    return 0;
}

int hash(const SimpleTidList& tids) {
    return std::accumulate(tids.begin(), tids.end(), 1);
}

int hash(const PartitionTidList& tids) {
    return std::accumulate(tids.tids.begin(), tids.tids.end(), 1) + (tids.sets_number - 1);
}
} //namespace algos
