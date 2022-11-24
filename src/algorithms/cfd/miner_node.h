#pragma once

// see ./LICENSE

#include <numeric>

#include "model/cfd_types.h"
namespace algos {

int support(const SimpleTidList& tids);
int support(const PartitionTidList& tids);
int hash(const SimpleTidList& tids);
int hash(const PartitionTidList& tids);

template<typename T>
struct MinerNode {
    MinerNode()
        : item(-1), node_supp(-1), node_hash(-1) {

    }

    MinerNode(const Item& item, int supp)
        : item(item), node_supp(supp), node_hash(0) {
            //tidmap.reserve(Supp);
    }

    MinerNode(const Item& item)
        : item(item), node_supp(0), node_hash(0) {

    }

    MinerNode(const Item& item, const T& tids)
        : item(item), tids(tids), node_supp(0), node_hash(0) {
        HashTids();
        node_supp = Supp();
    }

    MinerNode(const Item& item, const T& tids, int supp)
        : item(item), tids(tids), node_supp(supp), node_hash(0), prefix() {
        HashTids();
    }

    MinerNode(const Item& item, const T& tids, int supp, const Itemset& prefix)
        : item(item), tids(tids), node_supp(supp), node_hash(0), prefix(prefix) {
        HashTids();
    }

    MinerNode(const Item& item, const T& tids, const Itemset& prefix)
        : item(item), tids(tids), node_supp(0), node_hash(0), prefix(prefix) {
        HashTids();
        node_supp = Supp();
    }

    MinerNode(const Item& item, const T& tids, int supp, int hash)
        : item(item), tids(tids), node_supp(supp), node_hash(hash) {

    }

    bool operator< (const MinerNode<T>& rhs) const {
        //return node_supp() < rhs.node_supp() || (node_supp() == rhs.Supp() && tidmap.size() > rhs.tidmap.size());
        //return cands.size() > rhs.cands.size();
        return LessThan(tids, rhs.tids);
    }

    bool operator== (const MinerNode<T>& rhs) const {
        return item == rhs.item && prefix == rhs.prefix && Equals(tids, rhs.tids);
    }

    int Supp() const {
        if (node_supp) return node_supp;
        return support(tids);
    }

    int Resupp() {
        node_supp = support(tids);
        return node_supp;
    }

    void HashTids() {
        node_hash = hash(tids);
    }

    Item item;
    T tids;
    int node_supp;
    int node_hash;
    Itemset prefix;
    Itemset cands;
};

typedef std::vector<MinerNode<PartitionTidList> >::const_iterator MNIter;
} //namespace algos