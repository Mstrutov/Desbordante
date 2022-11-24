#pragma once

// see ../algorithms/cfd/LICENSE

#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <boost/functional/hash.hpp>

typedef std::vector<int> Transaction;
typedef int Item;
typedef std::vector<Item> Itemset;
typedef std::vector<Item> SimpleTidList;
typedef std::vector<int> Diffset;
typedef boost::hash<std::pair<int, std::string> > pairhash;

struct PartitionTidList {
    PartitionTidList() : sets_number(0) {}
    PartitionTidList(SimpleTidList tids, int nrSets) : tids(std::move(tids)), sets_number(nrSets) {}
    SimpleTidList tids;
    int sets_number;
    static const int SEP;  // = -1;
    bool operator==(const PartitionTidList&) const;
    bool operator!=(const PartitionTidList&) const;
    bool operator<(const PartitionTidList&) const;
};

PartitionTidList Convert(const SimpleTidList&);
SimpleTidList Convert(const PartitionTidList&);
bool LessThan(const PartitionTidList& lhs, const PartitionTidList& rhs);
bool Equals(const PartitionTidList&, const PartitionTidList&);
bool LessThan(const SimpleTidList& lhs, const SimpleTidList& rhs);
