#pragma once

#include <unordered_set>
#include <vector>

#include "order_utility.h"

namespace algos::order {

class ListLattice {
public:
    using LatticeLevel = std::unordered_set<Node, Hash>;
    using PrefixBlocks = std::unordered_map<AttributeList, std::vector<Node>, Hash>;

private:
    LatticeLevel lattice_level_;
    unsigned int level_num_ = 1;

    PrefixBlocks GetPrefixBlocks() const;
    Node JoinNodes(Node const& l, Node const& r) const;

public:
    ListLattice(CandidateSets& candidate_sets, std::vector<AttributeList> single_attributes);

    CandidatePairs ObtainCandidates(Node const& node) const;
    void Prune(CandidateSets& candidate_sets);
    void GenerateNextLevel(CandidateSets& candidate_sets);

    LatticeLevel const& GetLatticeLevel() const {
        return lattice_level_;
    }

    unsigned int GetLevelNumber() const {
        return level_num_;
    }

    bool IsEmpty() const {
        return lattice_level_.empty();
    }
};

}  // namespace algos::order
