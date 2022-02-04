#pragma once

#include <unordered_map>

#include "Vertical.h"
#include "../src/custom/CustomHashes.h"
#include "DFD/ColumnOrder/ColumnOrder.h"
#include "DFD/NodeCategory.h"

class LatticeObservations : public std::unordered_map<Vertical, NodeCategory> {
public:
    bool IsCandidate(Vertical const& node) const;
    bool IsVisited(Vertical const& node) const { return this->find(node) != this->end(); }

    NodeCategory UpdateDependencyCategory(Vertical const& node);
    NodeCategory UpdateNonDependencyCategory(Vertical const& node, unsigned int rhs_index);

    std::unordered_set<Vertical> GetUncheckedSubsets(const Vertical &node, ColumnOrder const&) const;
    std::unordered_set<Vertical> GetUncheckedSupersets(const Vertical &node, unsigned int rhs_index, ColumnOrder const&) const;
};
