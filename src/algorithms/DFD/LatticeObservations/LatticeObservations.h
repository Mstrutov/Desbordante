#pragma once

#include <unordered_map>

#include "Vertical.h"
#include "../src/custom/CustomHashes.h"
#include "DFD/ColumnOrder/ColumnOrder.h"
#include "DFD/NodeCategory.h"

class LatticeObservations : public std::unordered_map<Vertical, NodeCategory> {
public:
    bool isCandidate(Vertical const& vertical) const;
    bool isVisited(Vertical const& vertical) const { return this->find(vertical) != this->end(); }

    NodeCategory updateDependencyCategory(Vertical const& vertical);
    NodeCategory updateNonDependencyCategory(Vertical const& vertical, int rhsIndex);

    std::unordered_set<Vertical> getUncheckedSubsets(const Vertical &node, ColumnOrder const&) const;
    std::unordered_set<Vertical> getUncheckedSupersets(const Vertical &node, size_t rhsIndex, ColumnOrder const&) const;
};
