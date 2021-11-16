#pragma once

#include <unordered_map>

#include "Vertical.h"
#include "../src/custom/CustomHashes.h"
#include "DFD/ColumnOrder/ColumnOrder.h"
#include "DFD/NodeCategory.h"

class LatticeObservations : public std::unordered_map<Vertical, NodeCategory> {
public:
    bool isCandidate(Vertical const& node) const;
    bool isVisited(Vertical const& node) const { return this->find(node) != this->end(); }

    NodeCategory updateDependencyCategory(Vertical const& node);
    NodeCategory updateNonDependencyCategory(Vertical const& node, unsigned int rhsIndex);

    std::unordered_set<Vertical> getUncheckedSubsets(const Vertical &node, ColumnOrder const&) const;
    std::unordered_set<Vertical> getUncheckedSupersets(const Vertical &node, unsigned int rhsIndex, ColumnOrder const&) const;
};
