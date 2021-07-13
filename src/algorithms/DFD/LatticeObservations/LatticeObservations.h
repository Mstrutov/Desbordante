//
// Created by alexandrsmirn
//

#pragma once

#include <unordered_map>
#include "Vertical.h"
#include "../src/custom/CustomHashes.h"
#include "ColumnOrder.h"

enum class NodeCategory {
        dependency,
        minimalDependency,
        candidateMinimalDependency,
        nonDependency,
        maximalNonDependency,
        candidateMaximalNonDependency
};

class LatticeObservations : public std::unordered_map<Vertical, NodeCategory> {
public:
    //using vertical_set = std::unordered_set<shared_ptr<Vertical>, std::hash<shared_ptr<Vertical>>, custom_comparator>;
    //TODO ctors

    bool isCandidate(Vertical const& vertical);
    bool isVisited(Vertical const& vertical) { return this->find(vertical) != this->end();} //нужна ли?

    NodeCategory updateDependencyCategory(Vertical const& vertical);
    NodeCategory updateNonDependencyCategory(Vertical const& vertical, int rhsIndex);

    std::unordered_set<Vertical> getUncheckedSubsets(const Vertical &node, size_t rhsIndex, ColumnOrder const&);
    std::unordered_set<Vertical> getUncheckedSupersets(const Vertical &node, size_t rhsIndex, ColumnOrder const&);

    bool inferCategory(Vertical const& node);

    bool getCategory(Vertical const& node); //TODO чтобы не возиться с end

};


