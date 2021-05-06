//
// Created by alexandrsmirn
//

#pragma once

#include <unordered_map>
#include "Vertical.h"
#include "../src/custom/CustomHashes.h"
#include "CustomComparator.h"

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
    using vertical_set = std::unordered_set<shared_ptr<Vertical>, std::hash<shared_ptr<Vertical>>, custom_comparator>;
    //TODO ctors

    bool isCandidate(shared_ptr<Vertical> const& vertical);
    bool isVisited(shared_ptr<Vertical> const& vertical) { return this->find(*vertical) != this->end();} //нужна ли?

    NodeCategory updateDependencyCategory(shared_ptr<Vertical>  const& vertical);
    NodeCategory updateNonDependencyCategory(shared_ptr<Vertical>  const& vertical);

    vertical_set getUncheckedSubsets(shared_ptr<Vertical> const& node);
    vertical_set getUncheckedSupersets(shared_ptr<Vertical> const& node);

    bool inferCategory(shared_ptr<Vertical> const& node);

    bool getCategory(shared_ptr<Vertical> node); //TODO чтобы не возиться с end

};


