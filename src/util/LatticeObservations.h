//
// Created by alexandrsmirn
//

#pragma once

#include <unordered_map>
#include "Vertical.h"

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

    bool isCandidate(shared_ptr<Vertical> const& vertical);
    bool isVisited(shared_ptr<Vertical> const& vertical) { return this->find(*vertical) != this->end();} //нужна ли?

    NodeCategory updateDependencyCategory(shared_ptr<Vertical>  const& vertical);
    NodeCategory updateNonDependencyCategory(shared_ptr<Vertical>  const& vertical);

};
