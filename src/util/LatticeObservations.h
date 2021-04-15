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
    bool checkIfMinimalDependency(shared_ptr<Vertical> const& vertical);
    bool checkIfMaximalNonDependency(shared_ptr<Vertical>  const& vertical);

};
