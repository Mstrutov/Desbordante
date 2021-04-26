//
// Created by alexandrsmirn
//

#pragma once

#include <unordered_map>
#include "Vertical.h"
#include "../src/custom/CustomHashes.h"

/*namespace std {
    template<>
    struct hash<Vertical> {
        size_t operator()(Vertical const &k) const {
            return k.getColumnIndices().to_ulong();
        }
    };
}*/


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

    //TODO ctors

    bool isCandidate(shared_ptr<Vertical> const& vertical);
    bool isVisited(shared_ptr<Vertical> const& vertical) { return this->find(*vertical) != this->end();} //нужна ли?

    NodeCategory updateDependencyCategory(shared_ptr<Vertical>  const& vertical);
    NodeCategory updateNonDependencyCategory(shared_ptr<Vertical>  const& vertical);

    bool inferCategory(shared_ptr<Vertical> const& node);

    bool getCategory(shared_ptr<Vertical> node); //TODO чтобы не возиться с end

};


